import numpy as np
from pyomo.environ import *
import pickle
from rtolib.util.misc import get_hypercube_sampling
from rtolib.core import ProblemDescription

class QuadraticBoostedDCCPWLFunction(object):
    def __init__(self):
        self.dimension = 0
        self.no_pos_segment = 0
        self.no_neg_segment = 0
        self.input_variables = []
        self.parameters = []
        self.no_parameter = 0
        self.pos_cpwl_coeff = np.array([[]])
        self.neg_cpwl_coeff = np.array([[]])
        self.pos_quadratic_A = np.array([])
        self.pos_quadratic_b = np.array([])
        self.pos_quadratic_c = 0
        self.neg_quadratic_A = np.array([]) # should be positive semidefinite
        self.neg_quadratic_b = np.array([])
        self.neg_quadratic_c = 0

    def load_from_cpwl_model_file(self, cpwl_file):
        with open(cpwl_file, "rb") as cpwl_f:
            cpwl_data = pickle.load(cpwl_f)
        self.dimension = len(cpwl_data[0][0])-1
        self.no_pos_segment = len(cpwl_data[0])
        self.no_neg_segment = len(cpwl_data[1])
        self.pos_cpwl_coeff = np.array(cpwl_data[0])
        self.neg_cpwl_coeff = np.array(cpwl_data[1])
        self.pos_quadratic_A = np.eye(self.dimension) * 0
        self.pos_quadratic_b = np.zeros((self.dimension,))
        self.pos_quadratic_c = 0
        self.neg_quadratic_A = np.eye(self.dimension) * 0
        self.neg_quadratic_b = np.zeros((self.dimension,))
        self.neg_quadratic_c = 0

    def load_from_qcpwl_model_file(self, cpwl_file, quadr_file):
        with open(cpwl_file, "rb") as cpwl_f:
            cpwl_data = pickle.load(cpwl_f)
        with open(quadr_file, "rb") as quadr_f:
            quadr_data = pickle.load(quadr_f)
        self.dimension = len(cpwl_data[0][0])-1
        self.no_pos_segment = len(cpwl_data[0])
        self.no_neg_segment = len(cpwl_data[1])
        # print(cpwl_data[2])
        self.pos_cpwl_coeff = np.array(cpwl_data[0])
        self.neg_cpwl_coeff = np.array(cpwl_data[1])
        self.pos_quadratic_A = np.array(quadr_data[0])
        self.pos_quadratic_b = np.array(quadr_data[2])
        self.pos_quadratic_c = np.array(quadr_data[3])
        self.neg_quadratic_A = np.array(quadr_data[1])
        self.neg_quadratic_b = np.zeros((self.dimension,))
        self.neg_quadratic_c = 0

    def simulate(self, input, with_modifier=False, model=None):
        '''

        :param input: RTO inputs and parameters
        :param with_modifier:
        :param model:
        :return:
        '''
        input_vector = np.zeros((self.dimension, ))
        # input dictionary may have more input than the model has
        for i, k in enumerate(self.input_variables):
            input_vector[i] = input[k]
        max_val = -np.Inf
        max_no = None
        for i in range(self.no_pos_segment):
            val = np.dot(self.pos_cpwl_coeff[i,:self.dimension] \
                         , input_vector) + self.pos_cpwl_coeff[i,self.dimension]
            if val > max_val:
                max_no = i
                max_val = val
        pos_cpwl_value = max_val

        pos_quadr_value = np.transpose(input_vector) @ self.pos_quadratic_A @ input_vector + \
                          np.dot(input_vector.reshape((self.dimension, )), self.pos_quadratic_b) + self.pos_quadratic_c

        input_vector = np.zeros((self.dimension, ))
        for i, k in enumerate(self.input_variables):
            input_vector[i] = input[k]
        max_val = -np.Inf
        max_no = None
        for i in range(self.no_neg_segment):
            val = np.dot(self.neg_cpwl_coeff[i, :self.dimension] \
                         , input_vector) + self.neg_cpwl_coeff[i, self.dimension]
            if val > max_val:
                max_no = i
                max_val = val
        neg_cpwl_value = max_val

        neg_quadr_value = np.transpose(input_vector) @ self.neg_quadratic_A @ input_vector + \
                       np.dot(input_vector.reshape((self.dimension, )), self.neg_quadratic_b) + self.neg_quadratic_c

        ret_original_model_value = pos_cpwl_value+pos_quadr_value-neg_cpwl_value-neg_quadr_value
        if not with_modifier:
            return ret_original_model_value
        else:
            linear_correction_k, linear_correction_b = self.get_linear_correction(model)
            base_point = self.get_base_point(model)
            return ret_original_model_value+np.dot(input_vector-base_point, linear_correction_k)+linear_correction_b

    def build_optimizer_mainbody(self, model, input):
        '''

        :param input: the input variable of the parent model
        :param model: assumes to be a block
        :param output_name:
        :return:
        '''
        model.PosSeg = Set(initialize=range(self.no_pos_segment))
        model.Dimension = Set(initialize=range(self.dimension))

        model.cpwl_pos = Var()
        model.cpwl_neg = Var()
        # model.quadr_pos = Var()
        model.quadr_neg = Var()
        model.linear_correction = Var()
        model.linear_correction_b = Param(initialize=0, mutable=True)
        model.base_point = Param(model.Dimension, mutable=True, initialize=0)
        model.cpwl_neg_k = Param(model.Dimension, mutable=True)
        model.cpwl_neg_b = Param(mutable=True)
        model.quadr_neg_k = Param(model.Dimension, mutable=True)
        model.quadr_neg_b = Param(mutable=True)
        def quadr_pos_calc(m):
            return sum([self.pos_quadratic_A[i,j]*input[self.input_variables[i]]*\
                                       input[self.input_variables[j]] for j in m.Dimension for i in m.Dimension])+\
                           sum([self.pos_quadratic_b[i]* input[self.input_variables[i]] for i in m.Dimension])+\
                  self.pos_quadratic_c
        model.quadr_pos = Expression(rule=quadr_pos_calc)
        def cpwl_pos_calc(m, s):
            return m.cpwl_pos >= sum([self.pos_cpwl_coeff[s,i]* \
                                     input[self.input_variables[i]] for i in m.Dimension])+ \
                   self.pos_cpwl_coeff[s,self.dimension]
        model.cpwl_pos_calc=Constraint(model.PosSeg, rule=cpwl_pos_calc)
        def cpwl_neg_calc(m):
            return m.cpwl_neg == sum([m.cpwl_neg_k[i]* \
                                     input[self.input_variables[i]] for i in m.Dimension])+ \
                                    m.cpwl_neg_b
        model.cpwl_neg_calc = Constraint(rule=cpwl_neg_calc)
        def quadr_neg_calc(m):
            return m.quadr_neg == sum([m.quadr_neg_k[i] * \
                                      input[self.input_variables[i]] for i in m.Dimension]) + \
                   m.quadr_neg_b
        model.quadr_neg_calc = Constraint(rule=quadr_neg_calc)
        def f_con(m):
            return model.cpwl_pos - model.cpwl_neg + model.quadr_pos - \
                   model.quadr_neg + model.linear_correction
        model.f=Expression(rule=f_con)


    def build_optimizer_linear_correction(self, model, input):
        model.linear_correction_k = Param(model.Dimension, initialize=0, mutable=True)
        def linear_correction_calc(m):
            return m.linear_correction == sum([m.linear_correction_k[i] * \
                                               (input[self.input_variables[i]] - m.base_point[i]) for i in
                                               m.Dimension]) + \
                   m.linear_correction_b
        model.linear_correction_calc = Constraint(rule=linear_correction_calc)

    def get_linear_correction(self, model):
        linear_correction_k = np.zeros((len(model.Dimension),))
        linear_correction_b = 0
        for i in model.Dimension:
            linear_correction_k[i] = value(model.linear_correction_k[i])
        linear_correction_b = value(model.linear_correction_b)
        return linear_correction_k, linear_correction_b

    def get_base_point(self, model):
        base_point = np.zeros((len(model.Dimension),))
        for i in model.Dimension:
            base_point[i] = value(model.base_point[i])
        return base_point

    def set_optimizer_concave_majorization(self, model, input):
        input_vector = np.zeros((self.dimension,1))
        for i,k in enumerate(self.input_variables):
            input_vector[i,0] = input[k]
        max_val = -np.Inf
        max_no = None
        for i in range(self.no_neg_segment):
            val = np.dot(self.neg_cpwl_coeff[i,:self.dimension]\
                         , input_vector) + self.neg_cpwl_coeff[i,self.dimension]
            if val > max_val:
                max_no = i
                max_val = val
        for i in range(self.dimension):
            model.cpwl_neg_k[i] = self.neg_cpwl_coeff[max_no,i]
        model.cpwl_neg_b = self.neg_cpwl_coeff[max_no,self.dimension]
        
        quadr_grad = 2*(self.neg_quadratic_A @ input_vector).reshape((self.dimension,))+self.neg_quadratic_b
        val_at_input = np.transpose(input_vector)@self.neg_quadratic_A @ input_vector+\
                    np.dot(input_vector.reshape((self.dimension,)), self.neg_quadratic_b) + self.neg_quadratic_c
        quadr_majorization_b = val_at_input[0,0] - np.dot(quadr_grad,input_vector.reshape((self.dimension,)))
        for i in range(self.dimension):
            model.quadr_neg_k[i] = quadr_grad[i]
        model.quadr_neg_b = quadr_majorization_b

    def set_linear_correction(self, model, modifier_epsilon, modifier_lambda, base_point):
        modifier_lambda_vector = np.zeros((self.dimension, ))
        base_point_vector = np.zeros((self.dimension, ))
        for i, k in enumerate(self.input_variables):
            modifier_lambda_vector[i] = modifier_lambda[k]
            base_point_vector[i] = base_point[k]
        for i in range(self.dimension):
            model.linear_correction_k[i] = modifier_lambda_vector[i]
            model.base_point[i] = base_point_vector[i]
        model.linear_correction_b = modifier_epsilon
        
    def build_parameter_estimator(self, model, alpha, n_samples):
        '''

        :param model:
        :param alpha: list of pyomo obj, must be consistent with parameters stated
        in QuadraticBoostedDCCPWLFunction initialization
        :param n_samples:
        :return:
        '''
        model.Dimension = Set(initialize=range(self.dimension))
        model.Samples = Set(initialize=range(n_samples))
        model.PosSeg = Set(initialize=range(self.no_pos_segment))
        model.NegSeg = Set(initialize=range(self.no_neg_segment))
        model.Parameters = Set(initialize=range(self.no_parameter))
        # bounded alpha and beta needed for convergence
        model.alpha = Var(model.Parameters, initialize=0)  # , bounds=(-100,100)
        model.error_pos = Var(model.Samples)
        model.error_neg = Var(model.Samples)
        model.beta = Param(model.Parameters, mutable=True)
        model.ineq1_u_diff_u = Param(model.Samples, model.PosSeg, model.Parameters, mutable=True)
        model.ineq2_u_diff_v = Param(model.Samples, model.PosSeg, model.Parameters, mutable=True)
        model.ineq3_v_diff_u = Param(model.Samples, model.NegSeg, model.Parameters, mutable=True)
        model.ineq4_v_diff_v = Param(model.Samples, model.NegSeg, model.Parameters, mutable=True)
        def ineq1(m, sample, pos_seg):
            return m.error_pos[sample] >= \
                   sum([m.alpha[p] * m.ineq1_u_diff_u[sample, pos_seg, p] for p in m.Parameters])
        model.ineq1 = Constraint(model.Samples, model.PosSeg, rule=ineq1)
        def ineq2(m, sample, pos_seg):
            return m.error_pos[sample] >= \
                   sum([m.alpha[p] * m.ineq2_u_diff_v[sample, pos_seg, p] for p in m.Parameters]) - \
                   f_samples[sample]
        model.ineq2 = Constraint(model.Samples, model.PosSeg, rule=ineq2)
        def ineq3(m, sample, neg_seg):
            return m.error_neg[sample] >= \
                   sum([m.alpha[p] * m.ineq3_v_diff_u[sample, neg_seg, p] for p in m.Parameters]) + \
                   f_samples[sample]
        model.ineq3 = Constraint(model.Samples, model.NegSeg, rule=ineq3)
        def ineq4(m, sample, neg_seg):
            return m.error_neg[sample] >= \
                   sum([m.alpha[p] * m.ineq4_v_diff_v[sample, neg_seg, p] for p in m.Parameters])
        model.ineq4 = Constraint(model.Samples, model.NegSeg, rule=ineq4)
        def obj(m):
            return 2 * sum([m.error_pos[sample] ** 2 for sample in model.Samples]) + \
                   2 * sum([m.error_neg[sample] ** 2 for sample in model.Samples]) - \
                   sum([m.beta[p] * m.alpha[p] for p in m.Parameters])
        model.obj = Expression(rule=obj)

    def update_parameter_estimation_model(self, model, alpha_value, x_samples, f_samples):
        '''
        Set beta, ineq1_u_diff_u, ineq2_u_diff_v, ineq3_v_diff_u, ineq4_v_diff_v
        :param model:
        :param alpha_value: list
        :param x_samples: ndarray n_sample*n_dim
        :param f_samples: ndarray n_sample
        :return:
        '''
        n_pos = self.no_pos_segment
        n_neg = self.no_neg_segment
        n_dim = self.dimension
        n_samples = x_samples.shape[0]

        alpha = np.zeros(shape=(self.no_parameter,))
        for i in range(n_para):
            alpha[i] = value(alpha_value[i])
        u = np.zeros(shape=(n_samples, n_pos, n_para))
        v = np.zeros(shape=(n_samples, n_neg, n_para))
        f = np.zeros(shape=(n_samples))
        g = np.zeros(shape=(n_samples))
        h = np.zeros(shape=(n_samples))
        error = np.zeros(shape=(n_samples))  # p_i(\alpha)
        j = np.zeros(shape=(n_samples), dtype=int)
        q = np.zeros(shape=(n_samples), dtype=int)
        for s in range(n_samples):
            for sg in range(n_pos):
                u[s, sg, (sg * (n_dim + 1)):((sg + 1) * (n_dim + 1) - 1)] = x_samples[s, :]
                u[s, sg, (sg + 1) * (n_dim + 1) - 1] = 1
        for s in range(n_samples):
            for sg in range(n_neg):
                r = sg + n_pos
                v[s, sg, (r * (n_dim + 1)):((r + 1) * (n_dim + 1) - 1)] = x_samples[s, :]
                v[s, sg, (r + 1) * (n_dim + 1) - 1] = 1
        for s in range(n_samples):
            maximum = -np.Inf
            maximum_index = None
            for sg in range(n_pos):
                val = np.dot(alpha, u[s, sg, :].reshape((n_para,)))
                if val > maximum:
                    maximum = val
                    maximum_index = sg
            g[s] = maximum
            j[s] = maximum_index
        for s in range(n_samples):
            maximum = -np.Inf
            maximum_index = None
            for sg in range(n_neg):
                val = np.dot(alpha, v[s, sg, :].reshape((n_para,)))
                if val > maximum:
                    maximum = val
                    maximum_index = sg
            h[s] = maximum
            q[s] = maximum_index
        for s in range(n_samples):
            f[s] = g[s] - h[s]
            error[s] = f[s] - f_samples[s]
        for p in model.Parameters:
            model.beta[p] = 2 * sum([error[s] * (u[s, j[s], p] - v[s, q[s], p]) for s in model.Samples])
        # model.beta.pprint()
        for s in model.Samples:
            for sg in model.PosSeg:
                for p in model.Parameters:
                    model.ineq1_u_diff_u[s, sg, p] = u[s, sg, p] - u[s, j[s], p]
                    model.ineq2_u_diff_v[s, sg, p] = u[s, sg, p] - v[s, q[s], p]
        for s in model.Samples:
            for sg in model.NegSeg:
                for p in model.Parameters:
                    model.ineq3_v_diff_u[s, sg, p] = v[s, sg, p] - u[s, j[s], p]
                    model.ineq4_v_diff_v[s, sg, p] = v[s, sg, p] - v[s, q[s], p]
        return f, error


class QuadraticBoostedDCCPWLFunctionSubgrad(QuadraticBoostedDCCPWLFunction):
    def set_linear_correction_parameters(self, model, plant_y, plant_k, base_point, merit_gd_direction):
        # dictionary to ndarray
        plant_k_vector = np.zeros((self.dimension,))
        base_point_vector = np.zeros((self.dimension, ))
        merit_gd_direction_vector = np.zeros((self.dimension, ))
        for i, k in enumerate(self.input_variables):
            plant_k_vector[i] = plant_k[k]
            base_point_vector[i] = base_point[k]
            merit_gd_direction_vector[i] = merit_gd_direction[k]
        merit_gd_direction_vector/=np.linalg.norm(merit_gd_direction_vector)

        # set base point, plant_k and linear_correction_b
        self.plant_k = plant_k_vector
        self.base_point = base_point_vector
        self.linear_correction_b = plant_y - self.simulate(base_point) # correct

        # set subgradient of concave part
        input_vector = base_point_vector.reshape((self.dimension, 1))
        max_val = -np.Inf
        max_no = None
        for i in range(self.no_neg_segment):
            val = np.dot(self.neg_cpwl_coeff[i, :self.dimension] \
                         , input_vector) + self.neg_cpwl_coeff[i, self.dimension]
            if val > max_val:
                max_no = i
                max_val = val
        neg_cpwl_grad = self.neg_cpwl_coeff[max_no, :self.dimension].reshape((self.dimension,))
        neg_quadr_grad = 2*self.neg_quadratic_A @ base_point_vector + self.neg_quadratic_b
        neg_grad = neg_cpwl_grad+neg_quadr_grad
        self.concave_k = neg_grad

        # get the current active convex segments
        merit_gd_point = base_point_vector + merit_gd_direction_vector* 1e-3
        merit_gd_point = merit_gd_point.reshape((self.dimension, 1))
        max_val = -np.Inf
        max_no = None
        for i in range(self.no_pos_segment):
            val = np.dot(self.pos_cpwl_coeff[i, :self.dimension] \
                         , merit_gd_point) + self.pos_cpwl_coeff[i, self.dimension]
            if val > max_val:
                max_no = i
                max_val = val
        pos_cpwl_grad = self.pos_cpwl_coeff[max_no, :self.dimension].reshape((self.dimension,))
        pos_quadr_grad = 2 * self.pos_quadratic_A @ base_point_vector + self.pos_quadratic_b
        pos_grad = pos_cpwl_grad + pos_quadr_grad
        self.convex_k = pos_grad

        # set base point and modifier
        model.linear_correction_b = self.linear_correction_b
        for d in range(self.dimension):
            model.linear_correction_k[d] = plant_k_vector[d]-(self.convex_k-self.concave_k)[d]
            model.base_point[d] = base_point_vector[d]


class QuadraticBoostedDCCPWL_RTOObject(object):
    def __init__(self, dc_cpwl_functions, inputs, cvs, parameters):
        '''

        :param dc_cpwl_functions:
        :param inputs:  tuple or list, subset of all inputs, including parameters
        :param cvs: tuple or list
        :param parameters: tuple or list
        '''
        self.dc_cpwl_functions = dc_cpwl_functions
        self.inputs = inputs
        self.cvs = cvs
        self.parameters = parameters
        self.tol = 1e-5
        self.subproblem_max_iter = 100

    def set_input_bounds(self, input_bounds):
        self.input_bounds = input_bounds

    def build(self, problem_description, constraint_scaling=1):
        assert isinstance(problem_description, ProblemDescription)
        self.problem_description = problem_description
        self.problem_description.scaling_factors['validity_con'] = 1

        subproblem_model = ConcreteModel()
        subproblem_model.InputIndex = Set(initialize=self.inputs)
        subproblem_model.OutputIndex = Set(initialize=self.cvs)
        subproblem_model.input = Var(subproblem_model.InputIndex)
        subproblem_model.output = Block(subproblem_model.OutputIndex)
        subproblem_model.slack = Var(subproblem_model.OutputIndex, within=NonNegativeReals)
        for name, dc_cpwl in self.dc_cpwl_functions.items():
            dc_cpwl.build_optimizer_mainbody(subproblem_model.output[name], subproblem_model.input)
            dc_cpwl.build_optimizer_linear_correction(subproblem_model.output[name], subproblem_model.input)

        def optimization_constraint(m):
            for con_name in self.cvs[1:]:
                yield m.output[con_name].f <= 0
        subproblem_model.optimization_constraint = ConstraintList(rule=optimization_constraint)

        def optimization_obj(m):
            return m.output[self.cvs[0]].f
        subproblem_model.optimization_obj = Objective(rule=optimization_obj, sense=minimize)

        self.subproblem_model = subproblem_model

    def build_parameter_estimation(self, n_samples, scaling_factor):
        self.no_pe_samples = n_samples
        pe_model = ConcreteModel()
        pe_model.ParameterIndex = Set(initialize=self.parameters)
        pe_model.OutputIndex = Set(initialize=self.cvs)
        pe_model.parameter = Var(pe_model.ParameterIndex)
        pe_model.output = Block(pe_model.OutputIndex)
        for name, dc_cpwl in self.dc_cpwl_functions.items():
            dc_cpwl.build_parameter_estimator(pe_model.output[name], pe_model.parameter, n_samples)
        def obj(m):
            return sum([m.output[name].submodel_obj/scaling_factor[name]**2 \
                        for name in pe_model.OutputIndex])
        pe_model.obj = Objective(rule=obj, sense=minimize)
        self.pe_model = pe_model

    def estimate_parameter(self, plant_data):
        '''

        :param plant_data: list of dict, every dictionary includes both input variable value
        and output variable value
        :param scaling_factor:
        :return:
        '''
        last_rmse = np.Inf
        tol = 1e-4
        # initializing parameter
        for para_name in self.parameters:
            self.DC_CPWL_RTO_model.pe_model[para_name] = parameter_value[para_name]
        alpha_value, x_samples, f_samples

        # parameter estimation
        # TODO: update for each output
        f_fit, f_error = self.update_parameter_estimation_model(self.pe_model, \
                        alpha_value, x_samples, f_samples)
        results = self.solver.solve(self.pe_model, tee=self.tee)
        if not ((results.solver.status == SolverStatus.ok) and (
                results.solver.termination_condition == TerminationCondition.optimal)):
            print(results.solver.termination_condition)
            raise Exception("QCQP solving fails.")
        rmse = np.linalg.norm(f_error) / self.no_pe_samples

        iter_count = 0
        while abs(rmse - last_rmse) >= tol or iter_count < self.subproblem_max_iter:
            iter_count += 1
            last_rmse = rmse
            f_fit, f_error = self.update_parameter_estimation_model(self.pe_model, \
                        alpha_value, x_samples, f_samples)
            results = self.solver.solve(self.pe_model, tee=self.tee)
            if not ((results.solver.status == SolverStatus.ok) and (
                    results.solver.termination_condition == TerminationCondition.optimal)):
                print(results.solver.termination_condition)
                raise Exception("QCQP solving fails.")
            rmse = np.linalg.norm(f_error) / self.no_pe_samples
        return f_error


    def update_modifiers(self, modifiers, base_point):
        # plant_y_and_k and base_point may contain more input
        # and output than the model has

        # print(modifiers)
        for cv in self.cvs:
            if cv == "validity_con":
                continue
            output_model = self.subproblem_model.output[cv]
            modifier_epsilon = modifiers[(cv, None)]
            modifier_lambda = {}
            for mv in self.inputs:
                if (cv, mv) in modifiers.keys():
                    modifier_lambda[mv] = modifiers[(cv, mv)]
                else:
                    modifier_lambda[mv] = 0
            self.dc_cpwl_functions[cv].set_linear_correction(output_model, \
                                modifier_epsilon, modifier_lambda, base_point)

    def get_input_vector(self):
        input_vector = np.zeros((len(self.inputs),))
        for i,mv in enumerate(self.inputs):
            input_vector[i] = value(self.subproblem_model.input[mv])
        return input_vector

    def calculate_infeasibility(self, point):
        y,_ = self.simulate(point, with_modifier=True)
        infeasibility = 0
        for i,cv in enumerate(self.cvs):
            if i == 0:
                continue
            if y[cv] > -1e-3:
                infeasibility += (y[cv]+1e-3)**2
        return infeasibility

    def optimize(self, spec_values, starting_point):
        '''

        :param spec_values: specification and parameters
        :param starting_point: other inputs
        :return:
        '''
        print(spec_values)
        print(starting_point)
        whole_input = {}
        for k,v in spec_values.items():
            if k not in self.inputs:
                continue
            if k not in self.parameters:
                whole_input[k] = v
                self.subproblem_model.input[k].fix(v)
        for k,v in starting_point.items():
            if k not in self.inputs:
                continue
            whole_input[k] = v
            self.subproblem_model.input[k] = v

            # -------------- find a feasible point ---------------------
            variable_input_names = []
            for k in starting_point.keys():
                if k not in spec_values.keys():
                    if k in self.inputs:
                        variable_input_names.append(k)
            variable_input_dimensionality = len(variable_input_names)
            n_data = 100
            lb = np.zeros((variable_input_dimensionality,))
            ub = np.zeros((variable_input_dimensionality,))
            current_point = np.zeros((variable_input_dimensionality,))
            for i, mv in enumerate(variable_input_names):
                lb[i] = self.input_bounds[mv][0]
                ub[i] = self.input_bounds[mv][1]
                current_point[i] = starting_point[mv]
            candidate_input1 = get_hypercube_sampling(variable_input_dimensionality, n_data, lb, ub, seed=1)
            full_input_dimensionality = len(self.inputs)
            candidate_input = np.zeros((n_data, full_input_dimensionality))
            for i in range(n_data):
                for j, mv in enumerate(self.inputs):
                    if mv not in variable_input_names:
                        candidate_input[i, j] = spec_values[mv]
                    else:
                        candidate_input[i, j] = candidate_input1[i, variable_input_names.index(mv)]
            infeasibility = np.zeros((n_data,))
            for i in range(n_data):
                candidate_input_dict = {}
                for j, k in enumerate(self.inputs):
                    candidate_input_dict[k] = candidate_input[i, j]
                infeasibility[i] = self.calculate_infeasibility(candidate_input_dict)
            nearest_dist = np.inf
            nearest_point_no = None
            for i in range(n_data):
                if infeasibility[i] <= 1e-4:
                    dist = np.linalg.norm(candidate_input1[i, :] - current_point)
                    if dist < nearest_dist:
                        nearest_dist = dist
                        nearest_point_no = i
            for i, k in enumerate(self.inputs):
                whole_input[k] = candidate_input[nearest_point_no, i]
                self.subproblem_model.input[k] = whole_input[k]
            print("find a feasible point to start optimization")
            print(whole_input)

            # -------------- solve optimization problem ----------------
            prev_u = self.get_input_vector()
            # print(prev_u)
            solve_status = False
            for i in range(self.subproblem_max_iter):
                for name, dc_cpwl in self.dc_cpwl_functions.items():
                    dc_cpwl.set_optimizer_concave_majorization(self.subproblem_model.output[name], whole_input)
                with open('before-optimize.txt', 'w') as model_file:
                    self.subproblem_model.pprint(model_file)
                results = self.solver.solve(self.subproblem_model, tee=self.tee)  # self.tee
                if not ((results.solver.status == SolverStatus.ok) and (
                        results.solver.termination_condition == TerminationCondition.optimal)):
                    print(results.solver.termination_condition)
                    raise Exception("QCQP solving fails.")
                for mv in self.inputs:
                    whole_input[mv] = value(self.subproblem_model.input[mv])
                this_u = self.get_input_vector()
                if np.linalg.norm(this_u - prev_u) < self.tol:
                    solve_status = True
                    break
                prev_u = this_u
            print("after optimization:")
            print(whole_input)
            return whole_input, solve_status

    def simulate(self, point, with_modifier=False):
        ret = {}
        for cv in self.cvs:
            ret[cv]=self.dc_cpwl_functions[cv].simulate(point, with_modifier, self.subproblem_model.output[cv])
        solve_status = True
        return ret, solve_status

    def set_model_parameters(self, parameter_value):
        '''

        :param parameter_value: dict
        :return:
        '''
        for para_name in self.parameters:
            self.DC_CPWL_RTO_model.pe_model[para_name] = parameter_value[para_name]

    def set_solver(self, solver, tee, default_options):
        self.solver = solver
        self.tee = tee
        for k, v in default_options.items():
            self.solver.options[k] = v

class QuadraticBoostedDCCPWL_RTOObjectSubgrad(QuadraticBoostedDCCPWL_RTOObject):
    def build(self, problem_description, constraint_scaling=1):
        assert isinstance(problem_description, ProblemDescription)
        self.problem_description = problem_description
        self.problem_description.scaling_factors['validity_con'] = 1

        subproblem_model = ConcreteModel()
        subproblem_model.InputIndex = Set(initialize=self.inputs)
        subproblem_model.OutputIndex = Set(initialize=self.cvs)
        subproblem_model.input = Var(subproblem_model.InputIndex)
        subproblem_model.output = Block(subproblem_model.OutputIndex)
        subproblem_model.slack = Var(subproblem_model.OutputIndex, within=NonNegativeReals)
        setattr(subproblem_model, "infeasibility_sum", Var(initialize=0, within=NonNegativeReals))
        setattr(subproblem_model, "tr_penalty_coeff", Var(initialize=0))
        getattr(subproblem_model, "tr_penalty_coeff").fixed = True
        for name, dc_cpwl in self.dc_cpwl_functions.items():
            dc_cpwl.build_optimizer_mainbody(subproblem_model.output[name], subproblem_model.input)
            dc_cpwl.build_optimizer_linear_correction(subproblem_model.output[name], subproblem_model.input)

        # inequality violation measure
        for con_var_name in self.cvs[1:]:
            setattr(subproblem_model, "con_vio_" + con_var_name, Var(initialize=0, within=NonNegativeReals))

        def _violation_cons(m):
            for con_var_name in self.cvs[1:]:
                yield (m.output[con_var_name].f /
                       self.problem_description.scaling_factors[con_var_name] \
                       - getattr(subproblem_model, "con_vio_" + con_var_name)) * constraint_scaling <= 0
        subproblem_model.violation_cons = ConstraintList(rule=_violation_cons)

        # merit function
        def _tr_infeasibility_sum(m):
            ret = 0
            for con_var_name in self.cvs[1:]:
                ret += getattr(m, "con_vio_" + con_var_name) / \
                       self.problem_description.scaling_factors[con_var_name]
            return ret - m.infeasibility_sum == 0
        subproblem_model.tr_infeasibility_sum = Constraint(rule=_tr_infeasibility_sum)

        # def _tr_merit_function(m):
        #     ret = m.output[self.cvs[0]].f
        #     ret += m.infeasibility_sum * m.tr_penalty_coeff
        #     return ret
        # subproblem_model.tr_merit_function = Expression(rule=_tr_merit_function)
        #
        # def optimization_obj(m):
        #     return m.tr_merit_function
        # subproblem_model.optimization_obj = Objective(rule=optimization_obj, sense=minimize)

        def _ineq_con(m):
            return m.infeasibility_sum == 0
        subproblem_model.ineq_con = Constraint(rule=_ineq_con)

        def optimization_obj(m):
            return m.output[self.cvs[0]].f
        subproblem_model.optimization_obj = Objective(rule=optimization_obj, sense=minimize)

        self.subproblem_model = subproblem_model

    def set_penalty_coeff(self, sigma):
        self.subproblem_model.tr_penalty_coeff.fix(sigma)

    def update_modifiers(self, plant_y_and_k, base_point, sigma):
        # plant_y_and_k and base_point may contain more input
        # and output than the model has
        merit_gd_direction = {}
        for mv in self.inputs:
            if (self.cvs[0], mv) in plant_y_and_k.keys():
                merit_gd_direction[mv] = -plant_y_and_k[(self.cvs[0], mv)]\
                                         /self.problem_description.scaling_factors[self.cvs[0]]
            else:
                merit_gd_direction[mv] = 0
        for cv in self.cvs[1:]:
            if cv == "validity_con":
                continue
            for mv in self.inputs:
                if (cv, mv) in plant_y_and_k.keys():
                    if plant_y_and_k[(cv, None)] >= 0:
                        merit_gd_direction[mv] += -sigma*plant_y_and_k[(cv, mv)]\
                                         /self.problem_description.scaling_factors[cv]

        for cv in self.cvs:
            if cv == "validity_con":
                continue
            plant_y = plant_y_and_k[(cv, None)]
            plant_k = {}
            for mv in self.inputs:
                if (cv, mv) in plant_y_and_k.keys():
                    plant_k[mv] = plant_y_and_k[(cv, mv)]
                else:
                    plant_k[mv] = 0
            self.dc_cpwl_functions[cv].set_linear_correction_parameters(self.subproblem_model.output[cv],\
                                plant_y, plant_k, base_point, merit_gd_direction)

    def optimize(self, spec_values, starting_point):
        '''

        :param spec_values: specification and parameters
        :param starting_point: other inputs
        :return:
        '''
        print(spec_values)
        print(starting_point)
        whole_input = {}
        for k,v in spec_values.items():
            if k not in self.inputs:
                continue
            if k not in self.parameters:
                whole_input[k] = v
                self.subproblem_model.input[k].fix(v)
        for k,v in starting_point.items():
            if k not in self.inputs:
                continue
            whole_input[k] = v
            self.subproblem_model.input[k] = v

        prev_u = self.get_input_vector()
        # print(prev_u)
        solve_status = False
        for i in range(self.subproblem_max_iter):
            for name, dc_cpwl in self.dc_cpwl_functions.items():
                dc_cpwl.set_optimizer_concave_majorization(self.subproblem_model.output[name], whole_input)
            # with open('before-optimize.txt', 'w') as model_file:
            #     self.subproblem_model.pprint(model_file)
            results = self.solver.solve(self.subproblem_model, tee=self.tee) #self.tee
            if not ((results.solver.status == SolverStatus.ok) and (
                    results.solver.termination_condition == TerminationCondition.optimal)):
                print(results.solver.termination_condition)
                raise Exception("QCQP solving fails.")
            for mv in self.inputs:
                whole_input[mv] = value(self.subproblem_model.input[mv])
            this_u = self.get_input_vector()
            if np.linalg.norm(this_u - prev_u) < self.tol:
                solve_status = True
                break
            prev_u = this_u
        print("after optimization:")
        print(whole_input)
        return whole_input, solve_status


class QuadraticBoostedDCCPWL_PenaltyTR_Object(QuadraticBoostedDCCPWL_RTOObjectSubgrad):
    def build(self, problem_description, constraint_scaling=1):
        assert isinstance(problem_description, ProblemDescription)
        self.problem_description = problem_description
        self.problem_description.scaling_factors['validity_con']=1

        subproblem_model = ConcreteModel()
        subproblem_model.InputIndex = Set(initialize=self.inputs)
        subproblem_model.OutputIndex = Set(initialize=self.cvs)
        subproblem_model.input = Var(subproblem_model.InputIndex)
        subproblem_model.output = Block(subproblem_model.OutputIndex)
        setattr(subproblem_model, "infeasibility_sum", Var(initialize=0, within=NonNegativeReals))
        setattr(subproblem_model, "tr_penalty_coeff", Var(initialize=0))
        getattr(subproblem_model, "tr_penalty_coeff").fixed = True
        for name, dc_cpwl in self.dc_cpwl_functions.items():
            dc_cpwl.build_optimizer_mainbody(subproblem_model.output[name], subproblem_model.input)
            dc_cpwl.build_optimizer_linear_correction(subproblem_model.output[name], subproblem_model.input)

        # inequality violation measure
        for con_var_name in self.cvs[1:]:
            setattr(subproblem_model, "con_vio_" + con_var_name, Var(initialize=0, within=NonNegativeReals))

        def _violation_cons(m):
            for con_var_name in self.cvs[1:]:
                yield (m.output[con_var_name].f /
                       self.problem_description.scaling_factors[con_var_name] \
                       - getattr(subproblem_model, "con_vio_" + con_var_name)) * constraint_scaling <= 0
        subproblem_model.violation_cons = ConstraintList(rule=_violation_cons)

        # merit function
        def _tr_infeasibility_sum(m):
            ret = 0
            for con_var_name in self.cvs[1:]:
                ret += getattr(m, "con_vio_" + con_var_name)  / \
                       self.problem_description.scaling_factors[con_var_name]
            return ret - m.infeasibility_sum == 0
        subproblem_model.tr_infeasibility_sum = Constraint(rule=_tr_infeasibility_sum)

        def _tr_merit_function(m):
            ret = m.output[self.cvs[0]].f
            ret += m.infeasibility_sum * m.tr_penalty_coeff
            return ret
        subproblem_model.tr_merit_function = Expression(rule=_tr_merit_function)

        # set tr con
        for ip in subproblem_model.InputIndex:
            setattr(subproblem_model, ip + "_tr_base", Var(initialize=0))
            getattr(subproblem_model, ip + "_tr_base").fixed = True
        setattr(subproblem_model, "tr_radius", Var(initialize=0))
        getattr(subproblem_model, "tr_radius").fixed = True

        def _trust_region_con(m):
            r = 0
            for ip in m.InputIndex:
                r += (m.input[ip] - getattr(m, ip + "_tr_base")) ** 2 / \
                     (self.problem_description.scaling_factors[ip] ** 2)
            return r <= m.tr_radius ** 2
        subproblem_model.trust_region_cons = Constraint(rule=_trust_region_con)

        # set obj
        def _objective(m):
            return m.tr_merit_function
        subproblem_model.objective_function = Objective(rule=_objective, sense=minimize)
        self.subproblem_model = subproblem_model

    def optimize(self, spec_values, starting_point, tr_radius, tr_base):
        '''

        :param spec_values: specification and parameters
        :param starting_point: other inputs
        :return:
        '''
        # print(spec_values)
        # print(starting_point)
        whole_input = {}
        for k,v in spec_values.items():
            if k not in self.inputs:
                continue
            if k not in self.parameters:
                whole_input[k] = v
                self.subproblem_model.input[k].fix(v)
        for k,v in starting_point.items():
            if k not in self.inputs:
                continue
            whole_input[k] = v
            self.subproblem_model.input[k] = v

        # set trust-region
        for var_name in self.subproblem_model.InputIndex:
            if var_name in tr_base.keys():
                getattr(self.subproblem_model, var_name + "_tr_base").fix(tr_base[var_name])
            else:
                getattr(self.subproblem_model, var_name + "_tr_base").fixed = False
        getattr(self.subproblem_model, "tr_radius").fix(tr_radius)

        # -------------- solve optimization problem ----------------
        prev_u = self.get_input_vector()
        # print(prev_u)
        solve_status = False
        for i in range(self.subproblem_max_iter):
            for name, dc_cpwl in self.dc_cpwl_functions.items():
                dc_cpwl.set_optimizer_concave_majorization(self.subproblem_model.output[name], whole_input)
            # self.subproblem_model.display("temp.txt")
            # with open('before-optimize.txt', 'w') as model_file:
            #     self.subproblem_model.pprint(model_file)
            results = self.solver.solve(self.subproblem_model, tee=self.tee)#
            if not ((results.solver.status == SolverStatus.ok) and (
                    results.solver.termination_condition == TerminationCondition.optimal)):
                print(results.solver.termination_condition)
                raise Exception("QCQP solving fails.")
            for mv in self.inputs:
                whole_input[mv] = value(self.subproblem_model.input[mv])
            this_u = self.get_input_vector()
            if np.linalg.norm(this_u - prev_u) < self.tol:
                solve_status = True
                break
            prev_u = this_u
        print("after optimization:")
        print(whole_input)
        return whole_input, solve_status