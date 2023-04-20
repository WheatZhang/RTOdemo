import numpy as np
from pyomo.environ import *
import pickle

class QuadraticBoostedDCCPWLFunction(object):
    def __init__(self):
        self.dimension = 0
        self.no_pos_segment = 0
        self.no_neg_segment = 0
        self.input_variables = []
        self.parameters = []
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
        print(cpwl_data[2])
        self.pos_cpwl_coeff = np.array(cpwl_data[0])
        self.neg_cpwl_coeff = np.array(cpwl_data[1])
        self.pos_quadratic_A = np.array(quadr_data[0])
        self.pos_quadratic_b = np.array(quadr_data[2])
        self.pos_quadratic_c = np.array(quadr_data[3])
        self.neg_quadratic_A = np.array(quadr_data[1])
        self.neg_quadratic_b = np.zeros((self.dimension,))
        self.neg_quadratic_c = 0

    def simulate(self, input, with_modifier=False, model=None):
        input_vector = np.zeros((self.dimension, 1))
        for i, k in enumerate(self.input_variables):
            input_vector[i, 0] = input[k]
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

        input_vector = np.zeros((self.dimension, 1))
        for i, k in enumerate(self.input_variables):
            input_vector[i, 0] = input[k]
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

        ret_original_model_value = pos_cpwl_value[0]+pos_quadr_value[0,0]-neg_cpwl_value[0]-neg_quadr_value[0,0]
        if not with_modifier:
            return ret_original_model_value
        else:
            linear_correction_k, linear_correction_b = self.get_linear_correction(model)
            return ret_original_model_value+np.dot(input_vector, linear_correction_k)+linear_correction_b

    def build_optimizer(self, model, input):
        '''

        :param input: the input variable of the parent model
        :param model: assumes to be a block
        :param output_name:
        :return:
        '''
        model.PosSeg = Set(initialize=range(self.no_pos_segment))
        model.Dimension = Set(initialize=range(self.dimension))

        model.f = Var()
        model.cpwl_pos = Var()
        model.cpwl_neg = Var()
        model.quadr_pos = Var()
        model.quadr_neg = Var()
        model.linear_correction = Var()
        model.linear_correction_k = Param(model.Dimension, initialize=0, mutable=True)
        model.linear_correction_b = Param(initialize=0, mutable=True)
        model.cpwl_neg_k = Param(model.Dimension, mutable=True)
        model.cpwl_neg_b = Param(mutable=True)
        model.quadr_neg_k = Param(model.Dimension, mutable=True)
        model.quadr_neg_b = Param(mutable=True)
        def f_con(m):
            return model.f == model.cpwl_pos - model.cpwl_neg + model.quadr_pos - \
                   model.quadr_neg + model.linear_correction
        model.f_con=Constraint(rule=f_con)
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
        def quadr_pos_calc(m):
            return m.quadr_pos == sum([self.pos_quadratic_A[i,j]*input[self.input_variables[i]]*\
                                       input[self.input_variables[j]] for j in m.Dimension for i in m.Dimension])+\
                           sum([self.pos_quadratic_b[i]* input[self.input_variables[i]] for i in m.Dimension])+\
                  self.pos_quadratic_c
        model.quadr_pos_calc = Constraint(rule=quadr_pos_calc)
        def quadr_neg_calc(m):
            return m.quadr_neg == sum([m.quadr_neg_k[i] * \
                                      input[self.input_variables[i]] for i in m.Dimension]) + \
                   m.quadr_neg_b
        model.quadr_neg_calc = Constraint(rule=quadr_neg_calc)
        def linear_correction_calc(m):
            return m.linear_correction == sum([m.linear_correction_k[i] * \
                                       input[self.input_variables[i]] for i in m.Dimension]) + \
                   m.linear_correction_b
        model.linear_correction_calc = Constraint(rule=linear_correction_calc)

    def get_linear_correction(self, model):
        linear_correction_k = np.zeros((len(model.Dimension),))
        linear_correction_b = 0
        for i in model.Dimension:
            linear_correction_k[i] = value(model.linear_correction_k[i])
        linear_correction_b = value(model.linear_correction_b)
        return linear_correction_k, linear_correction_b

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
        model.linear_correction_b = modifier_epsilon - np.dot(modifier_lambda_vector,base_point_vector)
        
    def build_parameter_estimator(self, model):
        pass
    


class QuadraticBoostedDCCPWL_RTOObject(object):
    def __init__(self, dc_cpwl_functions, mvs, cvs):
        '''

        :param dc_cpwl_functions:
        :param mvs:  subset of all inputs
        :param cvs:
        '''
        self.dc_cpwl_functions = dc_cpwl_functions
        self.mvs = mvs
        self.cvs = cvs
        self.tol = 1e-5
        self.subproblem_max_iter = 100

    def build(self, problem_description):
        optimizer_model = ConcreteModel()
        optimizer_model.InputIndex = Set(initialize=self.mvs)
        optimizer_model.OutputIndex = Set(initialize=self.cvs)
        optimizer_model.input = Var(optimizer_model.InputIndex)
        optimizer_model.output = Block(optimizer_model.OutputIndex)
        for name, dc_cpwl in self.dc_cpwl_functions.items():
            dc_cpwl.build_optimizer(optimizer_model.output[name],optimizer_model.input)
        #TODO: direction of inequality_constraint, and scaling
        def inequality_constraint(m):
            for con_name in self.cvs[1:]:
                yield m.output[con_name].f <= 0
        optimizer_model.inequality_constraint=ConstraintList(rule=inequality_constraint)

        # TODO: direction of obj, and scaling
        def obj(m):
            return m.output[self.cvs[0]].f
        optimizer_model.obj = Objective(rule=obj, sense=minimize)
        self.optimizer_model = optimizer_model

        feasibility_problem_model = ConcreteModel()
        feasibility_problem_model.InputIndex = Set(initialize=self.mvs)
        feasibility_problem_model.OutputIndex = Set(initialize=self.cvs)
        feasibility_problem_model.input = Var(feasibility_problem_model.InputIndex)
        feasibility_problem_model.output = Block(feasibility_problem_model.OutputIndex)
        feasibility_problem_model.slack = Var(feasibility_problem_model.OutputIndex, within=NonNegativeReals)
        for name, dc_cpwl in self.dc_cpwl_functions.items():
            dc_cpwl.build_optimizer(feasibility_problem_model.output[name], feasibility_problem_model.input)

        # TODO: direction of inequality_constraint, and scaling
        def inequality_constraint(m):
            for con_name in self.cvs[1:]:
                yield m.output[con_name].f-m.slack[con_name] <= 0
        feasibility_problem_model.inequality_constraint = ConstraintList(rule=inequality_constraint)

        # TODO: direction of obj, and scaling
        def obj(m):
            return sum([m.slack[con_name] for con_name in self.cvs[1:]])
        feasibility_problem_model.obj = Objective(rule=obj, sense=minimize)
        self.feasibility_problem_model = feasibility_problem_model

    def update_modifiers(self, modifiers, base_point):
        print(modifiers)
        for cv in self.cvs:
            output_model = self.optimizer_model.output[cv]
            modifier_epsilon = modifiers[(cv, None)]
            modifier_lambda = {}
            for mv in self.mvs:
                modifier_lambda[mv] = modifiers[(cv, mv)]
            self.dc_cpwl_functions[cv].set_linear_correction(output_model, \
                                modifier_epsilon, modifier_lambda, base_point)
        for cv in self.cvs:
            output_model = self.feasibility_problem_model.output[cv]
            modifier_epsilon = modifiers[(cv, None)]
            modifier_lambda = {}
            for mv in self.mvs:
                modifier_lambda[mv] = modifiers[(cv, mv)]
            self.dc_cpwl_functions[cv].set_linear_correction(output_model, \
                                                             modifier_epsilon, modifier_lambda, base_point)

    def get_input_vector(self):
        input_vector = np.zeros((len(self.mvs),))
        for i,mv in enumerate(self.mvs):
            input_vector[i] = value(self.optimizer_model.input[mv])
        return input_vector

    def get_input_vector_from_feas_problem(self):
        input_vector = np.zeros((len(self.mvs),))
        for i,mv in enumerate(self.mvs):
            input_vector[i] = value(self.feasibility_problem_model.input[mv])
        return input_vector

    def optimize(self, spec_values, starting_point):
        '''

        :param input_values: fixed input
        :return:
        '''
        whole_input = {}
        for k,v in spec_values.items():
            whole_input[k] = v
            self.optimizer_model.input[k].fix(v)
            self.feasibility_problem_model.input[k].fix(v)
        for k,v in starting_point.items():
            whole_input[k] = v
            self.feasibility_problem_model.input[k] = v

        #-------------- solve feasibility problem ----------------
        if len(self.cvs) > 1:
            prev_u = self.get_input_vector_from_feas_problem()
            solve_status = False
            for i in range(self.subproblem_max_iter):
                for name, dc_cpwl in self.dc_cpwl_functions.items():
                    dc_cpwl.set_optimizer_concave_majorization(self.feasibility_problem_model.output[name], whole_input)
                results = self.solver.solve(self.feasibility_problem_model, tee=self.tee)
                if not ((results.solver.status == SolverStatus.ok) and (
                        results.solver.termination_condition == TerminationCondition.optimal)):
                    print(results.solver.termination_condition)
                    raise Exception("QCQP solving fails.")
                for mv in self.mvs:
                    whole_input[mv] = value(self.feasibility_problem_model.input[mv])
                this_u = self.get_input_vector_from_feas_problem()
                if np.linalg.norm(this_u - prev_u) < self.tol:
                    solve_status = True
                    break
                prev_u = this_u
        print(value(self.feasibility_problem_model.obj))
        # -------------- solve optimization problem ----------------
        for k,v in whole_input.items():
             self.optimizer_model.input[k] = v
        prev_u = self.get_input_vector()
        solve_status = False
        for i in range(self.subproblem_max_iter):
            for name, dc_cpwl in self.dc_cpwl_functions.items():
                dc_cpwl.set_optimizer_concave_majorization(self.optimizer_model.output[name], whole_input)
            results = self.solver.solve(self.optimizer_model, tee=self.tee)
            if not ((results.solver.status == SolverStatus.ok) and (
                    results.solver.termination_condition == TerminationCondition.optimal)):
                print(results.solver.termination_condition)
                raise Exception("QCQP solving fails.")
            for mv in self.mvs:
                whole_input[mv] = value(self.optimizer_model.input[mv])
            this_u = self.get_input_vector()
            if np.linalg.norm(this_u - prev_u) < self.tol:
                solve_status = True
                break
            prev_u = this_u
        return whole_input, solve_status

    def simulate(self, point, with_modifier=False):
        ret = {}
        for cv in self.cvs:
            ret[cv]=self.dc_cpwl_functions[cv].simulate(point, with_modifier, self.optimizer_model.output[cv])
        solve_status = True
        return ret, solve_status

    def set_model_parameters(self, parameter_value):
        pass

    def set_solver(self, solver, tee, default_options):
        self.solver = solver
        self.tee = tee
        for k, v in default_options.items():
            self.solver.options[k] = v