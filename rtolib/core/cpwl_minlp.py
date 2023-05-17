from rtolib.core.dc_cpwl_model import QuadraticBoostedDCCPWLFunctionSubgrad,\
    QuadraticBoostedDCCPWL_RTOObjectSubgrad
from pyomo.environ import *
from rtolib.core import ProblemDescription

class QCPWLFunction_MINLP(QuadraticBoostedDCCPWLFunctionSubgrad):
    def __init__(self):
        super().__init__()
        self.bigM = 10000

    def set_optimizer_concave_majorization(self, model, input):
        return NotImplementedError()

    def build_optimizer_mainbody(self, model, input):
        '''

        :param input: the input variable of the parent model
        :param model: assumes to be a block
        :param output_name:
        :return:
        '''
        model.PosSeg = Set(initialize=range(self.no_pos_segment))
        model.NegSeg = Set(initialize=range(self.no_neg_segment))
        model.Dimension = Set(initialize=range(self.dimension))

        model.cpwl_pos = Var()
        model.cpwl_neg = Var()
        model.linear_correction = Var()
        model.linear_correction_b = Param(initialize=0, mutable=True)
        model.base_point = Param(model.Dimension, mutable=True, initialize=0)
        model.cpwl_neg_k = Param(model.Dimension, mutable=True)
        model.cpwl_neg_b = Param(mutable=True)
        model.cpwl_pos_seg_active = Var(model.PosSeg, initialize=0, within=Boolean)
        model.cpwl_neg_seg_active = Var(model.NegSeg, initialize=0, within=Boolean)
        def quadr_pos_calc(m):
            return sum([self.pos_quadratic_A[i,j]*input[self.input_variables[i]]*\
                                       input[self.input_variables[j]] for j in m.Dimension for i in m.Dimension])+\
                           sum([self.pos_quadratic_b[i]* input[self.input_variables[i]] for i in m.Dimension])+\
                  self.pos_quadratic_c
        model.quadr_pos = Expression(rule=quadr_pos_calc)
        def quadr_neg_calc(m):
            return sum([self.neg_quadratic_A[i,j]*input[self.input_variables[i]]*\
                                       input[self.input_variables[j]] for j in m.Dimension for i in m.Dimension])+\
                           sum([self.neg_quadratic_b[i]* input[self.input_variables[i]] for i in m.Dimension])+\
                  self.neg_quadratic_c
        model.quadr_neg = Expression(rule=quadr_neg_calc)
        def cpwl_pos_calc1(m, s):
            return m.cpwl_pos >= sum([self.pos_cpwl_coeff[s,i]* \
                                     input[self.input_variables[i]] for i in m.Dimension])+ \
                   self.pos_cpwl_coeff[s,self.dimension]
        model.cpwl_pos_calc1=Constraint(model.PosSeg, rule=cpwl_pos_calc1)
        def cpwl_pos_calc2(m, s):
            return m.cpwl_pos <= sum([self.pos_cpwl_coeff[s,i]* \
                                     input[self.input_variables[i]] for i in m.Dimension])+ \
                   self.pos_cpwl_coeff[s,self.dimension] +\
                      self.bigM*(1-m.cpwl_pos_seg_active[s])
        model.cpwl_pos_calc2=Constraint(model.PosSeg, rule=cpwl_pos_calc2)
        def cpwl_neg_calc1(m, s):
            return m.cpwl_neg >= sum([self.neg_cpwl_coeff[s,i]* \
                                     input[self.input_variables[i]] for i in m.Dimension])+ \
                   self.neg_cpwl_coeff[s,self.dimension]
        model.cpwl_neg_calc1 = Constraint(model.NegSeg, rule=cpwl_neg_calc1)
        def cpwl_neg_calc2(m, s):
            return m.cpwl_neg <= sum([self.neg_cpwl_coeff[s,i]* \
                                     input[self.input_variables[i]] for i in m.Dimension])+ \
                   self.neg_cpwl_coeff[s,self.dimension] +\
                      self.bigM*(1-m.cpwl_neg_seg_active[s])
        model.cpwl_neg_calc2 = Constraint(model.NegSeg, rule=cpwl_neg_calc2)
        def pos_active_seg_con(m):
            return sum([m.cpwl_pos_seg_active[s] for s in m.PosSeg]) == 1
        model.pos_active_seg_con = Constraint(rule=pos_active_seg_con)
        def neg_active_seg_con(m):
            return sum([m.cpwl_neg_seg_active[s] for s in m.NegSeg]) == 1
        model.neg_active_seg_con = Constraint(rule=neg_active_seg_con)
        def f_con(m):
            return model.cpwl_pos - model.cpwl_neg + model.quadr_pos - \
                   model.quadr_neg + model.linear_correction
        model.f=Expression(rule=f_con)

class QCPWL_RTOObjectSubgradMINLP(QuadraticBoostedDCCPWL_RTOObjectSubgrad):
    def build(self, problem_description, constraint_scaling=1):
        assert isinstance(problem_description, ProblemDescription)
        self.problem_description = problem_description
        self.problem_description.scaling_factors['validity_con'] = 1

        subproblem_model = ConcreteModel()
        subproblem_model.InputIndex = Set(initialize=self.inputs)
        subproblem_model.OutputIndex = Set(initialize=self.cvs)
        def bounds_func(m, i):
            return self.problem_description.bounds[i]
        subproblem_model.input = Var(subproblem_model.InputIndex, bounds=bounds_func)
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
                if con_var_name == "validity_con":
                    continue
                yield (m.output[con_var_name].f /
                       self.problem_description.scaling_factors[con_var_name] \
                       - getattr(subproblem_model, "con_vio_" + con_var_name)) * constraint_scaling <= 0
        subproblem_model.violation_cons = ConstraintList(rule=_violation_cons)

        # merit function
        def _tr_infeasibility_sum(m):
            ret = 0
            for con_var_name in self.cvs[1:]:
                if con_var_name == "validity_con":
                    continue
                ret += getattr(m, "con_vio_" + con_var_name) / \
                       self.problem_description.scaling_factors[con_var_name]
            return ret - m.infeasibility_sum == 0
        subproblem_model.tr_infeasibility_sum = Constraint(rule=_tr_infeasibility_sum)

        def _tr_merit_function(m):
            ret = m.output[self.cvs[0]].f
            ret += m.infeasibility_sum * m.tr_penalty_coeff
            return ret
        subproblem_model.tr_merit_function = Expression(rule=_tr_merit_function)

        def _validity_con(m):
            return m.output["validity_con"].f <= 0
        subproblem_model.validity_con = Constraint(rule=_validity_con)

        def optimization_obj(m):
            return m.tr_merit_function
        subproblem_model.optimization_obj = Objective(rule=optimization_obj, sense=minimize)

        self.subproblem_model = subproblem_model

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

        results = self.solver.solve(self.subproblem_model, tee=self.tee) #self.tee
        if not ((results.solver.status == SolverStatus.ok) and (
                results.solver.termination_condition == TerminationCondition.optimal)):
            print(results.solver.termination_condition)
            raise Exception("QCQP solving fails.")
        for mv in self.inputs:
            whole_input[mv] = value(self.subproblem_model.input[mv])
        print("after optimization:")
        print(whole_input)
        solve_status = True

        output = {}
        for cv in self.cvs:
            output[cv] = value(self.subproblem_model.output[cv].f)
        print("minlp output:")
        print(output)
        output = self.simulate(whole_input, with_modifier=True)
        print("simulation output:")
        print(output)

        return whole_input, solve_status