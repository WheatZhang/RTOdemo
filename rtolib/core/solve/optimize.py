from pyomo.environ import ConcreteModel, Var, Objective, Expression,\
    TerminationCondition, value, SolverStatus, minimize, ConstraintList, Constraint,\
    NonNegativeReals, sqrt
from rtolib.core.pyomo_model import PyomoModel, PyomoModelWithModifiers,\
    PyomoModelSolvingStatus
from rtolib.core.basic import ProblemDescription
import rtolib.util.init_value as init_value
import numpy


class Optimizer():
    def optimize(self, param_values=None):
        return NotImplementedError()


class PyomoOptimizer(Optimizer):
    def __init__(self, pyomo_model):
        assert isinstance(pyomo_model, PyomoModel)
        self.pyomo_model = pyomo_model
        self.model = ConcreteModel()
        self.homotopy_var = []
        for var_name in self.pyomo_model.input_variables.keys():
            self.homotopy_var.append(var_name)

    def set_base_point(self, base_point):
        if not isinstance(self.pyomo_model, PyomoModelWithModifiers):
            raise AttributeError('this model does not have base point')
        self.pyomo_model.set_base_point(self.model, base_point)

    def set_solver(self, solver, tee, default_options):
        self.solver = solver
        self.tee = tee
        for k, v in default_options.items():
            self.solver.options[k] = v

    def build(self, problem_description, homotopy_factor=1000):
        assert isinstance(problem_description, ProblemDescription)
        self.problem_description = problem_description
        self.pyomo_model.build(self.model)
        # print(self.pyomo_model.parameters['Fb_profit_lam'].__call__(self.model))
        # set constraints
        for con_name in self.problem_description.symbol_list['CON']:
            con_var = self.pyomo_model.output_variables[con_name].__call__(self.model)
            setattr(self.model, "con_"+con_name, Constraint(expr=(con_var<=0)))
        # set obj
        def _objective(m):
            return self.pyomo_model.output_variables[self.problem_description.symbol_list['OBJ']].__call__(self.model)
        self.model.objective_function = Objective(rule=_objective, sense=minimize)

        # create variables for homotopy
        for ip in self.pyomo_model.input_variables.keys():
            setattr(self.model, ip + "_homo", Var(initialize=0))
            getattr(self.model, ip + "_homo").fixed = True
        for theta in self.pyomo_model.parameters.keys():
            setattr(self.model, theta + "_homo", Var(initialize=0))
            getattr(self.model, theta + "_homo").fixed = True
        def _homotopy_expr(m):
            ip_scaling_coeffs = self.problem_description.scaling_factors
            theta_scaling_coeffs = self.pyomo_model.parameter_scaling_factors
            return sum([(getattr(m, k + "_homo") - \
                         self.pyomo_model.input_variables[k].__call__(self.model)) ** 2 \
                        / ip_scaling_coeffs[k] / ip_scaling_coeffs[k] \
                        for k in self.pyomo_model.input_variables.keys()])+ \
                   sum([(getattr(m, k + "_homo") - \
                         self.pyomo_model.parameters[k].__call__(self.model)) ** 2 \
                        / theta_scaling_coeffs[k] / theta_scaling_coeffs[k] \
                        for k in self.pyomo_model.parameters.keys()])
        self.model.homotopy_expr = Expression(rule=_homotopy_expr)

        def _homotopy_optimization_obj(m):
            return self.model.homotopy_expr*homotopy_factor+ m.objective_function
        self.model.homotopy_optimization_obj = Objective(rule=_homotopy_optimization_obj, sense=minimize)
        # load initial value file
        init_value.load_init_from_template(self.model, self.pyomo_model.initial_value_file, \
                                                       ignore_init_mismatch=True)

    def optimize(self, input_values, param_values=None, use_homo=True):
        # autocomplete
        if param_values is None:
            param_values = {}
            for para_name in self.pyomo_model.parameters.keys():
                param_values[para_name] = self.pyomo_model.default_value[para_name]
        else:
            for para_name in self.pyomo_model.parameters.keys():
                if para_name not in param_values.keys():
                    param_values[para_name] = self.pyomo_model.default_value[para_name]

        model = self.pyomo_model.parameters
        # Fix inputs and parameters
        if use_homo:
            for var_name in self.pyomo_model.input_variables.keys():
                if var_name in input_values.keys():
                    if var_name in self.homotopy_var:
                        getattr(self.model, var_name + "_homo").fix(input_values[var_name])
                    else:
                        var = self.pyomo_model.input_variables[var_name].__call__(self.model)
                        var.fix(input_values[var_name])
                        getattr(self.model, var_name + "_homo").fix(input_values[var_name])
                else:
                    getattr(self.model, var_name + "_homo").fixed=False
            for var_name in self.pyomo_model.parameters.keys():
                var = self.pyomo_model.parameters[var_name].__call__(self.model)
                if var_name in self.homotopy_var:
                    var.fixed = False
                    getattr(self.model, var_name + "_homo").fix(param_values[var_name])
                else:
                    var.fix(param_values[var_name])
                    getattr(self.model, var_name + "_homo").fix(param_values[var_name])

        # TODO: backup initial value
        # homotopy optimization
            self.model.homotopy_optimization_obj.activate()
            self.model.objective_function.deactivate()

            results = self.solver.solve(self.model, tee=self.tee)
            if not ((results.solver.status == SolverStatus.ok) and (
                    results.solver.termination_condition == TerminationCondition.optimal)):
                print(results.solver.termination_condition)
                # TODO: restore init value

        # usual optimization
        for var_name in self.pyomo_model.input_variables.keys():
            var = self.pyomo_model.input_variables[var_name].__call__(self.model)
            if var_name in input_values.keys():
                var.fix(input_values[var_name])
            else:
                var.fixed=False
        for var_name in self.pyomo_model.parameters.keys():
            var = self.pyomo_model.parameters[var_name].__call__(self.model)
            var.fix(param_values[var_name])
        self.model.homotopy_optimization_obj.deactivate()
        self.model.objective_function.activate()

        results = self.solver.solve(self.model, tee=self.tee)

        if not ((results.solver.status == SolverStatus.ok) and (
                results.solver.termination_condition == TerminationCondition.optimal)):
            print(results.solver.termination_condition)
            solve_status = PyomoModelSolvingStatus.OPTIMIZATION_FAILED
        else:
            solve_status = PyomoModelSolvingStatus.OK
        inputs = {}
        for ip in self.pyomo_model.input_variables.keys():
            var = self.pyomo_model.input_variables[ip].__call__(self.model)
            inputs[ip] = value(var)
        return inputs, solve_status

class TrustRegionOptimizer(PyomoOptimizer):
    def build(self, problem_description, homotopy_factor=1000):
        assert isinstance(problem_description, ProblemDescription)
        self.problem_description = problem_description
        self.pyomo_model.build(self.model)
        # print(self.pyomo_model.parameters['Fb_profit_lam'].__call__(self.model))

        # set tr con
        for ip in self.pyomo_model.input_variables.keys():
            setattr(self.model, ip + "_tr_base", Var(initialize=0))
            getattr(self.model, ip + "_tr_base").fixed = True
        setattr(self.model, "tr_radius", Var(initialize=0))
        getattr(self.model, "tr_radius").fixed = True
        def _trust_region_con(m):
            r=0
            for ip in self.pyomo_model.input_variables.keys():
                r += (self.pyomo_model.input_variables[ip].__call__(m)-getattr(m, ip + "_tr_base"))**2/\
                     (self.problem_description.scaling_factors[ip]**2)
            return r <= m.tr_radius**2
        self.model.trust_region_cons = Constraint(rule=_trust_region_con)
        # set constraints
        for con_name in self.problem_description.symbol_list['CON']:
            con_var = self.pyomo_model.output_variables[con_name].__call__(self.model)
            setattr(self.model, "con_" + con_name, Constraint(expr=(con_var <= 0)))
        # set obj
        def _objective(m):
            return self.pyomo_model.output_variables[self.problem_description.symbol_list['OBJ']].__call__(self.model)
        self.model.objective_function = Objective(rule=_objective, sense=minimize)

        # load initial value file
        init_value.load_init_from_template(self.model, self.pyomo_model.initial_value_file, \
                                                       ignore_init_mismatch=True)

    def optimize(self, input_values, tr_radius, tr_base, param_values=None, use_homo=True):
        # autocomplete
        if param_values is None:
            param_values = {}
            for para_name in self.pyomo_model.parameters.keys():
                param_values[para_name] = self.pyomo_model.default_value[para_name]
        else:
            for para_name in self.pyomo_model.parameters.keys():
                if para_name not in param_values.keys():
                    param_values[para_name] = self.pyomo_model.default_value[para_name]

        model = self.pyomo_model.parameters
        # Fix inputs and parameters
        if use_homo:
            raise NotImplementedError("homotopy optimization not supported")

        # usual optimization
        for var_name in self.pyomo_model.input_variables.keys():
            var = self.pyomo_model.input_variables[var_name].__call__(self.model)
            if var_name in input_values.keys():
                var.fix(input_values[var_name])
            else:
                var.fixed=False
        for var_name in self.pyomo_model.parameters.keys():
            var = self.pyomo_model.parameters[var_name].__call__(self.model)
            var.fix(param_values[var_name])

        # trust-region optimization
        for var_name in self.pyomo_model.input_variables.keys():
            if var_name in tr_base.keys():
                getattr(self.model, var_name + "_tr_base").fix(tr_base[var_name])
            else:
                getattr(self.model, var_name + "_tr_base").fixed=False
        getattr(self.model, "tr_radius").fix(tr_radius)

        # self.model.homotopy_optimization_obj.deactivate()
        self.model.objective_function.activate()

        results = self.solver.solve(self.model, tee=self.tee)

        if not ((results.solver.status == SolverStatus.ok) and (
                results.solver.termination_condition == TerminationCondition.optimal)):
            print(results.solver.termination_condition)
            solve_status = PyomoModelSolvingStatus.OPTIMIZATION_FAILED
        else:
            solve_status = PyomoModelSolvingStatus.OK
        inputs = {}
        for ip in self.pyomo_model.input_variables.keys():
            var = self.pyomo_model.input_variables[ip].__call__(self.model)
            inputs[ip] = value(var)
        return inputs, solve_status


class PenaltyTrustRegionOptimizer(PyomoOptimizer):
    def build(self, problem_description, homotopy_factor=1000):
        assert isinstance(problem_description, ProblemDescription)
        self.problem_description = problem_description
        self.pyomo_model.build(self.model)
        # print(self.pyomo_model.parameters['Fb_profit_lam'].__call__(self.model))

        # set penalty coefficient
        setattr(self.model, "dummy_var_for_sqrt_c", Var(initialize=0, within=NonNegativeReals))
        setattr(self.model, "tr_penalty_coeff", Var(initialize=0))
        getattr(self.model, "tr_penalty_coeff").fixed = True

        # inequality violation measure
        for con_var_name in self.problem_description.symbol_list['CON']:
            setattr(self.model, "con_vio_"+con_var_name, Var(initialize=0, within=NonNegativeReals))
        def _violation_cons(m):
            for con_var_name in self.problem_description.symbol_list['CON']:
                yield self.pyomo_model.output_variables[con_var_name].__call__(m)-getattr(self.model, "con_vio_"+con_var_name) <=0
        self.model.violation_cons = ConstraintList(rule=_violation_cons)

        # merit function
        def _tr_infeasibility_sq(m):
            ret = 0
            for con_var_name in self.problem_description.symbol_list['CON']:
                ret += getattr(self.model, "con_vio_"+con_var_name)**2 / \
                      self.problem_description.scaling_factors[con_var_name]**2
            return ret
        self.model.tr_infeasibility_sq = Expression(rule=_tr_infeasibility_sq)

        def c_sqrt_cons(m):
            return m.dummy_var_for_sqrt_c**2 == m.tr_infeasibility_sq
        self.model.c_sqrt_cons = Constraint(rule=c_sqrt_cons)

        def _tr_merit_function(m):
            obj_var_name = self.problem_description.symbol_list['OBJ']
            ret=self.pyomo_model.output_variables[obj_var_name].__call__(m)/ \
                self.problem_description.scaling_factors[obj_var_name]
            ret += m.dummy_var_for_sqrt_c * m.tr_penalty_coeff
            # ret += sqrt(m.tr_infeasibility_sq) * m.tr_penalty_coeff
            return ret
        self.model.tr_merit_function = Expression(rule=_tr_merit_function)

        # set tr con
        for ip in self.pyomo_model.input_variables.keys():
            setattr(self.model, ip + "_tr_base", Var(initialize=0))
            getattr(self.model, ip + "_tr_base").fixed = True
        setattr(self.model, "tr_radius", Var(initialize=0))
        getattr(self.model, "tr_radius").fixed = True
        def _trust_region_con(m):
            r = 0
            for ip in self.pyomo_model.input_variables.keys():
                r += (self.pyomo_model.input_variables[ip].__call__(m) - getattr(m, ip + "_tr_base")) ** 2 / \
                     (self.problem_description.scaling_factors[ip] ** 2)
            return r <= m.tr_radius ** 2
        self.model.trust_region_cons = Constraint(rule=_trust_region_con)

        # set obj
        def _objective(m):
            return m.tr_merit_function
        self.model.objective_function = Objective(rule=_objective, sense=minimize)

        # load initial value file
        init_value.load_init_from_template(self.model, self.pyomo_model.initial_value_file, \
                                                       ignore_init_mismatch=True)

    def set_penalty_coeff(self, sigma):
        self.model.tr_penalty_coeff.fix(sigma)

    def optimize(self, input_values, tr_radius, tr_base, param_values=None, use_homo=True):
        # autocomplete
        if param_values is None:
            param_values = {}
            for para_name in self.pyomo_model.parameters.keys():
                param_values[para_name] = self.pyomo_model.default_value[para_name]
        else:
            for para_name in self.pyomo_model.parameters.keys():
                if para_name not in param_values.keys():
                    param_values[para_name] = self.pyomo_model.default_value[para_name]

        model = self.pyomo_model.parameters
        # Fix inputs and parameters
        if use_homo:
            raise NotImplementedError("homotopy optimization not supported")

        # usual optimization
        for var_name in self.pyomo_model.input_variables.keys():
            var = self.pyomo_model.input_variables[var_name].__call__(self.model)
            if var_name in input_values.keys():
                var.fix(input_values[var_name])
            else:
                var.fixed=False
        for var_name in self.pyomo_model.parameters.keys():
            var = self.pyomo_model.parameters[var_name].__call__(self.model)
            var.fix(param_values[var_name])

        # trust-region optimization
        for var_name in self.pyomo_model.input_variables.keys():
            if var_name in tr_base.keys():
                getattr(self.model, var_name + "_tr_base").fix(tr_base[var_name])
            else:
                getattr(self.model, var_name + "_tr_base").fixed = False
        getattr(self.model, "tr_radius").fix(tr_radius)
        # self.model.dummy_var_for_sqrt_c = numpy.sqrt(value(self.model.tr_infeasibility_sq))

        # self.model.homotopy_optimization_obj.deactivate()
        self.model.objective_function.activate()

        results = self.solver.solve(self.model, tee=self.tee)
        # print(value(self.model.con_vio_con))
        # print(value(self.model.tr_infeasibility_sq))
        # print(value(self.model.dummy_var_for_sqrt_c))

        if not ((results.solver.status == SolverStatus.ok) and (
                results.solver.termination_condition == TerminationCondition.optimal)):
            print(results.solver.termination_condition)
            solve_status = PyomoModelSolvingStatus.OPTIMIZATION_FAILED
        else:
            solve_status = PyomoModelSolvingStatus.OK
        inputs = {}
        for ip in self.pyomo_model.input_variables.keys():
            var = self.pyomo_model.input_variables[ip].__call__(self.model)
            inputs[ip] = value(var)
        return inputs, solve_status


class CompoStepTrustRegionOptimizer(PyomoOptimizer):
    def build(self, problem_description, homotopy_factor=1000):
        assert isinstance(problem_description, ProblemDescription)
        self.problem_description = problem_description
        self.pyomo_model.build(self.model)
        # print(self.pyomo_model.parameters['Fb_profit_lam'].__call__(self.model))

        # set base infeasibility
        setattr(self.model, "base_infeasibility", Var(initialize=0))
        getattr(self.model, "base_infeasibility").fixed = True

        # inequality violation measure
        for con_var_name in self.problem_description.symbol_list['CON']:
            setattr(self.model, "con_vio_"+con_var_name, Var(initialize=0, within=NonNegativeReals))
        def _violation_cons(m):
            for con_var_name in self.problem_description.symbol_list['CON']:
                yield self.pyomo_model.output_variables[con_var_name].__call__(m)-getattr(self.model, "con_vio_"+con_var_name) <=0
        self.model.violation_cons = ConstraintList(rule=_violation_cons)

        # c^2
        def _tr_infeasibility_sq(m):
            ret = 0
            for con_var_name in self.problem_description.symbol_list['CON']:
                ret += getattr(self.model, "con_vio_"+con_var_name)**2 / \
                      self.problem_description.scaling_factors[con_var_name]**2
            return ret
        self.model.tr_infeasibility_sq = Expression(rule=_tr_infeasibility_sq)

        # set tr con
        for ip in self.pyomo_model.input_variables.keys():
            setattr(self.model, ip + "_tr_base", Var(initialize=0))
            getattr(self.model, ip + "_tr_base").fixed = True
        setattr(self.model, "tr_radius", Var(initialize=0))
        getattr(self.model, "tr_radius").fixed = True
        def _trust_region_con(m):
            r=0
            for ip in self.pyomo_model.input_variables.keys():
                r += (self.pyomo_model.input_variables[ip].__call__(m)-getattr(m, ip + "_tr_base"))**2/\
                     (self.problem_description.scaling_factors[ip]**2)
            return r <= m.tr_radius**2
        self.model.trust_region_cons = Constraint(rule=_trust_region_con)

        # set tangential step constraints
        def _tangential_step_feasibility(m):
            return m.tr_infeasibility_sq <= m.base_infeasibility
        self.model.tangential_step_feasibility = Constraint(rule=_tangential_step_feasibility)

        # set obj
        def _normal_step_objective(m):
            ret = m.tr_infeasibility_sq
            return ret
        self.model.normal_step_objective = Objective(rule=_normal_step_objective, sense=minimize)

        def _tangential_step_objective(m):
            obj_var_name = self.problem_description.symbol_list['OBJ']
            ret = self.pyomo_model.output_variables[obj_var_name].__call__(m)/ \
                           self.problem_description.scaling_factors[obj_var_name]
            return ret
        self.model.tangential_step_objective = Objective(rule=_tangential_step_objective, sense=minimize)

        # load initial value file
        init_value.load_init_from_template(self.model, self.pyomo_model.initial_value_file, \
                                                       ignore_init_mismatch=True)


    def optimize(self, input_values, tr_radius, tr_base, xi_N, param_values=None, use_homo=True):
        # autocomplete
        if param_values is None:
            param_values = {}
            for para_name in self.pyomo_model.parameters.keys():
                param_values[para_name] = self.pyomo_model.default_value[para_name]
        else:
            for para_name in self.pyomo_model.parameters.keys():
                if para_name not in param_values.keys():
                    param_values[para_name] = self.pyomo_model.default_value[para_name]

        model = self.pyomo_model.parameters
        # Fix inputs and parameters
        if use_homo:
            raise NotImplementedError("homotopy optimization not supported")

        # usual optimization
        for var_name in self.pyomo_model.input_variables.keys():
            var = self.pyomo_model.input_variables[var_name].__call__(self.model)
            if var_name in input_values.keys():
                var.fix(input_values[var_name])
            else:
                var.fixed=False
        for var_name in self.pyomo_model.parameters.keys():
            var = self.pyomo_model.parameters[var_name].__call__(self.model)
            var.fix(param_values[var_name])

        # trust-region optimization - normal step
        for var_name in self.pyomo_model.input_variables.keys():
            if var_name in tr_base.keys():
                getattr(self.model, var_name + "_tr_base").fix(tr_base[var_name])
            else:
                getattr(self.model, var_name + "_tr_base").fixed = False
        getattr(self.model, "tr_radius").fix(tr_radius*xi_N)

        # self.model.homotopy_optimization_obj.deactivate()
        self.model.normal_step_objective.activate()
        self.model.tangential_step_objective.deactivate()
        self.model.tangential_step_feasibility.deactivate()

        results = self.solver.solve(self.model, tee=self.tee)
        # results = self.solver.solve(self.model, tee=True)

        if not ((results.solver.status == SolverStatus.ok) and (
                results.solver.termination_condition == TerminationCondition.optimal)):
            print(results.solver.termination_condition)
            solve_status = PyomoModelSolvingStatus.OPTIMIZATION_FAILED
        else:
            solve_status = PyomoModelSolvingStatus.OK

        # trust-region optimization - tangential step
        getattr(self.model, "tr_radius").fix(tr_radius)

        # self.model.homotopy_optimization_obj.deactivate()
        self.model.normal_step_objective.deactivate()
        self.model.tangential_step_objective.activate()
        self.model.tangential_step_feasibility.activate()
        self.model.base_infeasibility.fix(value(self.model.tr_infeasibility_sq))

        results = self.solver.solve(self.model, tee=self.tee)

        if not ((results.solver.status == SolverStatus.ok) and (
                results.solver.termination_condition == TerminationCondition.optimal)):
            print(results.solver.termination_condition)
            solve_status = PyomoModelSolvingStatus.OPTIMIZATION_FAILED
        else:
            solve_status = PyomoModelSolvingStatus.OK

        # return result
        inputs = {}
        for ip in self.pyomo_model.input_variables.keys():
            var = self.pyomo_model.input_variables[ip].__call__(self.model)
            inputs[ip] = value(var)
        return inputs, solve_status