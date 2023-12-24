from pyomo.environ import ConcreteModel, Var, Objective, Expression,\
    TerminationCondition, value, SolverStatus, minimize, ConstraintList, Constraint,\
    NonNegativeReals, sqrt
from rtolib.core.pyomo_model import PyomoModel, PyomoModelWithModifiers,\
    ModelSolvingStatus
from rtolib.core.black_box_model import BlackBoxModel, BlackBoxModelWithModifiers
from rtolib.core.basic import ProblemDescription
import rtolib.util.init_value as init_value
import numpy
from sko.PSO import PSO
from rtolib.util.misc import square_circle_mapping

class Optimizer():
    def optimize(self, param_values=None):
        return NotImplementedError()

    def adapt_to_bound_mv(self, optimized_input, mv_bounds, print_warning=False):
        bounded_input = {}
        # for k in self.problem_description.symbol_list['MV']:
        for k in optimized_input.keys():
            bounded_input[k] = optimized_input[k]
            if mv_bounds[k][0] is not None and optimized_input[k] <= mv_bounds[k][0]:
                if print_warning:
                    print("MV %s reaches its lower bound." % k)
                bounded_input[k] = mv_bounds[k][0]
            if mv_bounds[k][0] is not None and optimized_input[k] >= mv_bounds[k][1]:
                if print_warning:
                    print("MV %s reaches its upper bound." % k)
                bounded_input[k] = mv_bounds[k][1]
        return bounded_input



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

    def build(self, problem_description, homotopy_factor=100):
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
            init_value.load_init_from_template(self.model, self.pyomo_model.initial_value_file, \
                                               ignore_init_mismatch=True)
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
                print("in solving homotopy optimization problem")
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
            print("in solving original optimization problem")
            print(results.solver.termination_condition)
            solve_status = ModelSolvingStatus.OPTIMIZATION_FAILED
        else:
            solve_status = ModelSolvingStatus.OK
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
            solve_status = ModelSolvingStatus.OPTIMIZATION_FAILED
        else:
            solve_status = ModelSolvingStatus.OK
        inputs = {}
        for ip in self.pyomo_model.input_variables.keys():
            var = self.pyomo_model.input_variables[ip].__call__(self.model)
            inputs[ip] = value(var)
        return inputs, solve_status


class PenaltyTrustRegionOptimizer(PyomoOptimizer):
    def build(self, problem_description, constraint_scaling=1):
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
                yield (self.pyomo_model.output_variables[con_var_name].__call__(m)/self.problem_description.scaling_factors[con_var_name]\
                      -getattr(self.model, "con_vio_"+con_var_name))*constraint_scaling <=0
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
            return (m.dummy_var_for_sqrt_c**2 - m.tr_infeasibility_sq)*constraint_scaling==0
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
            solve_status = ModelSolvingStatus.OPTIMIZATION_FAILED
        else:
            solve_status = ModelSolvingStatus.OK
        inputs = {}
        for ip in self.pyomo_model.input_variables.keys():
            var = self.pyomo_model.input_variables[ip].__call__(self.model)
            inputs[ip] = value(var)
        return inputs, solve_status


class CompoStepTrustRegionOptimizer(PyomoOptimizer):
    def build(self, problem_description, constraint_scaling=1e3, objective_scaling=1e3):
        assert isinstance(problem_description, ProblemDescription)
        self.problem_description = problem_description
        self.pyomo_model.build(self.model)
        # print(self.pyomo_model.parameters['Fb_profit_lam'].__call__(self.model))

        # set base infeasibility
        setattr(self.model, "base_infeasibility_sq", Var(initialize=0))
        getattr(self.model, "base_infeasibility_sq").fixed = True

        # inequality violation measure
        for con_var_name in self.problem_description.symbol_list['CON']:
            setattr(self.model, "con_vio_"+con_var_name, Var(initialize=0, within=NonNegativeReals))
        def _violation_cons(m):
            for con_var_name in self.problem_description.symbol_list['CON']:
                yield (self.pyomo_model.output_variables[con_var_name].__call__(m)/self.problem_description.scaling_factors[con_var_name]\
                      -getattr(self.model, "con_vio_"+con_var_name))*constraint_scaling <=0
        self.model.violation_cons = ConstraintList(rule=_violation_cons)

        # c^2
        def _tr_infeasibility_sq(m):
            ret = 0
            for con_var_name in self.problem_description.symbol_list['CON']:
                ret += getattr(self.model, "con_vio_"+con_var_name)**2
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
            return (m.tr_infeasibility_sq - m.base_infeasibility_sq)*constraint_scaling**2<=0
        self.model.tangential_step_feasibility = Constraint(rule=_tangential_step_feasibility)

        # set obj
        def _normal_step_objective(m):
            ret = m.tr_infeasibility_sq*constraint_scaling
            return ret
        self.model.normal_step_objective = Objective(rule=_normal_step_objective, sense=minimize)

        def _tangential_step_objective(m):
            obj_var_name = self.problem_description.symbol_list['OBJ']
            ret = self.pyomo_model.output_variables[obj_var_name].__call__(m)/ \
                           self.problem_description.scaling_factors[obj_var_name]*objective_scaling
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
            solve_status = ModelSolvingStatus.OPTIMIZATION_FAILED
        else:
            solve_status = ModelSolvingStatus.OK

        input_after_normal_step = {}
        for ip in self.pyomo_model.input_variables.keys():
            var = self.pyomo_model.input_variables[ip].__call__(self.model)
            input_after_normal_step[ip] = value(var)

        # trust-region optimization - tangential step
        getattr(self.model, "tr_radius").fix(tr_radius)

        # self.model.homotopy_optimization_obj.deactivate()
        self.model.normal_step_objective.deactivate()
        self.model.tangential_step_objective.activate()
        self.model.tangential_step_feasibility.activate()
        # print(value(self.model.tr_infeasibility_sq))
        self.model.base_infeasibility_sq.fix(value(self.model.tr_infeasibility_sq))

        results = self.solver.solve(self.model, tee=self.tee)
        # print(value(self.model.tr_infeasibility_sq))

        if not ((results.solver.status == SolverStatus.ok) and (
                results.solver.termination_condition == TerminationCondition.optimal)):
            print(results.solver.termination_condition)
            solve_status = ModelSolvingStatus.OPTIMIZATION_FAILED
        else:
            solve_status = ModelSolvingStatus.OK

        # return result
        inputs = {}
        for ip in self.pyomo_model.input_variables.keys():
            var = self.pyomo_model.input_variables[ip].__call__(self.model)
            inputs[ip] = value(var)
        return inputs, input_after_normal_step, solve_status


class BlackBoxOptimizer(Optimizer):
    def __init__(self, black_box_model):
        assert isinstance(black_box_model, BlackBoxModel)
        self.black_box_model = black_box_model
    def build(self, problem_description):
        raise NotImplementedError()

    def optimize(self):
        raise NotImplementedError()

    def set_base_point(self, base_point):
        raise NotImplementedError()

    def set_solver_options(self, options):
        raise NotImplementedError()

class CompoStepTrustRegionBBMOptimizer(BlackBoxOptimizer):
    '''
    Black box optimizer
    '''
    def build(self, problem_description):
        assert isinstance(problem_description, ProblemDescription)
        self.problem_description = problem_description
        self.black_box_model.build()

    def set_base_point(self, base_point):
        if not isinstance(self.black_box_model, BlackBoxModelWithModifiers):
            raise AttributeError('this model does not have base point')
        self.black_box_model.set_base_point(base_point)

    def set_solver_options(self, solution_method_name, options={}):
        '''
        :param solution_method_name:
        :param options:

        :return:
        '''
        self.solution_method_name = solution_method_name
        if solution_method_name == "PSO":
            default_options = {}
            default_options["population"] = 20
            default_options["max_iter"] = 5
            default_options["inertia"] = 0.8
            default_options["c1"] = 0.5
            default_options["c2"] = 0.5
        for k,v in options.items():
            default_options[k] = v
        self.options = default_options

    def optimize_direct_solve_con_problem(self, tr_radius, tr_base, xi_N, modifiers_value, objective_scaling=1e3):
        # set modifiers
        assert(isinstance(self.black_box_model, BlackBoxModelWithModifiers))
        self.black_box_model.set_modifiers(modifiers_value)

        # get mv bounds
        lb = []
        ub = []
        for i, mv_name in enumerate(self.problem_description.symbol_list['MV']):
            lb.append(self.problem_description.bounds[mv_name][0])
            ub.append(self.problem_description.bounds[mv_name][1])

        def normal_step_problem_obj(x):
            input_dict = {}
            for i,mv_name in enumerate(self.problem_description.symbol_list['MV']):
                input_dict[mv_name] = x[i]
            output_dict = self.black_box_model.simulate(input_dict)
            infeasibility_sq = 0
            for con_var_name in self.problem_description.symbol_list['CON']:
                infeasibility_sq+=(output_dict[con_var_name]/max(0,self.problem_description.scaling_factors[con_var_name]))**2
            return infeasibility_sq

        def normal_step_problem_tr_con(x):
            radius_sq = 0
            for i, mv_name in enumerate(self.problem_description.symbol_list['MV']):
                radius_sq+=(x[i]-tr_base[mv_name])**2/(self.problem_description.scaling_factors[mv_name]**2)
            return radius_sq-(tr_radius*xi_N)**2

        # solve normal step optimization problem
        pso_normal_problem = PSO(func=normal_step_problem_obj, n_dim=len(lb), pop=self.options['population'], \
                                 max_iter=self.options["max_iter"], lb=lb, ub=ub, \
                                 w=self.options["inertia"], c1=self.options["c1"], \
                                 c2=self.options["c2"], constraint_ueq=(normal_step_problem_tr_con,))
        pso_normal_problem.run()
        x_after_normal_step = pso_normal_problem.gbest_x[0]
        tr_infeasibility_sq_after_normal_step = normal_step_problem_obj(x_after_normal_step)
        input_after_normal_step = {}
        for i, mv_name in enumerate(self.problem_description.symbol_list['MV']):
            input_after_normal_step[mv_name] = x_after_normal_step[i]
        print(tr_infeasibility_sq_after_normal_step)

        def tangential_step_problem_obj(x):
            input_dict = {}
            for i, mv_name in enumerate(self.problem_description.symbol_list['MV']):
                input_dict[mv_name] = x[i]
            output_dict = self.black_box_model.simulate(input_dict)
            return output_dict[self.problem_description.symbol_list['OBJ']]*objective_scaling

        def tangential_step_problem_c_con(x):
            input_dict = {}
            for i, mv_name in enumerate(self.problem_description.symbol_list['MV']):
                input_dict[mv_name] = x[i]
            output_dict = self.black_box_model.simulate(input_dict)
            infeasibility_sq = 0
            for con_var_name in self.problem_description.symbol_list['CON']:
                infeasibility_sq += (output_dict[con_var_name] / max(0,
                                 self.problem_description.scaling_factors[con_var_name]))**2
            return infeasibility_sq-tr_infeasibility_sq_after_normal_step

        def tangential_step_problem_tr_con(x):
            radius_sq = 0
            for i, mv_name in enumerate(self.problem_description.symbol_list['MV']):
                radius_sq += (x[i] - tr_base[mv_name]) ** 2 / (self.problem_description.scaling_factors[mv_name] ** 2)
            return max(radius_sq - tr_radius ** 2, 0)

        # solve tangential step optimization problem
        pso_normal_problem = PSO(func=tangential_step_problem_obj, n_dim=len(lb), pop=self.options['population'], \
                                 max_iter=self.options["max_iter"], lb=lb, ub=ub, \
                                 w=self.options["inertia"], c1=self.options["c1"], \
                                 c2=self.options["c2"], constraint_ueq=(tangential_step_problem_c_con,\
                                                                        tangential_step_problem_tr_con))
        pso_normal_problem.run()
        x_after_tangential_step = pso_normal_problem.gbest_x[0]
        inputs = {}
        for i, mv_name in enumerate(self.problem_description.symbol_list['MV']):
            inputs[mv_name] = x_after_tangential_step[i]

        solve_status = ModelSolvingStatus.OK
        return inputs, input_after_normal_step, solve_status

    def optimize(self, tr_radius, tr_base, xi_N, modifiers_value, objective_scaling=1):
        sigma = 1e6
        # set modifiers
        assert(isinstance(self.black_box_model, BlackBoxModelWithModifiers))
        self.black_box_model.set_modifiers(modifiers_value)

        # get mv bounds
        lb = []
        ub = []
        mv_scaling = []
        for i, mv_name in enumerate(self.problem_description.symbol_list['MV']):
            lb.append(-1)
            ub.append(1)
            mv_scaling.append(self.problem_description.scaling_factors[mv_name])

        def normal_step_problem_obj(x):
            x = square_circle_mapping(x, tr_radius*xi_N, mv_scaling)
            for i, mv_name in enumerate(self.problem_description.symbol_list['MV']):
                x[i]+=tr_base[mv_name]
            input_dict = {}
            for i,mv_name in enumerate(self.problem_description.symbol_list['MV']):
                input_dict[mv_name] = x[i]
            input_dict = self.adapt_to_bound_mv(input_dict, self.problem_description.bounds)
            output_dict = self.black_box_model.simulate(input_dict)
            infeasibility_sq = 0
            for con_var_name in self.problem_description.symbol_list['CON']:
                infeasibility_sq+=(max(0,output_dict[con_var_name])/self.problem_description.scaling_factors[con_var_name])**2
            return infeasibility_sq

        # solve normal step optimization problem
        pso_normal_problem = PSO(func=normal_step_problem_obj, n_dim=len(lb), pop=self.options['population'], \
                                 max_iter=self.options["max_iter"], lb=lb, ub=ub, \
                                 w=self.options["inertia"], c1=self.options["c1"], \
                                 c2=self.options["c2"])
        pso_normal_problem.run()
        gbest_x_after_normal_step = pso_normal_problem.gbest_x
        x_after_normal_step = square_circle_mapping(gbest_x_after_normal_step, tr_radius*xi_N, mv_scaling)
        for i, mv_name in enumerate(self.problem_description.symbol_list['MV']):
            x_after_normal_step[i] += tr_base[mv_name]
        tr_infeasibility_sq_after_normal_step = normal_step_problem_obj(gbest_x_after_normal_step)
        print("infeasibility after normal step: %.6e"%numpy.sqrt(tr_infeasibility_sq_after_normal_step))
        input_after_normal_step = {}
        for i, mv_name in enumerate(self.problem_description.symbol_list['MV']):
            input_after_normal_step[mv_name] = x_after_normal_step[i]
        input_after_normal_step = self.adapt_to_bound_mv(input_after_normal_step, self.problem_description.bounds)

        def tangential_step_problem_obj(x):
            x = square_circle_mapping(x, tr_radius, mv_scaling)
            for i, mv_name in enumerate(self.problem_description.symbol_list['MV']):
                x[i] += tr_base[mv_name]
            input_dict = {}
            for i, mv_name in enumerate(self.problem_description.symbol_list['MV']):
                input_dict[mv_name] = x[i]
            input_dict = self.adapt_to_bound_mv(input_dict, self.problem_description.bounds)
            output_dict = self.black_box_model.simulate(input_dict)
            infeasibility_sq = 0
            for con_var_name in self.problem_description.symbol_list['CON']:
                infeasibility_sq += (max(0,output_dict[con_var_name]) /
                         self.problem_description.scaling_factors[con_var_name]) ** 2
            infeasibility_penalty = max(numpy.sqrt(infeasibility_sq)-\
                                        numpy.sqrt(tr_infeasibility_sq_after_normal_step),0)
            return output_dict[self.problem_description.symbol_list['OBJ']]*objective_scaling+\
                sigma*infeasibility_penalty

        # solve tangential step optimization problem
        pso_tangential_problem = PSO(func=tangential_step_problem_obj, n_dim=len(lb), pop=self.options['population'], \
                                 max_iter=self.options["max_iter"], lb=lb, ub=ub, \
                                 w=self.options["inertia"], c1=self.options["c1"], \
                                 c2=self.options["c2"])
        pso_tangential_problem.run()
        gbest_after_tangential_step = pso_tangential_problem.gbest_x
        x_after_tangential_step = square_circle_mapping(gbest_after_tangential_step, tr_radius, mv_scaling)
        for i, mv_name in enumerate(self.problem_description.symbol_list['MV']):
            x_after_tangential_step[i] += tr_base[mv_name]
        inputs = {}
        for i, mv_name in enumerate(self.problem_description.symbol_list['MV']):
            inputs[mv_name] = x_after_tangential_step[i]
        inputs = self.adapt_to_bound_mv(inputs, self.problem_description.bounds)
        solve_status = ModelSolvingStatus.OK
        print(inputs)
        tr_infeasibility_sq_after_tangential_step = normal_step_problem_obj(gbest_after_tangential_step)
        print("infeasibility after tangential step: %.6e" % numpy.sqrt(tr_infeasibility_sq_after_tangential_step))
        return inputs, input_after_normal_step, solve_status