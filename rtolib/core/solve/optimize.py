from pyomo.environ import ConcreteModel, Var, Objective, Expression,\
    TerminationCondition, value, SolverStatus, minimize
from rtolib.core.pyomo_model import PyomoModel, PyomoModelWithModifiers,\
    PyomoModelSolvingStatus
from rtolib.core.basic import ProblemDescription
import rtolib.util.init_value as init_value


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
