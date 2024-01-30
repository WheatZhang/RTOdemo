from pyomo.environ import ConcreteModel, Var, Objective, Constraint, Expression,\
    TerminationCondition, value, SolverStatus
from rtolib.core.pyomo_model import PyomoModel, PyomoModelWithModifiers,\
    ModelSolvingStatus
from rtolib.core.basic import ProblemDescription
import rtolib.util.init_value as init_value


class Simulator():
    def simulate(self, input_values, param_values=None):
        return NotImplementedError()

class PyomoSimulator(Simulator):
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

    def build(self, problem_description):
        assert isinstance(problem_description, ProblemDescription)
        self.problem_description = problem_description
        self.pyomo_model.build(self.model)

        # create variables for homotopy
        for ip in self.pyomo_model.input_variables.keys():
            setattr(self.model, ip + "_homo", Var(initialize=0))
            getattr(self.model, ip + "_homo").fixed = True
        for theta in self.pyomo_model.parameters.keys():
            setattr(self.model, theta + "_homo", Var(initialize=0))
            getattr(self.model, theta + "_homo").fixed = True
        def _homotopy_simulation_obj(m):
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
        self.model.homotopy_simulation_obj = Objective(rule=_homotopy_simulation_obj)
        # deactivate inequality constraints
        for c in self.model.component_data_objects(
                Constraint, descend_into=True):
            if c.upper is None or c.lower is None:
                if c.name not in self.pyomo_model.model_inequality_constraints:
                    c.deactivate()
        # load initial value file
        init_value.load_init_from_template(self.model, self.pyomo_model.initial_value_file, \
                                                       ignore_init_mismatch=True)

    def simulate(self, input_values, param_values=None, use_homo=True):
        # autocomplete
        if param_values is None:
            param_values = {}
            for para_name in self.pyomo_model.parameters.keys():
                param_values[para_name] = self.pyomo_model.default_value[para_name]
        else:
            for para_name in self.pyomo_model.parameters.keys():
                if para_name not in param_values.keys():
                    param_values[para_name] = self.pyomo_model.default_value[para_name]
        for var_name in self.pyomo_model.input_variables.keys():
            if var_name not in input_values.keys():
                input_values[var_name] = self.pyomo_model.default_value[var_name]
        # Fix inputs and parameters
        if not use_homo:
            for var_name in self.pyomo_model.input_variables.keys():
                var = self.pyomo_model.input_variables[var_name].__call__(self.model)
                var.fix(input_values[var_name])
            for var_name in self.pyomo_model.parameters.keys():
                var = self.pyomo_model.parameters[var_name].__call__(self.model)
                var.fix(param_values[var_name])
        else:
            for var_name in self.pyomo_model.input_variables.keys():
                var = self.pyomo_model.input_variables[var_name].__call__(self.model)
                if var_name in self.homotopy_var:
                    var.fixed = False
                    getattr(self.model, var_name + "_homo").fix(input_values[var_name])
                else:
                    var.fix(input_values[var_name])
                    getattr(self.model, var_name + "_homo").fix(input_values[var_name])
            for var_name in self.pyomo_model.parameters.keys():
                var = self.pyomo_model.parameters[var_name].__call__(self.model)
                if var_name in self.homotopy_var:
                    var.fixed = False
                    getattr(self.model, var_name + "_homo").fix(param_values[var_name])
                else:
                    var.fix(param_values[var_name])
                    getattr(self.model, var_name + "_homo").fix(param_values[var_name])

        # homotopy simulation
        if use_homo:
            self.model.homotopy_simulation_obj.activate()
        else:
            self.model.homotopy_simulation_obj.deactivate()
        # solve the problem
        if self.solver.name == 'ipopt':
            self.solver.options['tol'] = 1e-10
        # results = self.solver.solve(self.model, tee=True)
        results = self.solver.solve(self.model, tee=self.tee)
        if not ((results.solver.status == SolverStatus.ok) and (
                results.solver.termination_condition == TerminationCondition.optimal)):
            print(results.solver.termination_condition)
            # self.model.display(r"F:\Research\RTOdemo\debug\WO\result\display.txt")
            solve_status = ModelSolvingStatus.OPTIMIZATION_FAILED
        else:
            if use_homo and value(self.model.homotopy_simulation_obj) >1e-4:
                solve_status = ModelSolvingStatus.HOMOTOPY_TARGET_NOT_REACHED
            else:
                solve_status = ModelSolvingStatus.OK
        outputs = {}
        for op in self.pyomo_model.output_variables.keys():
            var = self.pyomo_model.output_variables[op].__call__(self.model)
            outputs[op] = value(var)
        return outputs, solve_status