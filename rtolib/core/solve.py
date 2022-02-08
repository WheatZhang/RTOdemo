from enum import Enum
from rtolib.core.basic import Simulator, Optimizer, ParameterEstimator, ProblemDescription
from pyomo.environ import *
from pyomo.opt import SolverStatus, TerminationCondition
import rtolib.util.init_value as init_value
from rtolib.core.pyomo_model import PyomoModel, PyomoModelWithModifiers

__name__ = ['PyomoModelSolvingStatus','PyomoSimulator','PyomoOptimizer',
           'PyomoMultiDataPointParamEstimator',
            'PyomoGradientParamEstimator']

class PyomoModelSolvingStatus(Enum):
    OK = 1
    HOMOTOPY_TARGET_NOT_REACHED = 2
    OPTIMIZATION_FAILED = 3
    OTHER = 0


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
        results = self.solver.solve(self.model, tee=self.tee)
        if not ((results.solver.status == SolverStatus.ok) and (
                results.solver.termination_condition == TerminationCondition.optimal)):
            print(results.solver.termination_condition)
            # self.model.display(r"F:\Research\RTOdemo\debug\WO\result\display.txt")
            solve_status = PyomoModelSolvingStatus.OPTIMIZATION_FAILED
        else:
            if use_homo and value(self.model.homotopy_simulation_obj) >1e-4:
                solve_status = PyomoModelSolvingStatus.HOMOTOPY_TARGET_NOT_REACHED
            else:
                solve_status = PyomoModelSolvingStatus.OK
        outputs = {}
        for op in self.pyomo_model.output_variables.keys():
            var = self.pyomo_model.output_variables[op].__call__(self.model)
            outputs[op] = value(var)
        return outputs, solve_status


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


class PyomoMultiDataPointParamEstimator(ParameterEstimator):
    def __init__(self, pyomo_model, number_of_data_point, use_obj_and_con_info=False):
        '''
        The first data point is used to estimate output modifier.
        The other data points are used to estimate gradient modifiers.
        :param pyomo_model:
        :param number_of_data_point:
        '''
        assert isinstance(pyomo_model, PyomoModel)
        self.pyomo_model = pyomo_model
        self.number_of_data_point = number_of_data_point
        self.model = ConcreteModel()
        self.homotopy_var = []
        for var_name in self.pyomo_model.input_variables.keys():
            self.homotopy_var.append(var_name)
        self.use_obj_and_con_info = use_obj_and_con_info


    def set_base_point(self, base_point):
        if not isinstance(self.pyomo_model, PyomoModelWithModifiers):
            raise AttributeError('this model does not have base point')
        for case in self.model.CaseIndex:
            self.pyomo_model.set_base_point(self.model.Cases[case], base_point)

    def set_solver(self, solver, tee, default_options):
        self.solver = solver
        self.tee = tee
        for k, v in default_options.items():
            self.solver.options[k] = v

    def build(self, problem_description):
        assert isinstance(problem_description, ProblemDescription)
        self.problem_description = problem_description
        self.all_available_measurement = []
        for meas in self.problem_description.symbol_list['CV']:
            self.all_available_measurement.append(meas)
        if self.use_obj_and_con_info:
            self.all_available_measurement.append(self.problem_description.symbol_list['OBJ'])
            for meas in self.problem_description.symbol_list['CON']:
                self.all_available_measurement.append(meas)

        # set case index
        self.model.CaseIndex = RangeSet(1, self.number_of_data_point)

        # parameters for public use
        for pn in self.pyomo_model.parameters.keys():
            setattr(self.model, pn, Var(initialize=self.pyomo_model.default_value[pn]))

        # set each block
        def _block_rule(b, case):
            self.pyomo_model.build(b)

            # create variables for homotopy
            for ip in self.pyomo_model.input_variables.keys():
                setattr(b, ip + "_homo", Var(initialize=0))
                getattr(b, ip + "_homo").fixed = True
            for theta in self.pyomo_model.parameters.keys():
                setattr(b, theta + "_homo", Var(initialize=0))
                getattr(b, theta + "_homo").fixed = True

            def _homotopy_simulation_obj(m):
                ip_scaling_coeffs = self.problem_description.scaling_factors
                theta_scaling_coeffs = self.pyomo_model.parameter_scaling_factors
                return sum([(getattr(m, k + "_homo") - \
                             self.pyomo_model.input_variables[k].__call__(b)) ** 2 \
                            / ip_scaling_coeffs[k] / ip_scaling_coeffs[k] \
                            for k in self.pyomo_model.input_variables.keys()]) + \
                       sum([(getattr(m, k + "_homo") - \
                             self.pyomo_model.parameters[k].__call__(b)) ** 2 \
                            / theta_scaling_coeffs[k] / theta_scaling_coeffs[k] \
                            for k in self.pyomo_model.parameters.keys()])
            b.homotopy_simulation_obj = Objective(rule=_homotopy_simulation_obj)

        self.model.Cases = Block(self.model.CaseIndex, rule=_block_rule)

        # link constraints
        def _link_con_rule(m):
            for pn in self.pyomo_model.parameters.keys():
                scaling_factor = self.pyomo_model.parameter_scaling_factors[pn]
                for no in range(1,self.number_of_data_point+1):
                    yield (getattr(m, pn)-self.pyomo_model.parameters[pn].__call__(m.Cases[no]))/scaling_factor*100==0
        self.model.link_con_rule = ConstraintList(rule=_link_con_rule)

        for case in self.model.CaseIndex:
            # deactivate inequality constraints
            for c in self.model.Cases[case].component_data_objects(
                    Constraint, descend_into=True):
                if c.upper is None or c.lower is None:
                    c.deactivate()

            # load initial value file
            init_value.load_init_from_template(self.model.Cases[case], self.pyomo_model.initial_value_file, \
                                               ignore_init_mismatch=True)

        # TODO: how to scale parameter variable, like in WO
        # parameter estimation weight and expectation
        for pn in self.pyomo_model.parameters.keys():
            setattr(self.model, pn + "_weight", Var(initialize=0))
            getattr(self.model, pn + "_weight").fixed = True
            setattr(self.model, pn + "_exp", Var(initialize=0))
            getattr(self.model, pn + "_exp").fixed = True
        for op in self.all_available_measurement:
            setattr(self.model, op + "_weight", Var(initialize=0))
            getattr(self.model, op + "_weight").fixed = True

        # output measurement
        for op in self.all_available_measurement:
            setattr(self.model, op+'_measure', Var(self.model.CaseIndex, initialize=0))

        # parameter estimation objective
        def _pe_obj(m):
            return sum([getattr(self.model, k + "_weight") * ( \
                        sum([(self.pyomo_model.output_variables[k].__call__(self.model.Cases[c]) - \
                         getattr(self.model, k+'_measure')[c]) ** 2 \
                             for c in self.model.CaseIndex]) \
                ) for k in self.all_available_measurement]) + \
                   sum([getattr(self.model, k + "_weight") * ( \
                               (getattr(self.model, k) - getattr(self.model,k + "_exp")) ** 2 \
                       ) for k in self.pyomo_model.parameters.keys()])
        self.model.pe_obj = Objective(rule=_pe_obj)

    def set_measure_data(self, plant_data):
        '''

        :param plant_data: list of dicts
        :return:
        '''
        for i, tp in enumerate(plant_data):
            for k in self.all_available_measurement:
                var = getattr(self.model, k+"_measure")[i+1]
                var.fix(tp[k])

    def estimate_parameter(self, plant_data, fixed_param_values={}, use_homo=True, pre_simulation=True):
        # set measurement
        if len(plant_data) != self.number_of_data_point:
            raise ValueError("plant data length does not match")
        self.set_measure_data(plant_data)

        # pre simulation
        if pre_simulation:
            for case_no in self.model.CaseIndex:
                block = self.model.Cases[case_no]
                # autocomplete
                param_values = {}
                for para_name in self.pyomo_model.parameters.keys():
                    if para_name not in fixed_param_values.keys():
                        param_values[para_name] = self.pyomo_model.default_value[para_name]
                input_values = {}
                for var_name in self.pyomo_model.input_variables.keys():
                    if var_name not in plant_data[case_no-1].keys():
                        input_values[var_name] = self.pyomo_model.default_value[var_name]
                    else:
                        input_values[var_name] = plant_data[case_no-1][var_name]

                # Fix inputs and parameters
                if not use_homo:
                    for var_name in self.pyomo_model.input_variables.keys():
                        var = self.pyomo_model.input_variables[var_name].__call__(block)
                        var.fix(input_values[var_name])
                    for var_name in self.pyomo_model.parameters.keys():
                        var = self.pyomo_model.parameters[var_name].__call__(block)
                        var.fix(param_values[var_name])
                else:
                    for var_name in self.pyomo_model.input_variables.keys():
                        var = self.pyomo_model.input_variables[var_name].__call__(block)
                        if var_name in self.homotopy_var:
                            var.fixed = False
                            getattr(block, var_name + "_homo").fix(input_values[var_name])
                        else:
                            var.fix(input_values[var_name])
                            getattr(block, var_name + "_homo").fix(input_values[var_name])
                    for var_name in self.pyomo_model.parameters.keys():
                        var = self.pyomo_model.parameters[var_name].__call__(block)
                        if var_name in self.homotopy_var:
                            var.fixed = False
                            getattr(block, var_name + "_homo").fix(param_values[var_name])
                        else:
                            var.fix(param_values[var_name])
                            getattr(block, var_name + "_homo").fix(param_values[var_name])
                # homotopy simulation
                if use_homo:
                    block.homotopy_simulation_obj.activate()
                else:
                    block.homotopy_simulation_obj.deactivate()
                # solve the problem
                # TODO: tuning options
                self.solver.options["obj_scaling_factor"]=1
                results = self.solver.solve(block, tee=self.tee)
                if not ((results.solver.status == SolverStatus.ok) and (
                        results.solver.termination_condition == TerminationCondition.optimal)):
                    print(results.solver.termination_condition)
                    solve_status = PyomoModelSolvingStatus.OPTIMIZATION_FAILED
                else:
                    if use_homo and value(block.homotopy_simulation_obj) > 1e-4:
                        solve_status = PyomoModelSolvingStatus.HOMOTOPY_TARGET_NOT_REACHED
                    else:
                        solve_status = PyomoModelSolvingStatus.OK
                #TODO: deal with solve_status
                print(solve_status)

        # prepare pe
        for case_no in self.model.CaseIndex:
            block = self.model.Cases[case_no]
            block.homotopy_simulation_obj.deactivate()
            for para_name in self.pyomo_model.parameters.keys():
                self.pyomo_model.parameters[para_name].__call__(block).fixed=False

        # set active parameters
        for para_name in self.pyomo_model.parameters.keys():
            init_value = value(self.pyomo_model.parameters[para_name].__call__(\
                    self.model.Cases[1]))
            if para_name in fixed_param_values.keys():
                getattr(self.model, para_name).fix(init_value)
            else:
                getattr(self.model, para_name)._value = init_value
                getattr(self.model, para_name).fixed=False

        # fix input variable
        for case_no in self.model.CaseIndex:
            block = self.model.Cases[case_no]
            for var_name in self.pyomo_model.input_variables.keys():
                if var_name not in plant_data[case_no - 1].keys():
                    self.pyomo_model.input_variables[var_name].__call__(block).fix(\
                        self.pyomo_model.default_value[var_name])
                else:
                    self.pyomo_model.input_variables[var_name].__call__(block).fix(\
                        plant_data[case_no - 1][var_name])

        # TODO: back up model data

        # solve pe problem
        # if obj_scaling_factor is small, the chance of failure rises
        self.solver.options["obj_scaling_factor"] = 1e3
        results = self.solver.solve(self.model, tee=self.tee)
        # TODO: back up model data
        # TODO: fallback strategy

        # return result
        para_ret={}
        for para_name in self.pyomo_model.parameters.keys():
            para_ret[para_name] = value(getattr(self.model,para_name))
        para_ret["$residue"]=value(self.model.pe_obj)
        return para_ret

    def set_weight(self, output_weight, parameter_weight):
        for k,v in output_weight.items():
            getattr(self.model, k + "_weight").fix(v)
        for k,v in parameter_weight.items():
            getattr(self.model, k + "_weight").fix(v)

    def set_parameter_guess(self, parameter_guess={}):
        for para_name in self.pyomo_model.parameters.keys():
            if para_name in parameter_guess.keys():
                getattr(self.model, para_name + "_exp").fix(parameter_guess[para_name])
            else:
                getattr(self.model, para_name + "_exp").fix(self.pyomo_model.default_value[para_name])


class PyomoGradientParamEstimator(PyomoMultiDataPointParamEstimator):
    def build(self, problem_description):
        assert isinstance(problem_description, ProblemDescription)
        self.problem_description = problem_description

        # set case index
        self.model.CaseIndex = RangeSet(1, self.number_of_data_point)
        self.model.CaseIndexForGradient = RangeSet(2, self.number_of_data_point)

        # parameters for public use
        for pn in self.pyomo_model.parameters.keys():
            setattr(self.model, pn, Var(initialize=self.pyomo_model.default_value[pn]))

        # set each block
        def _block_rule(b, case):
            self.pyomo_model.build(b)

            # create variables for homotopy
            for ip in self.pyomo_model.input_variables.keys():
                setattr(b, ip + "_homo", Var(initialize=0))
                getattr(b, ip + "_homo").fixed = True
            for theta in self.pyomo_model.parameters.keys():
                setattr(b, theta + "_homo", Var(initialize=0))
                getattr(b, theta + "_homo").fixed = True

            def _homotopy_simulation_obj(m):
                ip_scaling_coeffs = self.problem_description.scaling_factors
                theta_scaling_coeffs = self.pyomo_model.parameter_scaling_factors
                return sum([(getattr(m, k + "_homo") - \
                             self.pyomo_model.input_variables[k].__call__(b)) ** 2 \
                            / ip_scaling_coeffs[k] / ip_scaling_coeffs[k] \
                            for k in self.pyomo_model.input_variables.keys()]) + \
                       sum([(getattr(m, k + "_homo") - \
                             self.pyomo_model.parameters[k].__call__(b)) ** 2 \
                            / theta_scaling_coeffs[k] / theta_scaling_coeffs[k] \
                            for k in self.pyomo_model.parameters.keys()])

            b.homotopy_simulation_obj = Objective(rule=_homotopy_simulation_obj)

        self.model.Cases = Block(self.model.CaseIndex, rule=_block_rule)

        # link constraints
        def _link_con_rule(m):
            for pn in self.pyomo_model.parameters.keys():
                scaling_factor = self.pyomo_model.parameter_scaling_factors[pn]
                for no in range(1, self.number_of_data_point + 1):
                    yield (getattr(m, pn) - self.pyomo_model.parameters[pn].__call__(
                        m.Cases[no])) / scaling_factor * 100 == 0

        self.model.link_con_rule = ConstraintList(rule=_link_con_rule)

        for case in self.model.CaseIndex:
            # deactivate inequality constraints
            for c in self.model.Cases[case].component_data_objects(
                    Constraint, descend_into=True):
                if c.upper is None or c.lower is None:
                    c.deactivate()

            # load initial value file
            init_value.load_init_from_template(self.model.Cases[case], self.pyomo_model.initial_value_file, \
                                               ignore_init_mismatch=True)

        # TODO: how to scale parameter variable, like in WO
        # parameter estimation weight and expectation
        for pn in self.pyomo_model.parameters.keys():
            setattr(self.model, pn + "_weight", Var(initialize=0))
            getattr(self.model, pn + "_weight").fixed = True
            setattr(self.model, pn + "_exp", Var(initialize=0))
            getattr(self.model, pn + "_exp").fixed = True
        for op in self.all_available_measurement:
            setattr(self.model, op + "_weight", Var(initialize=0))
            getattr(self.model, op + "_weight").fixed = True

        # output measurement
        for op in self.all_available_measurement:
            setattr(self.model, op + '_measure', Var(self.model.CaseIndex, initialize=0))

        # parameter estimation objective
        def _pe_obj(m):
            return sum([getattr(self.model, k + "_weight") * ( \
                        (self.pyomo_model.output_variables[k].__call__(self.model.Cases[1]) - \
                         getattr(self.model, k + '_measure')[1]) ** 2 + \
                        sum([(self.pyomo_model.output_variables[k].__call__(self.model.Cases[c]) - \
                              self.pyomo_model.output_variables[k].__call__(self.model.Cases[1]) - \
                              getattr(self.model, k + '_measure')[c] + \
                              getattr(self.model, k + '_measure')[1]) ** 2 / 2 \
                             for c in self.model.CaseIndexForGradient]) \
                ) for k in self.all_available_measurement]) + \
                   sum([getattr(self.model, k + "_weight") * ( \
                               (getattr(self.model, k) - getattr(self.model, k + "_exp")) ** 2 \
                       ) for k in self.pyomo_model.parameters.keys()])

        self.model.pe_obj = Objective(rule=_pe_obj)

    def set_measure_data(self, plant_data):
        '''

        :param plant_data: list of dicts
        :return:
        '''
        for i,tp in enumerate(plant_data):
            if i != 0:
                flag = False
                for k in self.problem_description.symbol_list['MV']:
                    if abs(plant_data[0][k] - plant_data[i][k])/self.problem_description.scaling_factors[k] > 1e-6:
                        flag = True
                        break
                if not flag:
                    raise ValueError("Perturbation stepsize too small")
        super().set_measure_data(plant_data)