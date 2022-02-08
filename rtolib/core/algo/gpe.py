from rtolib.core import Algorithm
from pyomo.environ import *
import rtolib.util.init_value as InitValueToolkit
from rtolib.util.misc import is_value_dict_same
import os
import numpy
import copy

class GPE_Algorithm(Algorithm):
    def set_weight(self, output_weight, param_weight):
        '''

        :param output_weight: objective function认为是一个特殊的output
        :param param_weight:
        :return:
        '''
        self.output_weight = output_weight
        self.param_weight = param_weight
        self.active_output_variables = self.output_weight.keys()
        self.active_parameters = self.param_weight.keys()

    def update_model_basepoint(self, basepoint):
        for k,v in basepoint.items():
            getattr(self.model_simulator,k+"_basepoint").fix(v)
            for case_no in self.model_pe_optimizer.CaseIndex:
                getattr(self.model_pe_optimizer.Cases[case_no],k+"_basepoint").fix(v)

    def set_simulation_model(self, rto_model, parameter_init_values, basepoint):
        '''

        :param rto_model:
        :param parameter_init_values: A dict
        :param basepoint: A dict
        :return:
        '''
        self.parameter_init_values = parameter_init_values
        self.rto_model = rto_model
        self.parameters = rto_model.parameters
        self.model_simulator = ConcreteModel()
        for pn in self.parameters:
            setattr(self.model_simulator, pn , Var(initialize=parameter_init_values[pn]))
            getattr(self.model_simulator, pn).fixed = True
        parameter = {}
        for pn in self.parameters:
            parameter[pn] = getattr(self.model_simulator, pn)
        rto_model.build(self.model_simulator, parameter, basepoint)
        if not hasattr(self, "mv_scaling_coeff"):
            self.mv_scaling_coeff = {}
            for mv in self.manipulated_variables:
                var = self.rto_model.manipulated_variables[mv].__call__(self.model_simulator)
                if self.mv_bounds[mv][0] is not None and self.mv_bounds[mv][1] is not None:
                    self.mv_scaling_coeff[mv] = self.mv_bounds[mv][1] - self.mv_bounds[mv][0]
                elif abs(value(var)) > 0.1:
                    self.mv_scaling_coeff[mv] = value(var)
                else:
                    self.mv_scaling_coeff[mv] = 1
        if not hasattr(self, "sv_scaling_coeff"):
            self.sv_scaling_coeff = {}
            for sv in self.specifications:
                var = self.rto_model.specifications[sv].__call__(self.model_simulator)
                if self.sv_bounds[sv][0] is not None and self.sv_bounds[sv][1] is not None:
                    self.sv_scaling_coeff[sv] = self.sv_bounds[sv][1] - self.sv_bounds[sv][0]
                elif abs(value(var)) > 0.1:
                    self.sv_scaling_coeff[sv] = value(var)
                else:
                    self.sv_scaling_coeff[sv] = 1

        for mv in self.manipulated_variables:
            setattr(self.model_simulator, mv + "_prosim", Var(initialize=0))
            getattr(self.model_simulator, mv + "_prosim").fixed = True
        for sv in self.specifications:
            setattr(self.model_simulator, sv + "_prosim", Var(initialize=0))
            getattr(self.model_simulator, sv + "_prosim").fixed = True

        def _simulation_promoting_obj(m):
            return sum([(getattr(m, k + "_prosim") - \
                         self.rto_model.manipulated_variables[k].__call__(self.model_simulator)) ** 2 \
                        / self.mv_scaling_coeff[k] / self.mv_scaling_coeff[k] \
                        for k in self.manipulated_variables])+\
                    sum([(getattr(m, k + "_prosim") - \
                          self.rto_model.specifications[k].__call__(self.model_simulator)) ** 2 \
                         / self.sv_scaling_coeff[k] / self.sv_scaling_coeff[k] \
                         for k in self.specifications])
        self.model_simulator.simulation_promoting_obj = Objective(rule=_simulation_promoting_obj)
        self.model_simulator.simulation_promoting_obj.deactivate()

    def set_model(self, rto_model, parameter_init_values, basepoint, parameter_bounds=None):
        '''

        :param rto_model:
        :param parameter_init_values: A dict
        :param basepoint: A dict
        :param parameter_bounds: A dict
        :return:
        '''
        self.parameter_init_values = parameter_init_values
        self.parameter_bounds = parameter_bounds
        self.rto_model = rto_model
        if not hasattr(self, "rto_plant"):
            raise Exception("Please set rto_plant before setting rto_model")
        for k in self.rto_plant.manipulated_variables.keys():
            if k not in self.rto_model.manipulated_variables.keys():
                raise KeyError("MV %s in rto_plant not available for rto_model")
        for k in self.rto_model.manipulated_variables.keys():
            if k not in self.rto_plant.manipulated_variables.keys():
                raise KeyError("MV %s in rto_model not available for rto_plant")
        for k in self.rto_plant.output_variables.keys():
            if k not in self.rto_model.output_variables.keys():
                raise KeyError("Output variable %s in rto_plant not available for rto_model")
        for k in self.rto_model.output_variables.keys():
            if k not in self.rto_plant.output_variables.keys():
                raise KeyError("Output variable %s in rto_model not available for rto_plant")
        self.parameters = rto_model.parameters
        self.model_simulator = ConcreteModel()
        for pn in self.parameters:
            if parameter_bounds is None:
                setattr(self.model_simulator, pn , Var(initialize=self.parameter_init_values[pn]))
            else:
                setattr(self.model_simulator, pn, Var(initialize=self.parameter_init_values[pn],bounds=self.parameter_bounds[pn]))
            getattr(self.model_simulator, pn).fixed = True
        parameter = {}
        for pn in self.parameters:
            parameter[pn] = getattr(self.model_simulator, pn)
        rto_model.build(self.model_simulator, parameter, basepoint)
        if not hasattr(self, "mv_scaling_coeff"):
            self.mv_scaling_coeff = {}
            for mv in self.manipulated_variables:
                var = self.rto_model.manipulated_variables[mv].__call__(self.model_simulator)
                if self.mv_bounds[mv][0] is not None and self.mv_bounds[mv][1] is not None:
                    self.mv_scaling_coeff[mv] = self.mv_bounds[mv][1] - self.mv_bounds[mv][0]
                elif abs(value(var)) > 0.1:
                    self.mv_scaling_coeff[mv] = value(var)
                else:
                    self.mv_scaling_coeff[mv] = 1
        if not hasattr(self, "sv_scaling_coeff"):
            self.sv_scaling_coeff = {}
            for sv in self.specifications:
                var = self.rto_model.specifications[sv].__call__(self.model_simulator)
                if self.sv_bounds[sv][0] is not None and self.sv_bounds[sv][1] is not None:
                    self.sv_scaling_coeff[sv] = self.sv_bounds[sv][1] - self.sv_bounds[sv][0]
                elif abs(value(var)) > 0.1:
                    self.sv_scaling_coeff[sv] = value(var)
                else:
                    self.sv_scaling_coeff[sv] = 1

        for mv in self.manipulated_variables:
            setattr(self.model_simulator, mv + "_prosim", Var(initialize=0))
            getattr(self.model_simulator, mv + "_prosim").fixed = True
        for sv in self.specifications:
            setattr(self.model_simulator, sv + "_prosim", Var(initialize=0))
            getattr(self.model_simulator, sv + "_prosim").fixed = True

        # prosim = simulation promotion technology
        def _simulation_promoting_obj(m):
            return sum([(getattr(m, k + "_prosim") - \
                         self.rto_model.manipulated_variables[k].__call__(self.model_simulator)) ** 2 \
                        / self.mv_scaling_coeff[k] / self.mv_scaling_coeff[k] \
                        for k in self.manipulated_variables])+\
                    sum([(getattr(m, k + "_prosim") - \
                          self.rto_model.specifications[k].__call__(self.model_simulator)) ** 2 \
                         / self.sv_scaling_coeff[k] / self.sv_scaling_coeff[k] \
                         for k in self.specifications])
        self.model_simulator.simulation_promoting_obj = Objective(rule=_simulation_promoting_obj)
        self.model_simulator.simulation_promoting_obj.deactivate()

        def _opt_promoted_obj(m):
            return self.rto_model.output_variables['Phi'].__call__(m) + \
                   self.obj_scaling_factor['Opt_promoted_multiplier'] * ( \
                               sum([(getattr(m, k + "_prosim") - \
                                     self.rto_model.specifications[k].__call__(m)) ** 2 \
                                    / self.sv_scaling_coeff[k] / self.sv_scaling_coeff[k] \
                                    for k in self.specifications]))

        self.model_simulator.opt_promoted_obj = Objective(rule=_opt_promoted_obj)
        self.model_simulator.opt_promoted_obj.deactivate()

        number_of_u = len(self.manipulated_variables)
        self.model_pe_optimizer = ConcreteModel()
        self.model_pe_optimizer.CaseIndex = RangeSet(1, number_of_u+1)
        self.model_pe_optimizer.CaseIndexForGradient = RangeSet(2, number_of_u + 1)
        for pn in self.parameters:
            if self.parameter_bounds is None:
                setattr(self.model_pe_optimizer, pn, Var(initialize=self.parameter_init_values[pn]))
            else:
                setattr(self.model_pe_optimizer, pn , Var(initialize=self.parameter_init_values[pn],bounds=self.parameter_bounds[pn]))
            if pn not in self.active_parameters:
                getattr(self.model_pe_optimizer, pn).fixed = True

        parameter = {}
        for pn in self.parameters:
            parameter[pn] = getattr(self.model_pe_optimizer, pn)
        def _block_rule(b, case):
            rto_model.build(b, parameter, basepoint)

            for mv in self.manipulated_variables:
                setattr(b, mv + "_prosim", Var(initialize=0))
                getattr(b, mv + "_prosim").fixed = True
            for sv in self.specifications:
                setattr(b, sv + "_prosim", Var(initialize=0))
                getattr(b, sv + "_prosim").fixed = True
            def _simulation_promoting_obj(m):
                return sum([(getattr(m, k + "_prosim") - \
                             self.rto_model.manipulated_variables[k].__call__(b)) ** 2 \
                            / self.mv_scaling_coeff[k] / self.mv_scaling_coeff[k] \
                            for k in self.manipulated_variables])+ \
                       sum([(getattr(m, k + "_prosim") - \
                             self.rto_model.specifications[k].__call__(b)) ** 2 \
                            / self.sv_scaling_coeff[k] / self.sv_scaling_coeff[k] \
                            for k in self.specifications])
            b.simulation_promoting_obj = Objective(rule=_simulation_promoting_obj)
            b.simulation_promoting_obj.deactivate()

        self.model_pe_optimizer.Cases = Block(self.model_pe_optimizer.CaseIndex, rule=_block_rule)
        for c in self.model_pe_optimizer.component_data_objects(
                Constraint, descend_into=True):
            if c.upper is None or c.lower is None:
                c.deactivate()

        for case in self.model_pe_optimizer.CaseIndex:
            obj = self.rto_model.objective.__call__(self.model_pe_optimizer.Cases[case])
            obj.deactivate()

        for pn in self.active_parameters:
            setattr(self.model_pe_optimizer, pn+"_weight" , Var(initialize=self.param_weight[pn]))
            getattr(self.model_pe_optimizer, pn+"_weight").fixed = True
            setattr(self.model_pe_optimizer, pn + "_expectation", Var(initialize=self.parameter_init_values[pn]))
            getattr(self.model_pe_optimizer, pn + "_expectation").fixed = True
        for op in self.active_output_variables:
            setattr(self.model_pe_optimizer, op+"_weight" , Var(initialize=self.output_weight[op]))
            getattr(self.model_pe_optimizer, op+"_weight").fixed = True
        def _pe_obj(m):
            return sum([getattr(self.model_pe_optimizer, k+"_weight")*(\
                (self.rto_model.output_variables[k].__call__(self.model_pe_optimizer.Cases[1])- \
                    self.rto_model.output_measurements[k].__call__(self.model_pe_optimizer.Cases[1]))**2+\
                sum([(self.rto_model.output_variables[k].__call__(self.model_pe_optimizer.Cases[c])- \
                    self.rto_model.output_variables[k].__call__(self.model_pe_optimizer.Cases[1])- \
                      self.rto_model.output_measurements[k].__call__(self.model_pe_optimizer.Cases[c]) + \
                       self.rto_model.output_measurements[k].__call__(self.model_pe_optimizer.Cases[1]))**2/2\
                    for c in self.model_pe_optimizer.CaseIndexForGradient])\
                ) for k in self.active_output_variables])+\
                   sum([getattr(self.model_pe_optimizer, k+"_weight")*( \
                       (getattr(self.model_pe_optimizer, k)-getattr(self.model_pe_optimizer, k+"_expectation"))**2\
                       ) for k in self.active_parameters])
        self.model_pe_optimizer.pe_obj = Objective(rule=_pe_obj)
        self.update_model_basepoint(basepoint)


    def set_model_init_value(self, init_value_filename):
        self.model_init_value_filename = init_value_filename
        InitValueToolkit.load_init_from_template(self.model_simulator, self.model_init_value_filename, \
                                                       ignore_init_mismatch=True)
        for case in self.model_pe_optimizer.CaseIndex:
            InitValueToolkit.load_init_from_template(self.model_pe_optimizer.Cases[case], self.model_init_value_filename, \
                                                           ignore_init_mismatch=True)


    def get_model_simulation(self, input, tee=False, simulation_promotion = True, print_level=2):
        for c in self.model_simulator.component_data_objects(
                Constraint, descend_into=True):
            if c.upper is None or c.lower is None:
                c.deactivate()
        for k in input.keys():
            var = self.rto_model.manipulated_variables[k].__call__(self.model_simulator)
            if simulation_promotion:
                getattr(self.model_simulator, k + "_prosim").fix(input[k])
                var.fixed = False
            else:
                var._value = input[k]
        for k in self.specifications:
            var = self.rto_model.specifications[k].__call__(self.model_simulator)
            if simulation_promotion:
                var.fixed = False
                getattr(self.model_simulator, k + "_prosim").fix(self.specification_values[k])
        if simulation_promotion:
            self.model_simulator.simulation_promoting_obj.activate()
            obj = self.rto_model.objective.__call__(self.model_simulator)
            obj.deactivate()

        if print_level >= 2:
            print("Solving Model Simulation")
        if simulation_promotion:
            self.solver.options['obj_scaling_factor'] = self.obj_scaling_factor['Simulation_promoted']
        else:
            self.solver.options['obj_scaling_factor'] = self.obj_scaling_factor['Simulation']
        results = self.solver.solve(self.model_simulator, tee=tee)
        self.check_solver_status(results)

        output_data = {}
        for output_key in self.output_variables:
            var = self.rto_model.output_variables[output_key].__call__(self.model_simulator)
            output_data[output_key] = value(var)
        for c in self.model_simulator.component_data_objects(
                Constraint, descend_into=True):
            if c.upper is None or c.lower is None:
                c.activate()
        if simulation_promotion:
            self.model_simulator.simulation_promoting_obj.deactivate()
            obj = self.rto_model.objective.__call__(self.model_simulator)
            obj.activate()
            for k in input.keys():
                var = self.rto_model.manipulated_variables[k].__call__(self.model_simulator)
                var.fixed = True
            for k in self.specifications:
                var = self.rto_model.specifications[k].__call__(self.model_simulator)
                var.fixed = True
        return output_data


    def set_measure_for_MultiPE(self, plant_data):
        for i,tp in enumerate(plant_data):
            for k in self.rto_model.output_measurements.keys():
                var = self.rto_model.output_measurements[k].__call__(self.model_pe_optimizer.Cases[i+1])
                var._value = tp[k]
            for k in self.rto_model.manipulated_variables.keys():
                var = self.rto_model.manipulated_variables[k].__call__(self.model_pe_optimizer.Cases[i+1])
                var._value = tp[k]
            if i != 0:
                flag = False
                for k in self.rto_model.manipulated_variables.keys():
                    if abs(plant_data[0][k] - plant_data[i][k]) > 1e-6:
                        flag = True
                        break
                if not flag:
                    raise Exception("Cannot find a reasonable FFD perturbation")

    def set_approximated_optimum(self, func):
        self.approximated_optimum_func = func


    def optimize_for_u(self, tee=False, print_level=2):
        # prepare init value
        # init_from_pe = {}
        # for v in InitValueToolkit.unfixed_variables_generator(self.model_pe_optimizer.Cases[1]):
        #     # print(v)
        #     if str(v).startswith("Cases[1]."):
        #         init_from_pe[str(v)[9:]] = value(v)
        # for v in InitValueToolkit.unfixed_variables_generator(self.model_simulator):
        #     if str(v) in init_from_pe.keys():
        #         v._value = init_from_pe[str(v)]
        # for c in self.model_simulator.component_data_objects(
        #         Constraint, descend_into=True):
        #     if c.upper is None or c.lower is None:
        #         c.deactivate()

        last_mv = {}
        for k in self.manipulated_variables:
            var = self.rto_model.manipulated_variables[k].__call__(self.model_simulator)
            last_mv[k] = value(var)
        approximated_optimum=self.approximated_optimum_func.__call__(last_mv, self.last_specification_values, \
                                                                     self.specification_values)
        # Using the approximated optimum as the init value
        # for k in self.manipulated_variables:
        #     var = self.rto_model.manipulated_variables[k].__call__(self.model_simulator)
        #     getattr(self.model_simulator, k + "_prosim").fix(approximated_optimum[k])
        #     var.fixed = False
        # for k in self.specifications:
        #     var = self.rto_model.specifications[k].__call__(self.model_simulator)
        #     var.fixed = False
        #     getattr(self.model_simulator, k + "_prosim").fix(self.specification_values[k])
        # self.model_simulator.simulation_promoting_obj.activate()
        # obj = self.rto_model.objective.__call__(self.model_simulator)
        # obj.deactivate()
        #
        # if print_level >= 2:
        #     print("Preparing initial value for optimization")
        # self.solver.options['obj_scaling_factor'] = self.obj_scaling_factor['Simulation_promoted']
        # results = self.solver.solve(self.model_simulator, tee=tee)
        # self.check_solver_status(results)
        #
        # for c in self.model_simulator.component_data_objects(
        #         Constraint, descend_into=True):
        #     if c.upper is None or c.lower is None:
        #         c.activate()
        # self.model_simulator.simulation_promoting_obj.deactivate()
        # obj = self.rto_model.objective.__call__(self.model_simulator)
        # obj.activate()
        # for k in self.manipulated_variables:
        #     var = self.rto_model.manipulated_variables[k].__call__(self.model_simulator)
        #     var.fixed = True
        # for k in self.specifications:
        #     var = self.rto_model.specifications[k].__call__(self.model_simulator)
        #     var.fixed = True

        # do homo optimization
        for k in self.specifications:
            var = self.rto_model.specifications[k].__call__(self.model_simulator)
            var.fixed = False
            getattr(self.model_simulator, k + "_prosim").fix(self.specification_values[k])
        self.model_simulator.opt_promoted_obj.activate()
        obj = self.rto_model.objective.__call__(self.model_simulator)
        obj.deactivate()

        for k in self.manipulated_variables:
            var = self.rto_model.manipulated_variables[k].__call__(self.model_simulator)
            var.fixed = False
        if print_level >= 2:
            print('Preparing initial value by homotopy optimization')
        self.solver.options['obj_scaling_factor'] = self.obj_scaling_factor['Optimization']
        results = self.solver.solve(self.model_simulator, tee=tee)
        self.check_solver_status(results)

        self.model_simulator.opt_promoted_obj.deactivate()
        obj = self.rto_model.objective.__call__(self.model_simulator)
        obj.activate()
        if not os.path.exists("init/debug_optimization/"):
            os.makedirs("init/debug_optimization/")
        InitValueToolkit.to_template(self.model_simulator, 'init/debug_optimization/homo_optimization%d.txt'%self._current_iter)

        # do optimization
        for k in self.manipulated_variables:
            var = self.rto_model.manipulated_variables[k].__call__(self.model_simulator)
            var.fixed = False
        for k in self.specifications:
            var = self.rto_model.specifications[k].__call__(self.model_simulator)
            var.fix(self.specification_values[k])
        if print_level >=2:
            print('Solving Optimization')
        self.solver.options['obj_scaling_factor'] = self.obj_scaling_factor['Optimization']
        results=None
        try:
            results = self.solver.solve(self.model_simulator, tee=tee)
            # InitValueToolkit.to_template(self.model_simulator, 'init/debug_optimization/%d.txt'%self._current_iter)
        except Exception as e:
            print(e)
            print("SSO problem restored.")
            InitValueToolkit.load_init_from_template(self.model_simulator, self.model_init_value_filename, \
                                                           ignore_init_mismatch=True)
        if results is not None:
            if not self.check_solver_status(results):
                print("SSO problem restored.")
                InitValueToolkit.load_init_from_template(self.model_simulator, self.model_init_value_filename, \
                                                               ignore_init_mismatch=True)
                return approximated_optimum
            else:
                optimized_input = {}
                for input_key in self.manipulated_variables:
                    var = self.rto_model.manipulated_variables[input_key].__call__(self.model_simulator)
                    optimized_input[input_key] = value(var)
                return optimized_input
        else:
            return approximated_optimum


    def set_new_parameter(self):
        for k in self.parameter_init_values.keys():
            getattr(self.model_simulator, k)._value = value(getattr(self.model_pe_optimizer, k))


    def get_damped_input(self, current_input, optimized_input):
        ret = {}
        for k in current_input.keys():
            ret[k] = current_input[k] + (optimized_input[k]-current_input[k])*(1-self.damping_factor[k])
            if self.mv_bounds[k][0] is not None and ret[k] <= self.mv_bounds[k][0]:
                print("MV %s reaches its lower bound."%k)
                ret[k] = self.mv_bounds[k][0]
            if self.mv_bounds[k][1] is not None and ret[k] >= self.mv_bounds[k][1]:
                print("MV %s reaches its upper bound." % k)
                ret[k] = self.mv_bounds[k][1]
        return ret

    def pe_promoted_solving(self, tee = False, print_level=2):
        for pn in self.parameter_init_values.keys():
            getattr(self.model_pe_optimizer, pn).fixed = True

        for case in self.model_pe_optimizer.CaseIndex:
            block = self.model_pe_optimizer.Cases[case]
            if print_level >= 2:
                print("Preparing initial value for pe.")
            # mv的值来自采样算法
            for k in self.manipulated_variables:
                var = self.rto_model.manipulated_variables[k].__call__(block)
                getattr(block, k + "_prosim").fix(value(var))
                var.fixed = False
            #spec的值来自设定
            for k in self.specifications:
                var = self.rto_model.specifications[k].__call__(block)
                getattr(block, k + "_prosim").fix(self.specification_values[k])
                var.fixed = False
            block.simulation_promoting_obj.activate()

            self.solver.options['obj_scaling_factor'] = self.obj_scaling_factor['Simulation_promoted']
            results = self.solver.solve(block, tee=tee)
            self.check_solver_status(results)

            block.simulation_promoting_obj.deactivate()
            for k in self.manipulated_variables:
                var = self.rto_model.manipulated_variables[k].__call__(block)
                var.fixed = True
            for k in self.specifications:
                var = self.rto_model.specifications[k].__call__(block)
                var.fixed = True

        if not os.path.exists("data/"):
            os.makedirs("data/")
        InitValueToolkit.back_up_var_value(self.model_pe_optimizer,"data/pe_backup%d.txt"%self._id)
        if print_level >= 2:
            print("Solving PE.")
        for pn in self.parameters:
            getattr(self.model_pe_optimizer, pn).fixed = True
        for k in self.active_parameters:
            getattr(self.model_pe_optimizer, k).fixed = False

        results = None
        try:
            self.solver.options['obj_scaling_factor'] = self.obj_scaling_factor['Parameter_estimation']
            results = self.solver.solve(self.model_pe_optimizer, tee=tee)
        except Exception as e:
            print(e)
            print("PE problem restored.")
            InitValueToolkit.load_backup_var_value(self.model_pe_optimizer, "data/pe_backup%d.txt"%self._id)
        if results is not None:
            if not self.check_solver_status(results):
                print("PE problem restored.")
                InitValueToolkit.load_backup_var_value(self.model_pe_optimizer, "data/pe_backup%d.txt"%self._id)
        return results

    def set_center_fixed_parameter(self, xi):
        self.center_fixed_parameter = xi

    def set_new_para_weight(self):
        for k in self.active_parameters:
            std=numpy.abs(value(getattr(self.model_pe_optimizer, k))-self.prev_parameter[k])
            if std < numpy.sqrt(1/self.param_weight[k])/10:
                std = numpy.sqrt(1/self.param_weight[k])/10
            cov_estimated = (3*std)**2
            getattr(self.model_pe_optimizer, k + "_weight").fix(1 / cov_estimated)


    def backup_prev_parameter(self):
        self.prev_parameter={}
        for k in self.active_parameters:
            self.prev_parameter[k] = value(getattr(self.model_pe_optimizer, k))

    def set_specifications(self, rto_iter):
        specifications = self.specification_func.__call__(rto_iter)
        self.last_specification_values = copy.deepcopy(self.specification_values)
        for k,v in specifications.items():
            # 不能直接定死，要写在spec里，否则容易求解失败
            # self.rto_plant.specifications[k].__call__(self.plant_simulator).fix(v)
            # self.rto_model.specifications[k].__call__(self.model_simulator).fix(v)
            # for c in self.model_pe_optimizer.CaseIndex:
            #     self.rto_model.specifications[k].__call__(self.model_pe_optimizer.Cases[c]).fix(v)
            self.specification_values[k] = v


    def do_test_Successive(self, tee = False, update_xi=True, update_para_weight=False, simulation_promotion=True, print_level=2):
        print("Testing GPE.")
        current_point = self.starting_point

        self.specification_values={}
        if self.specification_func is not None:
            self.set_specifications(0)
        for rto_iter in range(self.max_iter):
            self._current_iter=rto_iter

            if print_level >= 1:
                print("RTO Iteration %d:"%rto_iter)
            self.input_history_data[rto_iter] = {}
            for k, v in current_point.items():
                self.input_history_data[rto_iter][k] = v

            plant_data = self.get_noised_plant_data(rto_iter, current_point, tee=tee, simulation_promotion=simulation_promotion, print_level=print_level)
            self.plant_history_data[rto_iter] = {}
            for k,v in plant_data[0].items():
                self.plant_history_data[rto_iter][k] = v

            model_output_now = self.get_model_simulation(current_point, tee=tee, simulation_promotion=simulation_promotion, print_level=print_level)
            self.model_history_data[rto_iter] = {}

            for k, v in model_output_now.items():
                self.model_history_data[rto_iter][k] = v
            for k in self.parameter_init_values.keys():
                self.model_history_data[rto_iter][k] = value(getattr(self.model_pe_optimizer,k))

            self.set_measure_for_MultiPE(plant_data)

            if update_para_weight:
                self.backup_prev_parameter()

            basepoint = current_point
            self.update_model_basepoint(basepoint)

            if not simulation_promotion:
                if print_level >= 1:
                    print('Solving PE')
                self.solver.options['obj_scaling_factor'] = self.obj_scaling_factor['Parameter_estimation']
                results = self.solver.solve(self.model_pe_optimizer, tee=tee)
            else:
                results = self.pe_promoted_solving(tee=tee, print_level=print_level)

            self.set_new_parameter()

            if update_xi:
                if hasattr(self, 'center_fixed_parameter'):
                    for k in self.active_parameters:
                        if k not in self.center_fixed_parameter:
                            getattr(self.model_pe_optimizer, k + "_expectation").fix(value(getattr(self.model_pe_optimizer, k)))
                else:
                    for k in self.active_parameters:
                        getattr(self.model_pe_optimizer, k + "_expectation").fix(value(getattr(self.model_pe_optimizer, k)))

            if update_para_weight:
                self.set_new_para_weight()
                for k in self.active_parameters:
                    self.model_history_data[rto_iter][k + '_weight'] = value(getattr(self.model_pe_optimizer, k + "_weight"))

            if self.specification_func is not None:
                self.set_specifications(rto_iter)
            optimized_input = self.optimize_for_u(tee=tee, print_level=print_level)

            # if the working condition is changed, then the optimum need not to be damped
            if is_value_dict_same(self.last_specification_values, self.specification_values):
                current_point = self.get_damped_input(current_point, optimized_input)
            else:
                current_point = optimized_input

            self.post_iteration_callback.__call__(self, rto_iter)
