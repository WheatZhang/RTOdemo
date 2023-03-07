#!/usr/bin/env python
#-*- coding:utf-8 -*-
import numpy

from rtolib.core.algo.common import MA_type_Algorithm
from rtolib.core.pyomo_model import *
from rtolib.core.solve import PyomoSimulator,PyomoOptimizer,TrustRegionOptimizer, PenaltyTrustRegionOptimizer,\
CompoStepTrustRegionOptimizer
import copy


class ModifierAdaptation(MA_type_Algorithm):

    def set_problem(self, problem_description,
                    plant,
                    model,
                    perturbation_method,
                    noise_generator,
                    solver_executable,
                    spec_function,
                    modifier_type,
                    skipped_modifiers=[],
                    ):
        self.problem_description = problem_description
        self.plant_simulator = PyomoSimulator(plant)
        mvs = self.problem_description.symbol_list['MV']
        if modifier_type == ModifierType.RTO:
            cvs = [self.problem_description.symbol_list['OBJ']]
            for cv in self.problem_description.symbol_list['CON']:
                cvs.append(cv)
        elif modifier_type == ModifierType.OUTPUT:
            cvs = self.problem_description.symbol_list['CV']
        else:
            raise ValueError("Not applicable")
        self.model_simulator = PyomoSimulator(PyomoModelWithModifiers(model, modifier_type, mvs, cvs))
        self.model_optimizer = PyomoOptimizer(PyomoModelWithModifiers(model, modifier_type, mvs, cvs))
        self.perturbation_method = perturbation_method
        self.noise_generator = noise_generator
        self.solver_executable = solver_executable
        self.spec_function = spec_function
        self.skipped_modifiers = skipped_modifiers


    def set_homotopy_var(self, homotopy_var):
        '''

        :param homotopy_var: list of names
        :return:
        '''
        self.plant_simulator.homotopy_var = homotopy_var
        self.model_simulator.homotopy_var = homotopy_var
        self.model_optimizer.homotopy_var = homotopy_var


    def initialize_simulation(self, starting_point, initial_parameter_value):
        # build model
        self.plant_simulator.build(self.problem_description)
        self.model_simulator.build(self.problem_description)
        self.model_optimizer.build(self.problem_description)

        # set solver
        # TODO: solver parameter tuning
        # TODO: each solver should be seperate
        default_options={'max_iter':500,
                         "tol":1e-10}
        solver1 = SolverFactory('ipopt', executable=self.solver_executable)
        self.plant_simulator.set_solver(solver1, tee=False, default_options=default_options)
        solver2 = SolverFactory('ipopt', executable=self.solver_executable)
        self.model_simulator.set_solver(solver2, tee=False, default_options=default_options)
        solver3 = SolverFactory('ipopt', executable=self.solver_executable)
        self.model_optimizer.set_solver(solver3, tee=False, default_options=default_options)

        self.iter_count = 0
        self.input_history_data[0] = {}
        self.set_current_point(starting_point)

        self.iter_count = 1

        self.model_history_data[0] = {}
        self.set_initial_model_parameters(initial_parameter_value)

        # set basepoint
        self.model_simulator.set_base_point(self.current_point)
        self.model_optimizer.set_base_point(self.current_point)


    def one_step_simulation(self):
        print("Iteration %d"%self.iter_count)

        # initialize storage
        self.model_history_data[self.iter_count] = {}
        self.plant_history_data[self.iter_count] = {}
        self.input_history_data[self.iter_count] = {}

        # get trial point
        trial_points = self.get_trial_point()

        # get plant simulation result
        plant_output_data = self.get_plant_simulation_result(trial_points)

        # get model simulation result with and without modifiers
        model_output_data = self.get_model_simulation_result(trial_points)

        # update modifiers
        self.update_modifiers(plant_output_data, model_output_data)

        # set basepoint
        self.model_simulator.set_base_point(self.current_point)
        self.model_optimizer.set_base_point(self.current_point)

        self.store_model_adaptation_data()

        optimized_input, solve_status = self.optimize_for_u()

        mv_bounds = self.problem_description.bounds
        filtered_input = self.filter_mv(optimized_input, mv_bounds)

        # set specification and store input data
        self.set_current_point(filtered_input)

        # iter count
        self.iter_count += 1


class ModifierAdaptationTR(ModifierAdaptation):

    def set_problem(self, problem_description,
                    plant,
                    model,
                    perturbation_method,
                    noise_generator,
                    solver_executable,
                    spec_function,
                    modifier_type,
                    skipped_modifiers=[],
                    ):
        self.problem_description = problem_description
        self.plant_simulator = PyomoSimulator(plant)
        mvs = self.problem_description.symbol_list['MV']
        if modifier_type == ModifierType.RTO:
            cvs = [self.problem_description.symbol_list['OBJ']]
            for cv in self.problem_description.symbol_list['CON']:
                cvs.append(cv)
        elif modifier_type == ModifierType.OUTPUT:
            cvs = self.problem_description.symbol_list['CV']
        else:
            raise ValueError("Not applicable")
        self.model_simulator = PyomoSimulator(PyomoModelWithModifiers(model, modifier_type, mvs, cvs))
        self.model_optimizer = TrustRegionOptimizer(PyomoModelWithModifiers(model, modifier_type, mvs, cvs))
        self.perturbation_method = perturbation_method
        self.noise_generator = noise_generator
        self.solver_executable = solver_executable
        self.spec_function = spec_function
        self.skipped_modifiers = skipped_modifiers

    def register_option(self):
        super().register_option()
        self.available_options["eta1"] = ("[0,1]", 0.1)
        self.available_options["eta2"] = ("[0,1]", 0.95)
        self.available_options["gamma1"] = ("[0,1]", 0.5)
        self.available_options["gamma2"] = ("[0,1]", 1)
        self.available_options["gamma3"] = ("[1, inf]", 1.5)
        self.available_options["max_iter"] = ("[1, inf]", 100)
        self.available_options["feasibility_tol"] = ("positive float", 1e-8)
        self.available_options["stationarity_tol"] = ("positive float", 1e-6)
        self.available_options["max_trust_radius"] = ("positive float", 10)
        self.available_options["initial_trust_radius"] = ("positive float", 1)


    def optimize_for_u(self, tr_radius, tr_base):
        input_values = {}
        if self.spec_function is not None:
            for k, v in self.current_spec.items():
                input_values[k] = v
        optimized_input, solve_status = self.model_optimizer.optimize(input_values, tr_radius, tr_base,
                                                                      param_values=self.current_parameter_value,
                                                                      use_homo=self.options["homotopy_optimization"])
        # TODO:deal with solve status, fallback strategy
        return optimized_input, solve_status

    def initialize_simulation(self, starting_point, initial_parameter_value):
        super().initialize_simulation(starting_point, initial_parameter_value)
        self.trust_radius = self.options["initial_trust_radius"]


    def one_step_simulation(self):
        print("Iteration %d"%self.iter_count)

        # initialize storage
        self.model_history_data[self.iter_count] = {}
        self.plant_history_data[self.iter_count] = {}
        self.input_history_data[self.iter_count] = {}
        if self.iter_count == 1:
            # If we want to add some print information, do it here.
            self.model_history_data[1]['rho'] = ""
            self.model_history_data[1]['tr'] = ""
            self.model_history_data[1]['tr_adjusted'] = ""
            self.model_history_data[1]['event'] = ""

        # get trial point
        trial_points = self.get_trial_point()
        base_iteration = self.iter_count
        base_input = copy.deepcopy(trial_points[0])

        # get plant simulation result
        plant_output_data = self.get_plant_simulation_result(trial_points)

        # get model simulation result with and without modifiers
        model_output_data = self.get_model_simulation_result(trial_points)

        # update modifiers
        self.update_modifiers(plant_output_data, model_output_data)

        # set basepoint
        self.model_simulator.set_base_point(self.current_point)
        self.model_optimizer.set_base_point(self.current_point)

        self.store_model_adaptation_data()

        # culculate_base_obj_value
        base_model_output, solve_status = self.model_simulator.simulate(base_input,
                                                                        param_values=self.current_parameter_value,
                                                                        use_homo=self.options["homotopy_simulation"])
        obj_var_name = self.problem_description.symbol_list['OBJ']
        base_obj_value = base_model_output[obj_var_name]

        iter_successful_flag=False
        rho = 100
        tr_base={}
        for k,v in trial_points[0].items():
            if k in self.problem_description.symbol_list['MV']:
                tr_base[k]=v
        while(not iter_successful_flag):
            self.model_history_data[self.iter_count]['tr'] = self.trust_radius
            optimized_input, solve_status = self.optimize_for_u(self.trust_radius, tr_base)

            if solve_status == PyomoModelSolvingStatus.OPTIMIZATION_FAILED:
                # starting point infeasible
                self.set_current_point(self.current_point)
                self.iter_count += 1
                self.model_history_data[self.iter_count] = {}
                self.plant_history_data[self.iter_count] = {}
                self.input_history_data[self.iter_count] = {}
                plant_trial_point_output = self.get_plant_simulation_result([self.current_point])[0]
                model_trial_point_output = self.get_model_simulation_result([self.current_point])[0]
                self.model_history_data[self.iter_count - 1]['rho'] = -1
                self.model_history_data[self.iter_count - 1]['event'] = "subproblem infeasible"
                break

            mv_bounds = self.problem_description.bounds
            filtered_input = self.adapt_to_bound_mv(optimized_input, mv_bounds)

            # set specification and store input data
            self.set_current_point(filtered_input)

            # iter count
            if self.iter_count >= self.options['max_iter']:
                self.iter_count += 1
                break
            self.iter_count += 1
            self.model_history_data[self.iter_count] = {}
            self.plant_history_data[self.iter_count] = {}
            self.input_history_data[self.iter_count] = {}

            plant_trial_point_output = self.get_plant_simulation_result([filtered_input])[0]
            model_trial_point_output = self.get_model_simulation_result([filtered_input])[0]

            # accept the trial point
            flag_infeasible = False
            for con_name in self.problem_description.symbol_list['CON']:
                if plant_trial_point_output[con_name] > self.options['feasibility_tol']:
                    flag_infeasible = True
                    rho = -1
                    break
            obj_name = self.problem_description.symbol_list['OBJ']
            if base_obj_value < self.model_history_data[self.iter_count][obj_name]:
                # Because the solver is not a global optimization solver, it is possible
                # that the model obj function increases. In this case, we ask for backtracking.
                self.model_history_data[self.iter_count - 1]['event'] = "model merit increases"
                rho = -2
            else:
                if not flag_infeasible:
                    if abs(plant_output_data[0][obj_name]-plant_trial_point_output[obj_name])<self.options['stationarity_tol']:
                        iter_successful_flag=True
                        self.model_history_data[self.iter_count - 1]['event'] = "plant converges"
                        self.model_history_data[self.iter_count-1]['rho'] = 1
                        continue
                    else:
                        if base_obj_value - self.model_history_data[self.iter_count][obj_name] > 1e-6:
                            rho=(self.plant_history_data[base_iteration][obj_name] - self.plant_history_data[self.iter_count][obj_name])/ \
                                (base_obj_value - self.model_history_data[self.iter_count][obj_name])
                        else:
                            rho = (self.plant_history_data[base_iteration][obj_name] -
                                   self.plant_history_data[self.iter_count][obj_name]) / 1e-6

            self.model_history_data[self.iter_count-1]['rho'] = rho

            # update trust-region radius
            if rho < self.options["eta1"]:
                gamma = self.options['gamma1']
            elif rho < self.options["eta2"]:
                gamma = self.options['gamma2']
                iter_successful_flag = True
            else:
                gamma = self.options['gamma3']
                iter_successful_flag = True
            self.trust_radius *= gamma
            if self.trust_radius < 1e-8:
                self.trust_radius = 1e-8
                self.model_history_data[self.iter_count - 1]['event'] = "trust region too small"
            if self.trust_radius > self.options['max_trust_radius']:
                self.trust_radius = self.options['max_trust_radius']
                self.model_history_data[self.iter_count - 1]['event'] = "trust region too large"
            self.model_history_data[self.iter_count - 1]['tr_adjusted'] = self.trust_radius

class ModifierAdaptationMaxTR(ModifierAdaptation):

    def set_problem(self, problem_description,
                    plant,
                    model,
                    perturbation_method,
                    noise_generator,
                    solver_executable,
                    spec_function,
                    modifier_type,
                    skipped_modifiers=[],
                    ):
        self.problem_description = problem_description
        self.plant_simulator = PyomoSimulator(plant)
        mvs = self.problem_description.symbol_list['MV']
        if modifier_type == ModifierType.RTO:
            cvs = [self.problem_description.symbol_list['OBJ']]
            for cv in self.problem_description.symbol_list['CON']:
                cvs.append(cv)
        elif modifier_type == ModifierType.OUTPUT:
            cvs = self.problem_description.symbol_list['CV']
        else:
            raise ValueError("Not applicable")
        self.model_simulator = PyomoSimulator(PyomoModelWithModifiers(model, modifier_type, mvs, cvs))
        self.model_optimizer = TrustRegionOptimizer(PyomoModelWithModifiers(model, modifier_type, mvs, cvs))
        self.perturbation_method = perturbation_method
        self.noise_generator = noise_generator
        self.solver_executable = solver_executable
        self.spec_function = spec_function
        self.skipped_modifiers = skipped_modifiers

    def register_option(self):
        super().register_option()
        self.available_options["max_iter"] = ("[1, inf]", 100)
        self.available_options["max_trust_radius"] = ("positive float", 10)

    def optimize_for_u(self):
        input_values = {}
        if self.spec_function is not None:
            for k, v in self.current_spec.items():
                input_values[k] = v
        tr_base = {}
        for k,v in self.current_point.items():
            if k in self.problem_description.symbol_list['MV']:
                tr_base[k]=v
        optimized_input, solve_status = self.model_optimizer.optimize(input_values, self.options['max_trust_radius'], tr_base,
                                                                      param_values=self.current_parameter_value,
                                                                      use_homo=self.options["homotopy_optimization"])
        # TODO:deal with solve status, fallback strategy
        return optimized_input, solve_status



class ModifierAdaptationPenaltyTR(ModifierAdaptationTR):

    def set_problem(self, problem_description,
                    plant,
                    model,
                    perturbation_method,
                    noise_generator,
                    solver_executable,
                    spec_function,
                    modifier_type,
                    skipped_modifiers=[],
                    ):
        self.problem_description = problem_description
        self.plant_simulator = PyomoSimulator(plant)
        mvs = self.problem_description.symbol_list['MV']
        if modifier_type == ModifierType.RTO:
            cvs = [self.problem_description.symbol_list['OBJ']]
            for cv in self.problem_description.symbol_list['CON']:
                cvs.append(cv)
        elif modifier_type == ModifierType.OUTPUT:
            cvs = self.problem_description.symbol_list['CV']
        else:
            raise ValueError("Not applicable")
        self.model_simulator = PyomoSimulator(PyomoModelWithModifiers(model, modifier_type, mvs, cvs))
        self.model_optimizer = PenaltyTrustRegionOptimizer(PyomoModelWithModifiers(model, modifier_type, mvs, cvs))
        self.perturbation_method = perturbation_method
        self.noise_generator = noise_generator
        self.solver_executable = solver_executable
        self.spec_function = spec_function
        self.skipped_modifiers = skipped_modifiers

    def register_option(self):
        super().register_option()
        self.available_options["sigma"] = ("positive float", 10)

    def initialize_simulation(self, starting_point, initial_parameter_value):
        super().initialize_simulation(starting_point, initial_parameter_value)
        self.model_optimizer.set_penalty_coeff(self.options['sigma'])

    def one_step_simulation(self):
        print("Iteration %d" % self.iter_count)

        # initialize storage
        self.model_history_data[self.iter_count] = {}
        self.plant_history_data[self.iter_count] = {}
        self.input_history_data[self.iter_count] = {}
        if self.iter_count == 1:
            # If we want to add some print information, do it here.
            self.model_history_data[1]['rho'] = ""
            self.model_history_data[1]['tr'] = ""
            self.model_history_data[1]['tr_adjusted'] = ""
            self.model_history_data[1]['event'] = ""
            obj_var_name = self.problem_description.symbol_list['OBJ']
            self.model_history_data[1]['base_'+obj_var_name] = ""
            for con_var_name in self.problem_description.symbol_list['CON']:
                self.model_history_data[1]['base_'+con_var_name] = ""
            self.model_history_data[1]['merit'] = ""
            self.model_history_data[1]['base_merit'] = ""
            self.plant_history_data[1]['merit'] = ""
            self.plant_history_data[1]['base_merit'] = ""

        # get trial point
        trial_points = self.get_trial_point()
        base_iteration = self.iter_count
        base_input = copy.deepcopy(trial_points[0])

        # get plant simulation result
        plant_output_data = self.get_plant_simulation_result(trial_points)
        obj_var_name = self.problem_description.symbol_list['OBJ']
        ret = self.plant_history_data[self.iter_count][obj_var_name] / \
              self.problem_description.scaling_factors[obj_var_name]
        infeasibility_sq = 0
        for con_var_name in self.problem_description.symbol_list['CON']:
            infeasibility_sq += max(self.plant_history_data[self.iter_count][con_var_name],0)**2 / \
                   self.problem_description.scaling_factors[con_var_name]**2
        ret += numpy.sqrt(infeasibility_sq)* self.options['sigma']
        self.plant_history_data[self.iter_count]['merit'] = ret
        self.plant_history_data[self.iter_count]['base_merit'] = ret

        # get model simulation result with and without modifiers
        model_output_data = self.get_model_simulation_result(trial_points)
        obj_var_name = self.problem_description.symbol_list['OBJ']
        ret = self.model_history_data[self.iter_count][obj_var_name] / \
              self.problem_description.scaling_factors[obj_var_name]
        infeasibility_sq = 0
        for con_var_name in self.problem_description.symbol_list['CON']:
            infeasibility_sq += max(self.model_history_data[self.iter_count][con_var_name], 0) ** 2 / \
                                self.problem_description.scaling_factors[con_var_name] ** 2
        ret += numpy.sqrt(infeasibility_sq) * self.options['sigma']
        self.model_history_data[self.iter_count]['merit'] = ret

        # update modifiers
        self.update_modifiers(plant_output_data, model_output_data)

        # set basepoint
        self.model_simulator.set_base_point(self.current_point)
        self.model_optimizer.set_base_point(self.current_point)

        self.store_model_adaptation_data()

        # culculate_base_merit
        base_model_output, solve_status = self.model_simulator.simulate(base_input, param_values=self.current_parameter_value,
                                                                  use_homo=self.options["homotopy_simulation"])
        obj_var_name = self.problem_description.symbol_list['OBJ']
        ret = base_model_output[obj_var_name] / \
              self.problem_description.scaling_factors[obj_var_name]
        infeasibility_sq = 0
        for con_var_name in self.problem_description.symbol_list['CON']:
            infeasibility_sq += max(base_model_output[con_var_name], 0) ** 2 / \
                                self.problem_description.scaling_factors[con_var_name] ** 2
        ret += numpy.sqrt(infeasibility_sq) * self.options['sigma']
        base_merit = ret

        self.model_history_data[self.iter_count]['base_'+obj_var_name] = base_model_output[obj_var_name]
        for con_var_name in self.problem_description.symbol_list['CON']:
            self.model_history_data[self.iter_count]['base_' + con_var_name] = base_model_output[con_var_name]
        self.model_history_data[self.iter_count]['base_merit'] = base_merit

        iter_successful_flag = False
        rho = 100
        tr_base = {}
        for k, v in trial_points[0].items():
            if k in self.problem_description.symbol_list['MV']:
                tr_base[k] = v
        while (not iter_successful_flag):
            self.model_history_data[self.iter_count]['tr'] = self.trust_radius
            try:
                optimized_input, solve_status = self.optimize_for_u(self.trust_radius, tr_base)
            except Exception as e:
                optimized_input = self.current_point
                self.model_history_data[self.iter_count]['event'] = "optimization failed"
            if solve_status == PyomoModelSolvingStatus.OPTIMIZATION_FAILED:
                optimized_input = self.current_point
                self.model_history_data[self.iter_count]['event'] = "optimization failed"

            mv_bounds = self.problem_description.bounds
            filtered_input = self.adapt_to_bound_mv(optimized_input, mv_bounds)

            # set specification and store input data
            self.set_current_point(filtered_input)

            # iter count
            if self.iter_count >= self.options['max_iter']:
                self.iter_count += 1
                break

            self.iter_count += 1
            self.model_history_data[self.iter_count] = {}
            self.plant_history_data[self.iter_count] = {}
            self.input_history_data[self.iter_count] = {}

            plant_trial_point_output = self.get_plant_simulation_result([filtered_input])[0]
            obj_var_name = self.problem_description.symbol_list['OBJ']
            ret = self.plant_history_data[self.iter_count][obj_var_name] / \
                  self.problem_description.scaling_factors[obj_var_name]
            infeasibility_sq = 0
            for con_var_name in self.problem_description.symbol_list['CON']:
                infeasibility_sq += max(self.plant_history_data[self.iter_count][con_var_name], 0) ** 2 / \
                                    self.problem_description.scaling_factors[con_var_name] ** 2
            ret += numpy.sqrt(infeasibility_sq) * self.options['sigma']
            self.plant_history_data[self.iter_count]['merit'] = ret
            model_trial_point_output = self.get_model_simulation_result([filtered_input])[0]
            obj_var_name = self.problem_description.symbol_list['OBJ']
            ret = self.model_history_data[self.iter_count][obj_var_name] / \
                  self.problem_description.scaling_factors[obj_var_name]
            infeasibility_sq = 0
            for con_var_name in self.problem_description.symbol_list['CON']:
                infeasibility_sq += max(self.model_history_data[self.iter_count][con_var_name], 0) ** 2 / \
                                    self.problem_description.scaling_factors[con_var_name] ** 2
            ret += numpy.sqrt(infeasibility_sq) * self.options['sigma']
            self.model_history_data[self.iter_count]['merit'] = ret

            # accept the trial point
            flag_infeasible = False
            for con_name in self.problem_description.symbol_list['CON']:
                if plant_trial_point_output[con_name] > self.options['feasibility_tol']:
                    flag_infeasible = True
                    break
            if base_merit < self.model_history_data[self.iter_count]['merit']:
                # Because the solver is not a global optimization solver, it is possible
                # that the model merit function increases. In this case, we ask for backtracking.
                self.model_history_data[self.iter_count-1]['event'] = "model merit increases"
                rho = -2
                # if self.plant_history_data[base_iteration]['merit'] -\
                #                self.plant_history_data[self.iter_count]['merit'] >0:
                #     rho = 0.1
                # else:
                #     rho = -2
            else:
                if (not flag_infeasible) and\
                            abs(self.plant_history_data[base_iteration]['merit'] - self.plant_history_data[self.iter_count]['merit']) < self.options[
                    'stationarity_tol']:
                    iter_successful_flag = True
                    self.model_history_data[self.iter_count-1]['event'] = "plant converges"
                    self.model_history_data[self.iter_count-1]['rho'] = 1
                    continue
                else:
                    if base_merit - self.model_history_data[self.iter_count]['merit'] >1e-6:
                        rho = (self.plant_history_data[base_iteration]['merit'] - self.plant_history_data[self.iter_count]['merit']) / \
                              (base_merit - self.model_history_data[self.iter_count]['merit'])
                    else:
                        rho = (self.plant_history_data[base_iteration]['merit'] -
                               self.plant_history_data[self.iter_count]['merit']) / 1e-6
            self.model_history_data[self.iter_count-1]['rho'] = rho

            # update trust-region radius
            if rho < self.options["eta1"]:
                gamma = self.options['gamma1']
            elif rho < self.options["eta2"]:
                gamma = self.options['gamma2']
                iter_successful_flag = True
            else:
                gamma = self.options['gamma3']
                iter_successful_flag = True
            self.trust_radius *= gamma
            if self.trust_radius < 1e-8:
                self.trust_radius = 1e-8
                self.model_history_data[self.iter_count - 1]['event'] = "trust region too small"
            if self.trust_radius > self.options['max_trust_radius']:
                self.trust_radius = self.options['max_trust_radius']
                self.model_history_data[self.iter_count - 1]['event'] = "trust region too large"
            self.model_history_data[self.iter_count-1]['tr_adjusted'] = self.trust_radius


class ModifierAdaptationCompoStepTR(ModifierAdaptationPenaltyTR):

    def set_problem(self, problem_description,
                    plant,
                    model,
                    perturbation_method,
                    noise_generator,
                    solver_executable,
                    spec_function,
                    modifier_type,
                    skipped_modifiers=[],
                    ):
        self.problem_description = problem_description
        self.plant_simulator = PyomoSimulator(plant)
        mvs = self.problem_description.symbol_list['MV']
        if modifier_type == ModifierType.RTO:
            cvs = [self.problem_description.symbol_list['OBJ']]
            for cv in self.problem_description.symbol_list['CON']:
                cvs.append(cv)
        elif modifier_type == ModifierType.OUTPUT:
            cvs = self.problem_description.symbol_list['CV']
        else:
            raise ValueError("Not applicable")
        self.model_simulator = PyomoSimulator(PyomoModelWithModifiers(model, modifier_type, mvs, cvs))
        self.model_optimizer = CompoStepTrustRegionOptimizer(PyomoModelWithModifiers(model, modifier_type, mvs, cvs))
        self.perturbation_method = perturbation_method
        self.noise_generator = noise_generator
        self.solver_executable = solver_executable
        self.spec_function = spec_function
        self.skipped_modifiers = skipped_modifiers

    def register_option(self):
        super().register_option()
        self.available_options["xi_N"] = ([0,1], 0.5)

    def initialize_simulation(self, starting_point, initial_parameter_value):
        ModifierAdaptationTR.initialize_simulation(self, starting_point, initial_parameter_value)

    def optimize_for_u(self, tr_radius, tr_base):
        input_values = {}
        if self.spec_function is not None:
            for k, v in self.current_spec.items():
                input_values[k] = v
        optimized_input, solve_status = self.model_optimizer.optimize(input_values, tr_radius, tr_base, self.options['xi_N'],
                                                                      param_values=self.current_parameter_value,
                                                                      use_homo=self.options["homotopy_optimization"])
        return optimized_input, solve_status

