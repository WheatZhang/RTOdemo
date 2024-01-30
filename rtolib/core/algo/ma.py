#!/usr/bin/env python
#-*- coding:utf-8 -*-
import numpy
from rtolib.core.algo.common import MA_type_Algorithm
from rtolib.core.pyomo_model import *
from rtolib.core.solve import PyomoSimulator,PyomoOptimizer,TrustRegionOptimizer, PenaltyTrustRegionOptimizer,\
CompoStepTrustRegionOptimizer, CompoStepTrustRegionBBMOptimizer, BlackBoxOptimizer
from rtolib.core.black_box_model import BlackBoxModelWithModifiers, BlackBoxModel
import copy
from rtolib.util.init_value import to_template, load_init_from_template


# TODO: solver parameter tuning
# TODO: each solver should be seperate
# TODO: fallback policy
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
        default_options={'max_iter':500,
                         "tol":1e-10}
        solver1 = SolverFactory('ipopt', executable=self.solver_executable)
        self.plant_simulator.set_solver(solver1, tee=False, default_options=default_options)
        solver2 = SolverFactory('ipopt', executable=self.solver_executable)
        self.model_simulator.set_solver(solver2, tee=False, default_options=default_options)
        if isinstance(self.model_optimizer, PyomoOptimizer):
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

        # initialize model optimizer
        to_template(self.model_simulator.model, "temp_init_file.txt")
        load_init_from_template(self.model_optimizer.model, "temp_init_file.txt", \
                                ignore_init_mismatch=True)

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

            if solve_status == ModelSolvingStatus.OPTIMIZATION_FAILED:
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
        self.sigma = self.options['sigma']
        self.model_optimizer.set_penalty_coeff(self.sigma)

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
        ret += numpy.sqrt(infeasibility_sq)* self.sigma
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
        ret += numpy.sqrt(infeasibility_sq) * self.sigma
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
        ret += numpy.sqrt(infeasibility_sq) * self.sigma
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
            if solve_status == ModelSolvingStatus.OPTIMIZATION_FAILED:
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
            ret += numpy.sqrt(infeasibility_sq) * self.sigma
            self.plant_history_data[self.iter_count]['merit'] = ret
            model_trial_point_output = self.get_model_simulation_result([filtered_input])[0]
            obj_var_name = self.problem_description.symbol_list['OBJ']
            ret = self.model_history_data[self.iter_count][obj_var_name] / \
                  self.problem_description.scaling_factors[obj_var_name]
            infeasibility_sq = 0
            for con_var_name in self.problem_description.symbol_list['CON']:
                infeasibility_sq += max(self.model_history_data[self.iter_count][con_var_name], 0) ** 2 / \
                                    self.problem_description.scaling_factors[con_var_name] ** 2
            ret += numpy.sqrt(infeasibility_sq) * self.sigma
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
                    black_box_model=None,
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
        if black_box_model is None:
            self.model_optimizer = CompoStepTrustRegionOptimizer(PyomoModelWithModifiers(model, modifier_type, mvs, cvs))
        else:
            self.model_optimizer = CompoStepTrustRegionBBMOptimizer(BlackBoxModelWithModifiers(black_box_model, mvs, cvs))
        self.perturbation_method = perturbation_method
        self.noise_generator = noise_generator
        self.solver_executable = solver_executable
        self.spec_function = spec_function
        self.skipped_modifiers = skipped_modifiers

    def register_option(self):
        super().register_option()
        self.available_options["xi_N"] = ("[0,1]", 0.5)
        self.available_options["sigma_inc"] = ("[1,Inf]", 2)
        self.available_options["sigma_kappa"] = ("(0,1)", 0.25)
        self.available_options['adaptive_sigma'] = ("bool", False)

    def initialize_simulation(self, starting_point, initial_parameter_value):
        ModifierAdaptationTR.initialize_simulation(self, starting_point, initial_parameter_value)
        self.sigma = self.options['sigma']

    def optimize_for_u(self, tr_radius, tr_base):
        input_values = {}
        if self.spec_function is not None:
            for k, v in self.current_spec.items():
                input_values[k] = v
        if isinstance(self.model_optimizer, PyomoOptimizer):
            optimized_input, input_after_normal_step, solve_status = \
                self.model_optimizer.optimize(input_values, tr_radius, tr_base, self.options['xi_N'],
                                                                      param_values=self.current_parameter_value,
                                                                      use_homo=self.options["homotopy_optimization"])
        elif isinstance(self.model_optimizer, BlackBoxOptimizer):
            optimized_input, input_after_normal_step, solve_status = \
                self.model_optimizer.optimize(tr_radius, tr_base, self.options['xi_N'], \
                                              modifiers_value=self.current_parameter_value)
        else:
            raise TypeError("unknown optimizer type")
        return optimized_input, input_after_normal_step, solve_status

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
            self.model_history_data[1]['sigma'] = ""
            self.model_history_data[1]['merit'] = ""
            self.model_history_data[1]['infeasibility'] = ""
            self.model_history_data[1]['obj_for_merit'] = ""
            self.model_history_data[1]['base_merit'] = ""
            self.model_history_data[1]['base_infeasibility'] = ""
            self.model_history_data[1]['base_obj_for_merit'] = ""

            self.plant_history_data[1]['merit'] = ""
            self.plant_history_data[1]['infeasibility'] = ""
            self.plant_history_data[1]['obj_for_merit'] = ""
            self.plant_history_data[1]['base_merit'] = ""
            self.plant_history_data[1]['base_infeasibility'] = ""
            self.plant_history_data[1]['base_obj_for_merit'] = ""

        # get trial point
        trial_points = self.get_trial_point()
        base_iteration = self.iter_count
        base_input = copy.deepcopy(trial_points[0])

        # get plant simulation result
        plant_output_data = self.get_plant_simulation_result(trial_points)
        obj_var_name = self.problem_description.symbol_list['OBJ']
        obj = self.plant_history_data[self.iter_count][obj_var_name] / \
              self.problem_description.scaling_factors[obj_var_name]
        self.plant_history_data[self.iter_count]['obj_for_merit'] = obj
        self.plant_history_data[base_iteration]['base_obj_for_merit'] = obj
        infeasibility_sq = 0
        for con_var_name in self.problem_description.symbol_list['CON']:
            infeasibility_sq += max(self.plant_history_data[self.iter_count][con_var_name],0)**2 / \
                   self.problem_description.scaling_factors[con_var_name]**2
        infeasibility = numpy.sqrt(infeasibility_sq)
        self.plant_history_data[self.iter_count]['infeasibility'] = infeasibility
        self.plant_history_data[base_iteration]['base_infeasibility'] = infeasibility
        merit = obj+infeasibility* self.sigma
        self.plant_history_data[self.iter_count]['merit'] = merit
        self.plant_history_data[base_iteration]['base_merit'] = merit

        # get model simulation result with and without modifiers
        model_output_data = self.get_model_simulation_result(trial_points)
        obj_var_name = self.problem_description.symbol_list['OBJ']
        obj = self.model_history_data[self.iter_count][obj_var_name] / \
              self.problem_description.scaling_factors[obj_var_name]
        self.model_history_data[self.iter_count]['obj_for_merit'] = obj
        infeasibility_sq = 0
        for con_var_name in self.problem_description.symbol_list['CON']:
            infeasibility_sq += max(self.model_history_data[self.iter_count][con_var_name], 0) ** 2 / \
                                self.problem_description.scaling_factors[con_var_name] ** 2
        infeasibility = numpy.sqrt(infeasibility_sq)
        self.model_history_data[self.iter_count]['infeasibility'] = infeasibility
        merit = obj + infeasibility * self.sigma
        self.model_history_data[self.iter_count]['merit'] = merit

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
        base_obj = base_model_output[obj_var_name] / \
              self.problem_description.scaling_factors[obj_var_name]
        infeasibility_sq = 0
        for con_var_name in self.problem_description.symbol_list['CON']:
            infeasibility_sq += max(base_model_output[con_var_name], 0) ** 2 / \
                                self.problem_description.scaling_factors[con_var_name] ** 2
        base_infeasibility = numpy.sqrt(infeasibility_sq)

        self.model_history_data[base_iteration]['base_'+obj_var_name] = base_model_output[obj_var_name]
        self.model_history_data[base_iteration]['base_obj_for_merit'] =base_obj
        for con_var_name in self.problem_description.symbol_list['CON']:
            self.model_history_data[base_iteration]['base_' + con_var_name] = base_model_output[con_var_name]
        self.model_history_data[base_iteration]['base_infeasibility'] = base_infeasibility
        merit = base_obj + base_infeasibility * self.sigma
        self.model_history_data[base_iteration]['base_merit'] = merit

        iter_successful_flag = False
        rho = 100
        tr_base = {}
        for k, v in trial_points[0].items():
            if k in self.problem_description.symbol_list['MV']:
                tr_base[k] = v
        while (not iter_successful_flag):
            self.model_history_data[self.iter_count]['tr'] = self.trust_radius
            try:
                optimized_input, input_after_normal_step, solve_status = self.optimize_for_u(self.trust_radius, tr_base)
            except Exception as e:
                # raise e
                optimized_input = self.current_point
                self.model_history_data[self.iter_count]['event'] = "optimization failed"
            if solve_status == ModelSolvingStatus.OPTIMIZATION_FAILED:
                optimized_input = self.current_point
                self.model_history_data[self.iter_count]['event'] = "optimization failed"

            self.model_history_data[self.iter_count]['sigma'] = self.sigma
            # bound the input
            mv_bounds = self.problem_description.bounds
            print("optimized_input calculated by the tr problem")
            print(optimized_input)
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

            if self.options['adaptive_sigma']:
                # check sigma is large enough
                model_trial_point_output = self.get_model_simulation_result([filtered_input])[0]
                obj_var_name = self.problem_description.symbol_list['OBJ']
                obj_after_optimization = self.model_history_data[self.iter_count][obj_var_name] / \
                      self.problem_description.scaling_factors[obj_var_name]
                infeasibility_sq = 0
                for con_var_name in self.problem_description.symbol_list['CON']:
                    infeasibility_sq += max(self.model_history_data[self.iter_count][con_var_name], 0) ** 2 / \
                                        self.problem_description.scaling_factors[con_var_name] ** 2
                infeasibility = numpy.sqrt(infeasibility_sq)
                merit = obj_after_optimization + infeasibility * self.sigma
                base_infeasibility = self.model_history_data[base_iteration]['base_infeasibility']
                base_obj = self.model_history_data[base_iteration]['base_obj_for_merit']
                base_merit =  base_obj+base_infeasibility * self.sigma

                output_after_normal_step = self.get_model_simulation_result([input_after_normal_step])[0]
                obj_after_normal_step = output_after_normal_step[obj_var_name] / \
                        self.problem_description.scaling_factors[obj_var_name]

                if (base_infeasibility-infeasibility < self.options['feasibility_tol']):
                    # in this case, the feasibility progress is insignificant
                    pass
                elif (base_merit - merit)/(base_infeasibility-infeasibility)/self.sigma < self.options["sigma_kappa"]:
                    # in this case, if sigma increases,
                    # it must be that (base_obj-obj_after_normal_step)<0
                    self.sigma = max(self.sigma,(base_obj-obj_after_normal_step)/(self.options["sigma_kappa"]-1)/\
                                 (base_infeasibility - infeasibility)+1)
                    print("self.sigma=%f"%self.sigma)

            plant_trial_point_output = self.get_plant_simulation_result([filtered_input])[0]
            obj_var_name = self.problem_description.symbol_list['OBJ']
            obj = self.plant_history_data[self.iter_count][obj_var_name] / \
                  self.problem_description.scaling_factors[obj_var_name]
            self.plant_history_data[self.iter_count]['obj_for_merit'] = obj
            infeasibility_sq = 0
            for con_var_name in self.problem_description.symbol_list['CON']:
                infeasibility_sq += max(self.plant_history_data[self.iter_count][con_var_name], 0) ** 2 / \
                                    self.problem_description.scaling_factors[con_var_name] ** 2
            infeasibility = numpy.sqrt(infeasibility_sq)
            self.plant_history_data[self.iter_count]['infeasibility'] = infeasibility
            merit = obj + infeasibility * self.sigma
            self.plant_history_data[self.iter_count]['merit'] = merit
            base_merit = self.plant_history_data[base_iteration]['base_obj_for_merit'] +\
                         self.plant_history_data[base_iteration]['base_infeasibility'] * self.sigma
            self.plant_history_data[self.iter_count]['base_merit'] = base_merit

            model_trial_point_output = self.get_model_simulation_result([filtered_input])[0]
            print("model_trial_point_output")
            print(model_trial_point_output)
            obj_var_name = self.problem_description.symbol_list['OBJ']
            obj = self.model_history_data[self.iter_count][obj_var_name] / \
                  self.problem_description.scaling_factors[obj_var_name]
            self.model_history_data[self.iter_count]['obj_for_merit'] = obj
            infeasibility_sq = 0
            for con_var_name in self.problem_description.symbol_list['CON']:
                infeasibility_sq += max(self.model_history_data[self.iter_count][con_var_name], 0) ** 2 / \
                                    self.problem_description.scaling_factors[con_var_name] ** 2
            infeasibility = numpy.sqrt(infeasibility_sq)
            self.model_history_data[self.iter_count]['infeasibility'] = infeasibility
            merit = obj + infeasibility * self.sigma
            self.model_history_data[self.iter_count]['merit'] = merit
            base_merit = self.model_history_data[base_iteration]['base_obj_for_merit'] + \
                         self.model_history_data[base_iteration]['base_infeasibility'] * self.sigma
            self.model_history_data[self.iter_count]['base_merit'] = base_merit

            # accept the trial point
            flag_infeasible = False
            for con_name in self.problem_description.symbol_list['CON']:
                if plant_trial_point_output[con_name] > self.options['feasibility_tol']:
                    flag_infeasible = True
                    break
            if self.model_history_data[self.iter_count]['base_merit'] < self.model_history_data[self.iter_count]['merit']:
                # Because the solver is not a global optimization solver, it is possible
                # that the model merit function increases. In this case, we ask for backtracking.
                self.model_history_data[self.iter_count-1]['event'] = "model merit increases"
                print("model merit increases")
                rho = -2
                # if self.plant_history_data[base_iteration]['merit'] -\
                #                self.plant_history_data[self.iter_count]['merit'] >0:
                #     rho = 0.1
                # else:
                #     rho = -2
            else:
                if (not flag_infeasible) and\
                            abs(self.plant_history_data[self.iter_count]['base_merit'] - self.plant_history_data[self.iter_count]['merit']) < self.options[
                    'stationarity_tol']:
                    iter_successful_flag = True
                    self.model_history_data[self.iter_count-1]['event'] = "plant converges"
                    self.model_history_data[self.iter_count-1]['rho'] = 1
                    continue
                else:
                    if self.model_history_data[self.iter_count]['base_merit'] - self.model_history_data[self.iter_count]['merit'] >1e-6:
                        rho = (self.plant_history_data[self.iter_count]['base_merit'] - self.plant_history_data[self.iter_count]['merit']) / \
                              (self.model_history_data[self.iter_count]['base_merit'] - self.model_history_data[self.iter_count]['merit'])
                    else:
                        rho = (self.plant_history_data[self.iter_count]['base_merit'] -
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

    def project_to_trust_region(self, optimized_input, current_point, trust_radius):
        norm = self.length_of_the_trial_step(optimized_input, current_point)
        if norm < trust_radius:
            ret=optimized_input
        else:
            coeff = trust_radius/norm
            ret = {}
            for mv in current_point.keys():
                ret[mv]=current_point[mv]+(optimized_input[mv]-current_point[mv])*coeff
        return ret

    def length_of_the_trial_step(self, optimized_input, current_point):
        mvs = self.problem_description.symbol_list['MV']
        norm = 0
        for mv in mvs:
            norm += ((optimized_input[mv] - current_point[mv]) / self.problem_description.scaling_factors[mv]) ** 2
        norm = numpy.sqrt(norm)
        return norm


class MACompoStepTRBackupModel(ModifierAdaptationCompoStepTR):
    def set_problem(self, problem_description,
                    plant,
                    model,
                    backup_model,
                    perturbation_method,
                    noise_generator,
                    solver_executable,
                    spec_function,
                    modifier_type,
                    skipped_modifiers=[],
                    black_box_model=None):
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
        # TODO: FOR MODELS OTHER THAN PYOMO TYPE
        self.model_simulator = PyomoSimulator(PyomoModelWithModifiers(model, modifier_type, mvs, cvs))
        if black_box_model is None:
            self.model_optimizer = CompoStepTrustRegionOptimizer(PyomoModelWithModifiers(model, modifier_type, mvs, cvs))
        else:
            self.model_optimizer = CompoStepTrustRegionBBMOptimizer(BlackBoxModelWithModifiers(black_box_model, mvs, cvs))
        self.backup_model_simulator = PyomoSimulator(PyomoModelWithModifiers(backup_model, modifier_type, mvs, cvs))
        self.backup_model_optimizer = CompoStepTrustRegionOptimizer(PyomoModelWithModifiers(backup_model, modifier_type, mvs, cvs))
        self.perturbation_method = perturbation_method
        self.noise_generator = noise_generator
        self.solver_executable = solver_executable
        self.spec_function = spec_function
        self.skipped_modifiers = skipped_modifiers

    def initialize_simulation(self, starting_point, initial_parameter_value):
        ModifierAdaptationCompoStepTR.initialize_simulation(self, starting_point, initial_parameter_value)
        # build model
        self.backup_model_simulator.build(self.problem_description)
        self.backup_model_optimizer.build(self.problem_description)

        # set solver
        default_options={'max_iter':500,
                         "tol":1e-10}
        solver4 = SolverFactory('ipopt', executable=self.solver_executable)
        self.backup_model_simulator.set_solver(solver4, tee=False, default_options=default_options)
        solver5 = SolverFactory('ipopt', executable=self.solver_executable)
        self.backup_model_optimizer.set_solver(solver5, tee=False, default_options=default_options)

        self.set_initial_backup_model_parameters(initial_parameter_value)

        # set basepoint
        self.backup_model_simulator.set_base_point(self.current_point)
        self.backup_model_optimizer.set_base_point(self.current_point)

        # set trust radius for the backup problem
        self.trust_radius_backup = self.options["initial_trust_radius"]
    def register_option(self):
        super().register_option()
        self.available_options["kappa_b"] = ("[0,1]", 0.5)
        self.available_options["separate_tr_management"] = ("bool", True)
        self.available_options["skip_backup"] = ("bool", False)
        self.available_options["use_premature_solution"] = ("bool", True)

    def optimize_for_u(self, tr_radius, tr_base):
        print("primary model opt")
        return super().optimize_for_u(tr_radius, tr_base)

    def optimize_for_u_backup_model(self, tr_radius, tr_base):
        '''
        used in self.one_step_simulation
        :return:
        '''
        input_values = {}
        if self.spec_function is not None:
            for k, v in self.current_spec.items():
                input_values[k] = v
        print("backup model opt")
        optimized_input, input_after_normal_step, solve_status = self.backup_model_optimizer.optimize(input_values, tr_radius,
                                                                                               tr_base,
                                                                                               self.options['xi_N'],
                                                                                               param_values=self.current_backup_parameter_value,
                                                                                               use_homo=self.options[
                                                                                                   "homotopy_optimization"])
        return optimized_input, input_after_normal_step, solve_status
    def one_step_simulation(self):
        print("Iteration %d" % self.iter_count)
        flag_primary_model_fatal_error = False

        # initialize storage
        self.model_history_data[self.iter_count] = {}
        self.plant_history_data[self.iter_count] = {}
        self.input_history_data[self.iter_count] = {}
        if self.iter_count == 1:
            # If we want to add some print information, do it here.
            self.model_history_data[1]['rho'] = ""
            self.model_history_data[1]['rho_b'] = ""
            self.model_history_data[1]['tr'] = ""
            self.model_history_data[1]['tr_b'] = ""
            self.model_history_data[1]['event'] = ""
            obj_var_name = self.problem_description.symbol_list['OBJ']
            self.model_history_data[1]['base_'+obj_var_name] = ""
            for con_var_name in self.problem_description.symbol_list['CON']:
                self.model_history_data[1]['base_'+con_var_name] = ""

            self.model_history_data[1]['m_sigma'] = ""
            self.model_history_data[1]['m_merit'] = ""
            self.model_history_data[1]['m_infeasibility'] = ""
            self.model_history_data[1]['m_obj_for_merit'] = ""
            self.model_history_data[1]['m_base_merit'] = ""
            self.model_history_data[1]['m_base_infeasibility'] = ""
            self.model_history_data[1]['m_base_obj_for_merit'] = ""

            self.model_history_data[1]['b_sigma'] = ""
            self.model_history_data[1]['b_merit'] = ""
            self.model_history_data[1]['b_infeasibility'] = ""
            self.model_history_data[1]['b_obj_for_merit'] = ""
            self.model_history_data[1]['b_base_merit'] = ""
            self.model_history_data[1]['b_base_infeasibility'] = ""
            self.model_history_data[1]['b_base_obj_for_merit'] = ""

            self.model_history_data[1]['selected_model']=0
            self.model_history_data[1]['merit'] = ""
            self.model_history_data[1]['base_merit'] = ""

            self.plant_history_data[1]['merit'] = ""
            self.plant_history_data[1]['infeasibility'] = ""
            self.plant_history_data[1]['obj_for_merit'] = ""
            self.plant_history_data[1]['base_merit'] = ""
            self.plant_history_data[1]['base_infeasibility'] = ""
            self.plant_history_data[1]['base_obj_for_merit'] = ""

        # get trial point
        trial_points = self.get_trial_point()
        base_iteration = self.iter_count
        base_input = copy.deepcopy(trial_points[0])

        # get plant simulation result
        plant_output_data = self.get_plant_simulation_result(trial_points)
        obj_var_name = self.problem_description.symbol_list['OBJ']
        obj = self.plant_history_data[self.iter_count][obj_var_name] / \
              self.problem_description.scaling_factors[obj_var_name]
        self.plant_history_data[self.iter_count]['obj_for_merit'] = obj
        self.plant_history_data[base_iteration]['base_obj_for_merit'] = obj
        infeasibility_sq = 0
        for con_var_name in self.problem_description.symbol_list['CON']:
            infeasibility_sq += max(self.plant_history_data[self.iter_count][con_var_name], 0) ** 2 / \
                                self.problem_description.scaling_factors[con_var_name] ** 2
        infeasibility = numpy.sqrt(infeasibility_sq)
        self.plant_history_data[self.iter_count]['infeasibility'] = infeasibility
        self.plant_history_data[base_iteration]['base_infeasibility'] = infeasibility
        merit = obj + infeasibility * self.sigma
        self.plant_history_data[self.iter_count]['merit'] = merit
        self.plant_history_data[base_iteration]['base_merit'] = merit

        # get model simulation result with and without modifiers
        model_output_data = self.get_model_simulation_result(trial_points)
        if model_output_data is None:
            print("model simulation fails when updating modifiers")
            flag_primary_model_fatal_error=True
        else:
            obj_var_name = self.problem_description.symbol_list['OBJ']
            obj = self.model_history_data[self.iter_count][obj_var_name] / \
                  self.problem_description.scaling_factors[obj_var_name]
            self.model_history_data[self.iter_count]['m_obj_for_merit'] = obj
            infeasibility_sq = 0
            for con_var_name in self.problem_description.symbol_list['CON']:
                infeasibility_sq += max(self.model_history_data[self.iter_count][con_var_name], 0) ** 2 / \
                                    self.problem_description.scaling_factors[con_var_name] ** 2
            infeasibility = numpy.sqrt(infeasibility_sq)
            self.model_history_data[self.iter_count]['m_infeasibility'] = infeasibility
            merit = obj + infeasibility * self.sigma
            self.model_history_data[self.iter_count]['m_merit'] = merit

            # update modifiers
            self.update_modifiers(plant_output_data, model_output_data)

        # get backup model simulation result with and without modifiers
        model_output_data = self.get_backup_model_simulation_result(trial_points)
        obj_var_name = self.problem_description.symbol_list['OBJ']
        obj = self.model_history_data[self.iter_count]['b_'+obj_var_name] / \
              self.problem_description.scaling_factors[obj_var_name]
        self.model_history_data[self.iter_count]['b_obj_for_merit'] = obj
        infeasibility_sq = 0
        for con_var_name in self.problem_description.symbol_list['CON']:
            infeasibility_sq += max(self.model_history_data[self.iter_count]['b_'+con_var_name], 0) ** 2 / \
                                self.problem_description.scaling_factors[con_var_name] ** 2
        infeasibility = numpy.sqrt(infeasibility_sq)
        self.model_history_data[self.iter_count]['b_infeasibility'] = infeasibility
        merit = obj + infeasibility * self.sigma
        self.model_history_data[self.iter_count]['b_merit'] = merit

        # update modifiers for the backup model
        self.update_modifiers_backup_model(plant_output_data, model_output_data)

        # set basepoint
        if not flag_primary_model_fatal_error:
            self.model_simulator.set_base_point(self.current_point)
            self.model_optimizer.set_base_point(self.current_point)
        self.backup_model_simulator.set_base_point(self.current_point)
        self.backup_model_optimizer.set_base_point(self.current_point)

        if not flag_primary_model_fatal_error:
            self.store_model_adaptation_data()
        self.store_backup_model_adaptation_data()

        # culculate_base_merit
        try:
            base_model_output, solve_status = self.model_simulator.simulate(base_input,
                                                                            param_values=self.current_parameter_value,
                                                                            use_homo=self.options["homotopy_simulation"])
        except Exception as e:
            print(e)
            print("model simulation failure when calculating base merit")
            flag_primary_model_fatal_error = True
        else:
            if solve_status != ModelSolvingStatus.OK:
                print("model simulation failure when calculating base merit")
                flag_primary_model_fatal_error = True
        if not flag_primary_model_fatal_error:
            obj_var_name = self.problem_description.symbol_list['OBJ']
            base_obj = base_model_output[obj_var_name] / \
                       self.problem_description.scaling_factors[obj_var_name]
            infeasibility_sq = 0
            for con_var_name in self.problem_description.symbol_list['CON']:
                infeasibility_sq += max(base_model_output[con_var_name], 0) ** 2 / \
                                    self.problem_description.scaling_factors[con_var_name] ** 2
            base_infeasibility = numpy.sqrt(infeasibility_sq)

            self.model_history_data[base_iteration]['m_base_' + obj_var_name] = base_model_output[obj_var_name]
            self.model_history_data[base_iteration]['m_base_obj_for_merit'] = base_obj
            for con_var_name in self.problem_description.symbol_list['CON']:
                self.model_history_data[base_iteration]['m_base_' + con_var_name] = base_model_output[con_var_name]
            self.model_history_data[base_iteration]['m_base_infeasibility'] = base_infeasibility
            merit = base_obj + base_infeasibility * self.sigma
            self.model_history_data[base_iteration]['m_base_merit'] = merit

        # culculate_base_merit for the backup model
        base_model_output, solve_status = self.backup_model_simulator.simulate(base_input,
                                                                        param_values=self.current_backup_parameter_value,
                                                                        use_homo=self.options["homotopy_simulation"])
        obj_var_name = self.problem_description.symbol_list['OBJ']
        base_obj = base_model_output[obj_var_name] / \
                   self.problem_description.scaling_factors[obj_var_name]
        infeasibility_sq = 0
        for con_var_name in self.problem_description.symbol_list['CON']:
            infeasibility_sq += max(base_model_output[con_var_name], 0) ** 2 / \
                                self.problem_description.scaling_factors[con_var_name] ** 2
        base_infeasibility = numpy.sqrt(infeasibility_sq)

        self.model_history_data[base_iteration]['b_base_' + obj_var_name] = base_model_output[obj_var_name]
        self.model_history_data[base_iteration]['b_base_obj_for_merit'] = base_obj
        for con_var_name in self.problem_description.symbol_list['CON']:
            self.model_history_data[base_iteration]['b_base_' + con_var_name] = base_model_output[con_var_name]
        self.model_history_data[base_iteration]['b_base_infeasibility'] = base_infeasibility
        merit = base_obj + base_infeasibility * self.sigma
        self.model_history_data[base_iteration]['b_base_merit'] = merit

        iter_successful_flag = False
        rho = 100
        tr_base = {}
        for k, v in trial_points[0].items():
            if k in self.problem_description.symbol_list['MV']:
                tr_base[k] = v

        while (not iter_successful_flag):
            flag_primary_model_iter_calc_failure = flag_primary_model_fatal_error
            self.model_history_data[self.iter_count]['tr'] = self.trust_radius
            self.model_history_data[self.iter_count]['tr_b'] = self.trust_radius_backup
            # calculate the optimal input
            # tr_base equals self.current_point
            try:
                optimized_input_m, input_after_normal_step_m, solve_status_m = self.optimize_for_u(self.trust_radius, tr_base)
                if solve_status != ModelSolvingStatus.OK:
                    if self.options["use_premature_solution"]:
                        optimized_input_m = self.project_to_trust_region(optimized_input_m, self.current_point, self.trust_radius)
                        input_after_normal_step_m = self.project_to_trust_region(input_after_normal_step_m, self.current_point, self.trust_radius)
                        # from rtolib.util.misc import distance_of_two_dicts
                        # print("in optimization")
                        # print(self.trust_radius)
                        # raise Exception("for testing")
                        # print(distance_of_two_dicts(optimized_input_m, self.current_point))
                    else:
                        optimized_input_m = base_input  # self.current_point
                        input_after_normal_step_m = base_input  # self.current_point
                        self.model_history_data[self.iter_count]['event'] = "primary optimization failed"
                        print("primary model optimization failure")
                        flag_primary_model_iter_calc_failure = True
            except Exception as e:
                optimized_input_m = base_input #self.current_point
                input_after_normal_step_m = base_input #self.current_point
                self.model_history_data[self.iter_count]['event'] = "primary optimization failed"
                print("primary model optimization failure")
                flag_primary_model_iter_calc_failure = True
            # if solve_status_m == ModelSolvingStatus.OPTIMIZATION_FAILED:
            #     optimized_input_m = self.current_point
            #     input_after_normal_step_m = self.current_point
            #     self.model_history_data[self.iter_count]['event'] = "primary optimization failed"
            try:
                optimized_input_b, input_after_normal_step_b, solve_status_b = self.optimize_for_u_backup_model(self.trust_radius_backup, tr_base)
            except Exception as e:
                input_after_normal_step_b = base_input #self.current_point
                optimized_input_b = base_input #self.current_point
                self.model_history_data[self.iter_count]['event'] = "backup optimization failed"
            if solve_status_b == ModelSolvingStatus.OPTIMIZATION_FAILED:
                input_after_normal_step_b = base_input #self.current_point
                optimized_input_b = base_input #self.current_point
                self.model_history_data[self.iter_count]['event'] = "backup optimization failed"

            self.model_history_data[self.iter_count]['sigma'] = self.sigma
            # bound the input
            mv_bounds = self.problem_description.bounds
            filtered_input_m = self.adapt_to_bound_mv(optimized_input_m, mv_bounds)
            filtered_input_b = self.adapt_to_bound_mv(optimized_input_b, mv_bounds)

            # iter count
            if self.iter_count >= self.options['max_iter']+1:
                self.iter_count += 1
                break
            self.iter_count += 1
            self.model_history_data[self.iter_count] = {}
            self.plant_history_data[self.iter_count] = {}
            self.input_history_data[self.iter_count] = {}
            # print(self.input_history_data)

            if self.options['adaptive_sigma']:
                # check sigma is large enough
                # the primary model part
                # from rtolib.util.misc import distance_of_two_dicts
                # print("in adapting sigma")
                # print(distance_of_two_dicts(filtered_input_m, self.current_point))
                if not flag_primary_model_iter_calc_failure:
                    model_trial_point_output_m = self.get_model_simulation_result([filtered_input_m])[0]
                    if model_trial_point_output_m is None:
                        flag_primary_model_iter_calc_failure = False
                        print("primary model simulation failure in evaluating trial point")

                if not flag_primary_model_iter_calc_failure:
                    obj_var_name = self.problem_description.symbol_list['OBJ']
                    obj_after_optimization_m = self.model_history_data[self.iter_count][obj_var_name] / \
                                             self.problem_description.scaling_factors[obj_var_name]
                    infeasibility_sq_m = 0
                    for con_var_name in self.problem_description.symbol_list['CON']:
                        infeasibility_sq_m += max(self.model_history_data[self.iter_count][con_var_name], 0) ** 2 / \
                                            self.problem_description.scaling_factors[con_var_name] ** 2
                    infeasibility_m = numpy.sqrt(infeasibility_sq_m)
                    merit_m = obj_after_optimization_m + infeasibility_m * self.sigma
                    print("base_iteration: %d"%base_iteration)
                    base_infeasibility_m = self.model_history_data[base_iteration]['m_base_infeasibility']
                    base_obj_m = self.model_history_data[base_iteration]['m_base_obj_for_merit']
                    base_merit_m = base_obj_m + base_infeasibility_m * self.sigma

                if not flag_primary_model_iter_calc_failure:
                    output_after_normal_step_m = self.get_model_simulation_result([input_after_normal_step_m])[0]
                    if output_after_normal_step_m is None:
                        flag_primary_model_iter_calc_failure = False
                        print("primary model simulation failure in evaluating trial point")

                if not flag_primary_model_iter_calc_failure:
                    obj_after_normal_step_m = output_after_normal_step_m[obj_var_name] / \
                                            self.problem_description.scaling_factors[obj_var_name]

                if not flag_primary_model_iter_calc_failure:
                    if (base_infeasibility_m - infeasibility_m < self.options['feasibility_tol']):
                        # in this case, the feasibility progress is insignificant
                        sigma_m = self.sigma
                    elif (base_merit_m - merit_m) / (base_infeasibility_m - infeasibility_m) / self.sigma < self.options[
                        "sigma_kappa"]:
                        # in this case, if sigma increases,
                        # it must be that (base_obj-obj_after_normal_step)<0
                        sigma_m = max(self.sigma,
                                         (base_obj_m - obj_after_normal_step_m) / (self.options["sigma_kappa"] - 1) / \
                                         (base_infeasibility_m - infeasibility_m) + 1)
                        print("sigma_m=%f" % sigma_m)
                    else:
                        sigma_m = self.sigma
                else:
                    sigma_m = self.sigma

                # the backup model part
                model_trial_point_output_b = self.get_backup_model_simulation_result([filtered_input_b])[0]
                obj_var_name = self.problem_description.symbol_list['OBJ']
                obj_after_optimization_b = self.model_history_data[self.iter_count]['b_'+obj_var_name] / \
                                           self.problem_description.scaling_factors[obj_var_name]
                infeasibility_sq_b = 0
                for con_var_name in self.problem_description.symbol_list['CON']:
                    infeasibility_sq_b += max(self.model_history_data[self.iter_count]['b_'+con_var_name], 0) ** 2 / \
                                          self.problem_description.scaling_factors[con_var_name] ** 2
                infeasibility_b = numpy.sqrt(infeasibility_sq_b)
                merit_b = obj_after_optimization_b + infeasibility_b * self.sigma
                base_infeasibility_b = self.model_history_data[base_iteration]['b_base_infeasibility']
                base_obj_b = self.model_history_data[base_iteration]['b_base_obj_for_merit']
                base_merit_b = base_obj_b + base_infeasibility_b * self.sigma

                output_after_normal_step_b = self.get_backup_model_simulation_result([input_after_normal_step_b])[0]
                obj_after_normal_step_b = output_after_normal_step_b[obj_var_name] / \
                                          self.problem_description.scaling_factors[obj_var_name]

                if (base_infeasibility_b - infeasibility_b < self.options['feasibility_tol']):
                    # in this case, the feasibility progress is insignificant
                    sigma_b = self.sigma
                elif (base_merit_b - merit_b) / (base_infeasibility_b - infeasibility_b) / self.sigma < self.options[
                    "sigma_kappa"]:
                    # in this case, if sigma increases,
                    # it must be that (base_obj-obj_after_normal_step)<0
                    sigma_b = max(self.sigma,
                                  (base_obj_b - obj_after_normal_step_b) / (self.options["sigma_kappa"] - 1) / \
                                  (base_infeasibility_b - infeasibility_b) + 1)
                    print("sigma_b=%f" % sigma_b)
                else:
                    sigma_b = self.sigma
            else:
                sigma_m = self.sigma
                sigma_b = self.sigma

            # model selection
            # calculate relevent data from the primary model
            if not flag_primary_model_iter_calc_failure:
                model_trial_point_output_m = self.get_model_simulation_result([filtered_input_m])[0]
                if model_trial_point_output_m is None:
                    flag_primary_model_iter_calc_failure = False
                    print("primary model simulation failure in preparing data for model selection")
            if not flag_primary_model_iter_calc_failure:
                obj_var_name = self.problem_description.symbol_list['OBJ']
                obj_m = self.model_history_data[self.iter_count][obj_var_name] / \
                      self.problem_description.scaling_factors[obj_var_name]
                # print(base_input)
                # print(filtered_input_m)
                self.model_history_data[self.iter_count]['m_obj_for_merit'] = obj_m
                infeasibility_sq_m = 0
                for con_var_name in self.problem_description.symbol_list['CON']:
                    infeasibility_sq_m += max(self.model_history_data[self.iter_count][con_var_name], 0) ** 2 / \
                                        self.problem_description.scaling_factors[con_var_name] ** 2
                infeasibility_m = numpy.sqrt(infeasibility_sq_m)
                self.model_history_data[self.iter_count]['m_infeasibility'] = infeasibility_m
                merit_m = obj_m + infeasibility_m * sigma_m
                self.model_history_data[self.iter_count]['m_merit'] = merit_m
                base_merit_m = self.model_history_data[base_iteration]['m_base_obj_for_merit'] + \
                             self.model_history_data[base_iteration]['m_base_infeasibility'] * sigma_m
                self.model_history_data[self.iter_count]['m_base_merit'] = base_merit_m

            # calculate relevent data from the backup model
            model_trial_point_output_b = self.get_backup_model_simulation_result([filtered_input_b])[0]
            obj_var_name = self.problem_description.symbol_list['OBJ']
            obj_b = self.model_history_data[self.iter_count]['b_'+obj_var_name] / \
                    self.problem_description.scaling_factors[obj_var_name]
            self.model_history_data[self.iter_count]['b_obj_for_merit'] = obj_b
            infeasibility_sq_b = 0
            for con_var_name in self.problem_description.symbol_list['CON']:
                infeasibility_sq_b += max(self.model_history_data[self.iter_count]['b_'+con_var_name], 0) ** 2 / \
                                      self.problem_description.scaling_factors[con_var_name] ** 2
            infeasibility_b = numpy.sqrt(infeasibility_sq_b)
            self.model_history_data[self.iter_count]['b_infeasibility'] = infeasibility_b
            merit_b = obj_b + infeasibility_b * sigma_b
            self.model_history_data[self.iter_count]['b_merit'] = merit_b
            base_merit_b = self.model_history_data[base_iteration]['b_base_obj_for_merit'] + \
                           self.model_history_data[base_iteration]['b_base_infeasibility'] * sigma_b
            self.model_history_data[self.iter_count]['b_base_merit'] = base_merit_b

            # compare the two models
            selected = 'm'
            if flag_primary_model_iter_calc_failure:
                selected = 'b'
            else:
                c_improvement_m = self.model_history_data[base_iteration]['m_base_infeasibility']-\
                            infeasibility_m
                c_improvement_b = self.model_history_data[base_iteration]['b_base_infeasibility'] - \
                                  infeasibility_b
                f_improvement_m = base_merit_m - merit_m
                f_improvement_b = base_merit_b - merit_b
                print("c_improvement_m = %.6e"%c_improvement_m)
                print("c_improvement_b = %.6e"%c_improvement_b)
                print("f_improvement_m = %.6e"%f_improvement_m)
                print("f_improvement_b = %.6e"%f_improvement_b)
                # if self.iter_count > 20:
                    # print("here")
                if (c_improvement_m < self.options['feasibility_tol'] and c_improvement_b > \
                        self.options['feasibility_tol']) or \
                        (c_improvement_m < -self.options['feasibility_tol'] and c_improvement_b >0):
                    selected = 'b'
                else:
                    if c_improvement_m < c_improvement_b*self.options['kappa_b']:
                        selected='b'
                if (f_improvement_m < self.options['stationarity_tol'] and f_improvement_b > \
                        self.options['stationarity_tol']) or \
                        (f_improvement_m < -self.options['stationarity_tol'] and f_improvement_b >0):
                    selected = 'b'
                else:
                    if f_improvement_m/sigma_m < f_improvement_b/sigma_b*self.options['kappa_b']:
                        selected='b'
                if self.options['skip_backup']:
                    selected='m'
            if selected == 'm':
                filtered_input = filtered_input_m
                self.model_history_data[self.iter_count]['base_merit'] = \
                    self.model_history_data[self.iter_count]['m_base_merit']
                self.model_history_data[self.iter_count]['merit'] = \
                    self.model_history_data[self.iter_count]['m_merit']
                self.model_history_data[self.iter_count-1]['selected_model'] = 1
                self.sigma = sigma_m
            elif selected == 'b':
                filtered_input = filtered_input_b
                self.model_history_data[self.iter_count]['base_merit'] = \
                    self.model_history_data[self.iter_count]['b_base_merit']
                self.model_history_data[self.iter_count]['merit'] = \
                    self.model_history_data[self.iter_count]['b_merit']
                self.model_history_data[self.iter_count-1]['selected_model'] = 2
                self.sigma = sigma_b

            # set specification and store input data
            len_trial_step = self.length_of_the_trial_step(filtered_input, self.current_point)
            self.current_point = filtered_input
            for k, v in self.current_point.items():
                self.input_history_data[self.iter_count-1][k] = v

            if self.spec_function is not None:
                self.current_spec = self.spec_function.__call__(self.iter_count)
                for k, v in self.current_spec.items():
                    self.input_history_data[self.iter_count-1][k] = v

            # implement the trial point
            plant_trial_point_output = self.get_plant_simulation_result([filtered_input])[0]
            obj_var_name = self.problem_description.symbol_list['OBJ']
            obj = self.plant_history_data[self.iter_count][obj_var_name] / \
                  self.problem_description.scaling_factors[obj_var_name]
            self.plant_history_data[self.iter_count]['obj_for_merit'] = obj
            infeasibility_sq = 0
            for con_var_name in self.problem_description.symbol_list['CON']:
                infeasibility_sq += max(self.plant_history_data[self.iter_count][con_var_name], 0) ** 2 / \
                                    self.problem_description.scaling_factors[con_var_name] ** 2
            infeasibility = numpy.sqrt(infeasibility_sq)
            self.plant_history_data[self.iter_count]['infeasibility'] = infeasibility
            merit = obj + infeasibility * self.sigma
            self.plant_history_data[self.iter_count]['merit'] = merit
            base_merit = self.plant_history_data[base_iteration]['base_obj_for_merit'] + \
                         self.plant_history_data[base_iteration]['base_infeasibility'] * self.sigma
            self.plant_history_data[self.iter_count]['base_merit'] = base_merit

            if selected == 'm':
                # calculate the merit function for the two models using the new sigma
                model_trial_point_output_m = self.get_model_simulation_result([filtered_input])[0]
                if model_trial_point_output_m is None:
                    flag_primary_model_iter_calc_failure = False
                    print("primary model simulation failure in preparing data for calculating rho")
                if not flag_primary_model_iter_calc_failure:
                    obj_var_name = self.problem_description.symbol_list['OBJ']
                    obj_m = self.model_history_data[self.iter_count][obj_var_name] / \
                            self.problem_description.scaling_factors[obj_var_name]
                    infeasibility_sq_m = 0
                    for con_var_name in self.problem_description.symbol_list['CON']:
                        infeasibility_sq_m += max(self.model_history_data[self.iter_count][con_var_name], 0) ** 2 / \
                                              self.problem_description.scaling_factors[con_var_name] ** 2
                    infeasibility_m = numpy.sqrt(infeasibility_sq_m)
                    merit_m = obj_m + infeasibility_m * self.sigma
                    self.model_history_data[self.iter_count]['merit'] = merit_m
                    self.model_history_data[self.iter_count]['m_merit'] = merit_m
                    base_merit_m = self.model_history_data[base_iteration]['m_base_obj_for_merit'] + \
                                   self.model_history_data[base_iteration]['m_base_infeasibility'] * self.sigma
                    self.model_history_data[self.iter_count]['base_merit'] = base_merit_m
                    self.model_history_data[self.iter_count]['m_base_merit'] = base_merit_m
            elif selected == 'b':
                # calculate relevent data from the backup model
                model_trial_point_output_b = self.get_backup_model_simulation_result([filtered_input])[0]
                obj_var_name = self.problem_description.symbol_list['OBJ']
                obj_b = self.model_history_data[self.iter_count]['b_' + obj_var_name] / \
                        self.problem_description.scaling_factors[obj_var_name]
                infeasibility_sq_b = 0
                for con_var_name in self.problem_description.symbol_list['CON']:
                    infeasibility_sq_b += max(self.model_history_data[self.iter_count]['b_' + con_var_name], 0) ** 2 / \
                                          self.problem_description.scaling_factors[con_var_name] ** 2
                infeasibility_b = numpy.sqrt(infeasibility_sq_b)
                merit_b = obj_b + infeasibility_b * self.sigma
                self.model_history_data[self.iter_count]['b_merit'] = merit_b
                self.model_history_data[self.iter_count]['merit'] = merit_b
                base_merit_b = self.model_history_data[base_iteration]['b_base_obj_for_merit'] + \
                               self.model_history_data[base_iteration]['b_base_infeasibility'] * self.sigma
                self.model_history_data[self.iter_count]['b_base_merit'] = base_merit_b
                self.model_history_data[self.iter_count]['base_merit'] = base_merit_b

            # accept the trial point
            flag_infeasible = False
            for con_name in self.problem_description.symbol_list['CON']:
                if plant_trial_point_output[con_name] > self.options['feasibility_tol']:
                    flag_infeasible = True
                    break

            if not self.options["separate_tr_management"]:
                if self.model_history_data[self.iter_count]['base_merit'] < self.model_history_data[self.iter_count][
                    'merit']:
                    self.model_history_data[self.iter_count - 1]['event'] = "model merit increases"
                    rho = -2
                else:
                    if (not flag_infeasible) and \
                            abs(self.plant_history_data[self.iter_count]['base_merit'] -
                                self.plant_history_data[self.iter_count]['merit']) < self.options[
                        'stationarity_tol']:
                        iter_successful_flag = True
                        self.model_history_data[self.iter_count - 1]['event'] = "plant converges"
                        self.model_history_data[self.iter_count - 1]['rho'] = 1
                        continue
                    else:
                        if self.model_history_data[self.iter_count]['base_merit'] - \
                                self.model_history_data[self.iter_count]['merit'] > 1e-6:
                            rho = (self.plant_history_data[self.iter_count]['base_merit'] -
                                   self.plant_history_data[self.iter_count]['merit']) / \
                                  (self.model_history_data[self.iter_count]['base_merit'] -
                                   self.model_history_data[self.iter_count]['merit'])
                        else:
                            rho = (self.plant_history_data[self.iter_count]['base_merit'] -
                                   self.plant_history_data[self.iter_count]['merit']) / 1e-6
                self.model_history_data[self.iter_count - 1]['rho'] = rho

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
                self.trust_radius_backup = self.trust_radius
            else:
                if self.model_history_data[self.iter_count]['base_merit'] < self.model_history_data[self.iter_count][
                    'merit']:
                    self.model_history_data[self.iter_count - 1]['event'] = "final merit increases"
                    rho = -2
                else:
                    if (not flag_infeasible) and \
                            abs(self.plant_history_data[self.iter_count]['base_merit'] -
                                self.plant_history_data[self.iter_count]['merit']) < self.options[
                        'stationarity_tol']:
                        iter_successful_flag = True
                        self.model_history_data[self.iter_count - 1]['event'] = "plant converges"
                        self.model_history_data[self.iter_count - 1]['rho'] = 1
                        continue
                # calculate rho for the primary model optimization problem
                if not flag_primary_model_iter_calc_failure:
                    model_improvement = self.model_history_data[self.iter_count]['m_base_merit'] - \
                                self.model_history_data[self.iter_count]['m_merit']
                    if model_improvement >= 0 and model_improvement < 1e-4:
                        model_improvement = 1e-4
                    elif model_improvement <= 0 and model_improvement > -1e-4:
                        model_improvement = -1e-4
                    rho_m = (self.plant_history_data[self.iter_count]['base_merit'] -
                                   self.plant_history_data[self.iter_count]['merit']) / model_improvement
                    self.model_history_data[self.iter_count - 1]['rho'] = rho_m
                else:
                    rho_m = -3
                    self.model_history_data[self.iter_count - 1]['rho'] = rho_m

                # calculate rho for the backup model optimization problem
                model_improvement = self.model_history_data[self.iter_count]['b_base_merit'] - \
                                    self.model_history_data[self.iter_count]['b_merit']
                if model_improvement >= 0 and model_improvement < 1e-4:
                    model_improvement = 1e-4
                elif model_improvement <= 0 and model_improvement > -1e-4:
                    model_improvement = -1e-4
                rho_b = (self.plant_history_data[self.iter_count]['base_merit'] -
                         self.plant_history_data[self.iter_count]['merit']) / model_improvement
                self.model_history_data[self.iter_count - 1]['rho_b'] = rho_b

                if selected == 'm':
                    rho = rho_m
                elif selected == 'b':
                    rho = rho_b

                # update trust-region radius
                if not flag_primary_model_iter_calc_failure:
                    if rho_m < self.options["eta1"] and self.trust_radius >= len_trial_step-1e-5:
                        gamma_m = self.options['gamma1']
                    elif rho_m > self.options["eta2"] and (self.trust_radius <= len_trial_step+1e-5 or selected == 'm'):
                        gamma_m = self.options['gamma3']
                    else:
                        gamma_m = 1
                    self.trust_radius *= gamma_m

                if rho_b < self.options["eta1"] and self.trust_radius_backup >= len_trial_step-1e-5:
                    gamma_b = self.options['gamma1']
                elif rho_b > self.options["eta2"] and (self.trust_radius_backup <= len_trial_step+1e-5 or selected == 'b'):
                    gamma_b = self.options['gamma3']
                else:
                    gamma_b = 1
                self.trust_radius_backup *= gamma_b

                if self.trust_radius < 1e-8:
                    self.trust_radius = 1e-8
                    self.model_history_data[self.iter_count - 1]['event'] = "model trust region too small"
                if self.trust_radius > self.options['max_trust_radius']:
                    self.trust_radius = self.options['max_trust_radius']
                    self.model_history_data[self.iter_count - 1]['event'] = "model trust region too large"
                if self.trust_radius_backup < 1e-8:
                    self.trust_radius_backup = 1e-8
                    self.model_history_data[self.iter_count - 1]['event'] = "backup trust region too small"
                if self.trust_radius_backup > self.options['max_trust_radius']:
                    self.trust_radius_backup = self.options['max_trust_radius']
                    self.model_history_data[self.iter_count - 1]['event'] = "backup trust region too large"

                # set successful flag
                if rho >= self.options["eta1"]:
                    iter_successful_flag = True