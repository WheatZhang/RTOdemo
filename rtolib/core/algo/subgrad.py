#!/usr/bin/env python
#-*- coding:utf-8 -*-
import numpy
from rtolib.core.algo.dccpwl_ma import DCCPWL_ModifierAdaptation
from rtolib.core.dc_cpwl_model import QuadraticBoostedDCCPWL_RTOObjectSubgrad,\
        QuadraticBoostedDCCPWL_PenaltyTR_Object
from rtolib.core.solve import PyomoSimulator
import copy
from rtolib.core import ModifierType
from pyomo.environ import SolverFactory


class DCCPWL_ModifierAdaptationSubgrad(DCCPWL_ModifierAdaptation):
    def set_problem(self, problem_description,
                    plant,
                    model_dc_cpwl_functions,
                    perturbation_method,
                    noise_generator,
                    nlp_solver_executable,
                    qcqp_solver_executable,
                    spec_function,
                    modifier_type,
                    parameter_set,
                    ):
        self.problem_description = problem_description
        self.plant_simulator = PyomoSimulator(plant)
        mvs = self.problem_description.symbol_list['MV']
        if modifier_type == ModifierType.RTO:
            cvs = [self.problem_description.symbol_list['OBJ']]
            for cv in self.problem_description.symbol_list['CON']:
                cvs.append(cv)
            cvs.append("validity_con")
        else:
            raise ValueError("Not applicable")
        self.DC_CPWL_RTO_model = QuadraticBoostedDCCPWL_RTOObjectSubgrad(model_dc_cpwl_functions, mvs, cvs, [])
        bounds = {}
        for k, v in self.problem_description.bounds.items():
            if k in mvs:
                bounds[k] = v
        self.DC_CPWL_RTO_model.set_input_bounds(bounds)
        self.perturbation_method = perturbation_method
        self.noise_generator = noise_generator
        self.nlp_solver_executable = nlp_solver_executable
        self.qcqp_solver_executable = qcqp_solver_executable
        self.spec_function = spec_function
        self.parameter_set = parameter_set

    def update_modifiers(self, plant_output_data):
        print(plant_output_data)
        plant_y_and_k = self.perturbation_method.calculate_plant_y_and_k(plant_output_data)
        self.DC_CPWL_RTO_model.update_modifiers(plant_y_and_k, self.current_point)


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
        self.update_modifiers(plant_output_data)

        self.store_model_adaptation_data()

        optimized_input, solve_status = self.optimize_for_u()

        mv_bounds = self.problem_description.bounds
        filtered_input = self.filter_mv(optimized_input, mv_bounds)

        # set specification and store input data
        self.set_current_point(filtered_input)

        # iter count
        self.iter_count += 1

class DCCPWL_ModifierAdaptationTRPenalty(DCCPWL_ModifierAdaptationSubgrad):
    def set_problem(self, problem_description,
                    plant,
                    model_dc_cpwl_functions,
                    perturbation_method,
                    noise_generator,
                    nlp_solver_executable,
                    qcqp_solver_executable,
                    spec_function,
                    modifier_type,
                    parameter_set,
                    model_mvs=None,
                    ):
        self.problem_description = problem_description
        self.plant_simulator = PyomoSimulator(plant)
        if model_mvs is None:
            mvs = self.problem_description.symbol_list['MV']
        else:
            mvs = model_mvs
        if modifier_type == ModifierType.RTO:
            cvs = [self.problem_description.symbol_list['OBJ']]
            for cv in self.problem_description.symbol_list['CON']:
                cvs.append(cv)
            cvs.append("validity_con")
        else:
            raise ValueError("Not applicable")
        self.DC_CPWL_RTO_model = QuadraticBoostedDCCPWL_PenaltyTR_Object(model_dc_cpwl_functions, mvs, cvs, [])
        bounds = {}
        for k, v in self.problem_description.bounds.items():
            if k in mvs:
                bounds[k] = v
        self.DC_CPWL_RTO_model.set_input_bounds(bounds)
        self.perturbation_method = perturbation_method
        self.noise_generator = noise_generator
        self.nlp_solver_executable = nlp_solver_executable
        self.qcqp_solver_executable = qcqp_solver_executable
        self.spec_function = spec_function
        self.parameter_set = parameter_set

    def register_option(self):
        super().register_option()
        self.available_options["sigma"] = ("positive float", 10)

    def initialize_simulation(self, starting_point, initial_parameter_value):
        super().initialize_simulation(starting_point, initial_parameter_value)
        self.sigma = self.options['sigma']
        self.model_optimizer.set_penalty_coeff(self.sigma)

    def optimize_for_u(self, tr_radius, tr_base):
        spec_values = {}
        if self.spec_function is not None:
            for k, v in self.current_spec.items():
                spec_values[k] = v
        optimized_input, solve_status = self.DC_CPWL_RTO_model.optimize(spec_values, self.current_point, tr_radius, tr_base)
        # TODO:deal with solve status, fallback strategy
        return optimized_input, solve_status

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
            self.model_history_data[1]['base_' + obj_var_name] = ""
            for con_var_name in self.problem_description.symbol_list['CON']:
                self.model_history_data[1]['base_' + con_var_name] = ""
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
            infeasibility_sq += max(self.plant_history_data[self.iter_count][con_var_name], 0) ** 2 / \
                                self.problem_description.scaling_factors[con_var_name] ** 2
        ret += numpy.sqrt(infeasibility_sq) * self.sigma
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
        self.update_modifiers(plant_output_data)

        self.store_model_adaptation_data()

        # culculate_base_merit
        base_model_output, solve_status = self.DC_CPWL_RTO_model.simulate(base_input)
        obj_var_name = self.problem_description.symbol_list['OBJ']
        ret = base_model_output[obj_var_name] / \
              self.problem_description.scaling_factors[obj_var_name]
        infeasibility_sq = 0
        for con_var_name in self.problem_description.symbol_list['CON']:
            infeasibility_sq += max(base_model_output[con_var_name], 0) ** 2 / \
                                self.problem_description.scaling_factors[con_var_name] ** 2
        ret += numpy.sqrt(infeasibility_sq) * self.sigma
        base_merit = ret

        self.model_history_data[self.iter_count]['base_' + obj_var_name] = base_model_output[obj_var_name]
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
                self.model_history_data[self.iter_count - 1]['event'] = "model merit increases"
                rho = -2
                # if self.plant_history_data[base_iteration]['merit'] -\
                #                self.plant_history_data[self.iter_count]['merit'] >0:
                #     rho = 0.1
                # else:
                #     rho = -2
            else:
                if (not flag_infeasible) and \
                        abs(self.plant_history_data[base_iteration]['merit'] - self.plant_history_data[self.iter_count][
                            'merit']) < self.options[
                    'stationarity_tol']:
                    iter_successful_flag = True
                    self.model_history_data[self.iter_count - 1]['event'] = "plant converges"
                    self.model_history_data[self.iter_count - 1]['rho'] = 1
                    continue
                else:
                    if base_merit - self.model_history_data[self.iter_count]['merit'] > 1e-6:
                        rho = (self.plant_history_data[base_iteration]['merit'] -
                               self.plant_history_data[self.iter_count]['merit']) / \
                              (base_merit - self.model_history_data[self.iter_count]['merit'])
                    else:
                        rho = (self.plant_history_data[base_iteration]['merit'] -
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
