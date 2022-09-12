#!/usr/bin/env python
#-*- coding:utf-8 -*-
from rtolib.core.algo.gpe import GeneralizedParameterEstimation
from pyomo.environ import SolverFactory
from rtolib.core.pyomo_model import *
from rtolib.core.solve import PyomoSimulator,PyomoOptimizer,PyomoGradientParamEstimator
import copy

class GPEAdaptedWeight(GeneralizedParameterEstimation):
    def initialize_simulation(self, starting_point, initial_parameter_value):
        super(GPEAdaptedWeight, self).initialize_simulation(starting_point, initial_parameter_value)

        self.candidate_pe_results=[]
        self.candidate_para_weights=[]


    def prepare_weight_candidates(self):
        self.candidate_para_weights=[]
        factor=2
        self.candidate_para_weights.append(copy.deepcopy(self.parameter_weight))
        self.candidate_para_weights.append(copy.deepcopy(self.parameter_weight))
        for modifier in self.model_simulator.pyomo_model.modifier_names_iterator():
            self.candidate_para_weights[1][modifier] *= factor
        self.candidate_para_weights.append(copy.deepcopy(self.parameter_weight))
        for modifier in self.model_simulator.pyomo_model.modifier_names_iterator():
            self.candidate_para_weights[2][modifier] /= factor
        return self.candidate_para_weights


    def one_step_simulation(self):
        print("Iteration %d"%self.iter_count)

        # initialize storage
        self.model_history_data[self.iter_count] = {}
        self.plant_history_data[self.iter_count] = {}
        self.input_history_data[self.iter_count] = {}

        # get trial point
        trial_points = self.get_trial_point()

        # get plant simulation result

        plant_data = self.get_plant_simulation_result(trial_points)
        for i, p in enumerate(trial_points):
            for k,v in p.items():
                plant_data[i][k]=v

        candidates_error=[]
        len_trial_point = len(trial_points)
        if len(self.candidate_pe_results) >0:
            for xi_candidates in self.candidate_pe_results:
                self.current_parameter_value = xi_candidates
                model_data = self.get_model_simulation_result(trial_points)
                error=0
                for meas in self.pe_estimator.all_available_measurement:
                    error += self.output_weight[meas] * (plant_data[0][meas]-model_data[0][meas])**2
                    for case in range(len_trial_point-1):
                        error += self.output_weight[meas] * (model_data[case+1][meas]-model_data[0][meas]-\
                                                             plant_data[case+1][meas]+plant_data[0][meas])**2/2
                candidates_error.append(error)
            smallest_error_index = candidates_error.index(min(candidates_error))
            self.parameter_weight = self.candidate_para_weights[smallest_error_index]
            self.model_history_data[self.iter_count]['smallest_error_index']=smallest_error_index
        else:
            self.model_history_data[self.iter_count]['smallest_error_index'] = 0

        self.current_point_model_simulation(trial_points[0])

        # set basepoint
        self.model_simulator.set_base_point(self.current_point)
        self.model_optimizer.set_base_point(self.current_point)
        self.pe_estimator.set_base_point(self.current_point)

        # parameter estimation
        self.prepare_weight_candidates()
        self.candidate_pe_results=[]
        for w_candidates in self.candidate_para_weights:
            self.parameter_weight = w_candidates
            self.pe_estimator.set_weight(self.output_weight, self.parameter_weight)
            self.candidate_pe_results.append(self.pe_estimator.estimate_parameter(plant_data,
                                             fixed_param_values=self.fixed_parameter_values,
                                             use_homo=self.options["homotopy_simulation"],
                                             pre_simulation=self.options["pre-simulation_before_pe"]))
        self.current_parameter_value = self.candidate_pe_results[0]

        # set prior theta
        if self.options["prior_theta_strategy"] == "Adapted":
            self.pe_estimator.set_parameter_guess(self.candidate_pe_results[0])

        self.store_model_adaptation_data()

        optimized_input, solve_status = self.optimize_for_u()

        mv_bounds = self.problem_description.bounds
        filtered_input = self.filter_mv(optimized_input, mv_bounds)

        # set specification and store input data
        self.set_current_point(filtered_input)

        # iter count
        self.iter_count += 1
