#!/usr/bin/env python
#-*- coding:utf-8 -*-
from rtolib.core.algorithm import Algorithm

class PE_type_Algorithm(Algorithm):
    def set_initial_model_parameters(self, set_pyomo_model_parameters):
        super().set_initial_model_parameters(set_pyomo_model_parameters)

        # self.pe_estimator.model.display("result/display.txt")
        self.pe_estimator.set_weight(self.output_weight, self.parameter_weight)
        self.pe_estimator.set_parameter_guess(self.parameter_initial_guess)

    def set_weight(self, output_weight, parameter_weight, parameter_initial_guess,
                    fixed_parameter_values):
        if output_weight is None:
            output_weight = {}
            for cv in self.problem_description.symbol_list['CV']:
                output_weight[cv] = 1 / self.problem_description.scaling_factors[cv]
        self.output_weight = output_weight
        if parameter_weight is None:
            parameter_weight = {}
        self.parameter_weight = parameter_weight
        if parameter_initial_guess is None:
            parameter_initial_guess = {}
        self.parameter_initial_guess = parameter_initial_guess
        if fixed_parameter_values is None:
            fixed_parameter_values = {}
        self.fixed_parameter_values = fixed_parameter_values


class MA_type_Algorithm(Algorithm):

    def get_model_simulation_result(self, trial_points):
        model_output_data = [{} for i in range(len(trial_points))]
        for i, p in enumerate(trial_points):
            outputs, solve_status = self.model_simulator.simulate(p, param_values=self.current_parameter_value,
                                                                  use_homo=self.options["homotopy_simulation"])
            # TODO:deal with solve status
            if i == 0:
                for k, v in outputs.items():
                    self.model_history_data[self.iter_count][k] = v
            for k in self.problem_description.symbol_list['CV']:
                if k + '_unmodified' in outputs.keys():
                    model_output_data[i][k] = outputs[k + '_unmodified']
        return model_output_data

    def update_modifiers(self, plant_output_data, model_output_data):
        modifiers = self.perturbation_method.calculate_modifiers(plant_output_data, model_output_data)
        for k, v in modifiers.items():
            if k[1] is None:
                if k[0]+"_eps" not in self.skipped_modifiers:
                    self.current_parameter_value[k[0] + "_eps"] = v
            else:
                if k[1] + "_" + k[0] + "_lam" not in self.skipped_modifiers:
                    self.current_parameter_value[k[1] + "_" + k[0] + "_lam"] = v



