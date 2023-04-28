#!/usr/bin/env python
#-*- coding:utf-8 -*-
import numpy
from rtolib.core.algo.dccpwl_ma import DCCPWL_ModifierAdaptation
from rtolib.core.dc_cpwl_model import QuadraticBoostedDCCPWL_RTOObjectSubgrad
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

