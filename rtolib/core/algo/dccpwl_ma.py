#!/usr/bin/env python
#-*- coding:utf-8 -*-
import numpy
from rtolib.core.algo.ma import ModifierAdaptation
from rtolib.core.dc_cpwl_model import QuadraticBoostedDCCPWL_RTOObject,\
        QuadraticBoostedDCCPWL_PenaltyTR_Object
from rtolib.core.solve import PyomoSimulator
import copy
from rtolib.core import ModifierType
from pyomo.environ import SolverFactory


class DCCPWL_ModifierAdaptation(ModifierAdaptation):
    def __init__(self):
        self.model_history_data = {}
        self.plant_history_data = {}
        self.input_history_data = {}
        self.available_options = {}
        self.register_option()
        self.iter_count = -1
        self.spec_function=None
        self.plant_simulator=None
        self.DC_CPWL_RTO_model=None
        self.problem_description=None
        self.noise_generator=None
        self.perturbation_method=None

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
        self.DC_CPWL_RTO_model = QuadraticBoostedDCCPWL_RTOObject(model_dc_cpwl_functions, mvs, cvs, [])
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

    def set_initial_model_parameters(self, initial_parameter_value):
        # set parameter
        self.current_parameter_value = {}
        for p in initial_parameter_value.keys():
            self.current_parameter_value[p] = initial_parameter_value[p]
        for k, v in self.current_parameter_value.items():
            self.model_history_data[0][k] = v

    def initialize_simulation(self, starting_point, initial_parameter_value):
        # build model
        self.plant_simulator.build(self.problem_description)

        default_options={'max_iter':500,
                         "tol":1e-10}
        solver1 = SolverFactory('ipopt', executable=self.nlp_solver_executable)
        self.plant_simulator.set_solver(solver1, tee=False, default_options=default_options)

        self.iter_count = 0
        self.input_history_data[0] = {}
        self.set_initial_model_parameters(initial_parameter_value)
        self.set_current_point(starting_point)

        self.iter_count = 1

        self.model_history_data[0] = {}

        self.DC_CPWL_RTO_model.build(self.problem_description)
        solver2 = SolverFactory('cplex', executable=self.qcqp_solver_executable)
        self.DC_CPWL_RTO_model.set_solver(solver2, tee=False,\
                                          default_options={'logfile':"qcqp.log"})

        # solver2 = SolverFactory('ipopt', executable=self.nlp_solver_executable)
        # self.DC_CPWL_RTO_model.set_solver(solver2, tee=False, default_options=default_options)


    def get_model_simulation_result(self, trial_points, with_modifier=False):
        model_output_data = [{} for i in range(len(trial_points))]
        for i, p in enumerate(trial_points):
            outputs, solve_status = self.DC_CPWL_RTO_model.simulate(p, with_modifier)
            # TODO:deal with solve status
            if i == 0:
                for k, v in outputs.items():
                    self.model_history_data[self.iter_count][k] = v

            for k in self.available_measurements():
                model_output_data[i][k] = outputs[k]
        return model_output_data

    def update_modifiers(self, plant_output_data, model_output_data):
        print(plant_output_data)
        print(model_output_data)
        modifiers = self.perturbation_method.calculate_modifiers(plant_output_data, model_output_data)
        self.DC_CPWL_RTO_model.update_modifiers(modifiers, self.current_point)
        for k, v in modifiers.items():
            if k[1] is None:
                self.current_parameter_value[k[0] + "_eps"] = v
            else:
                self.current_parameter_value[k[1] + "_" + k[0] + "_lam"] = v

    def optimize_for_u(self):
        spec_values = {}
        if self.spec_function is not None:
            for k, v in self.current_spec.items():
                spec_values[k] = v
        optimized_input, solve_status = self.DC_CPWL_RTO_model.optimize(spec_values, self.current_point)
        # TODO:deal with solve status, fallback strategy
        return optimized_input, solve_status


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

        self.store_model_adaptation_data()

        optimized_input, solve_status = self.optimize_for_u()

        mv_bounds = self.problem_description.bounds
        filtered_input = self.filter_mv(optimized_input, mv_bounds)

        # set specification and store input data
        self.set_current_point(filtered_input)

        # iter count
        self.iter_count += 1


