#!/usr/bin/env python
#-*- coding:utf-8 -*-
from rtolib.core.algo.common import PE_type_Algorithm, PE_type_Algorithm
from pyomo.environ import SolverFactory
from rtolib.core.pyomo_model import *
from rtolib.core.solve import PyomoSimulator,PyomoOptimizer,PyomoGradientParamEstimator

class GeneralizedParameterEstimation(PE_type_Algorithm):
    def register_option(self):
        super(GeneralizedParameterEstimation, self).register_option()
        self.available_options["prior_theta_strategy"]=("[Fixed,Adapted]", "Fixed")

    def set_homotopy_var(self, homotopy_var):
        '''

        :param homotopy_var: list of names
        :return:
        '''
        self.plant_simulator.homotopy_var = homotopy_var
        self.model_simulator.homotopy_var = homotopy_var
        self.model_optimizer.homotopy_var = homotopy_var
        self.pe_estimator.homotopy_var = homotopy_var

    def set_problem(self, problem_description,
                    plant,
                    model,
                    perturbation_method,
                    noise_generator,
                    solver_executable,
                    spec_function,
                    modifier_type,
                    output_weight=None,
                    parameter_weight=None,
                    parameter_initial_guess=None,
                    fixed_parameter_values=None,
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
        elif modifier_type == ModifierType:
            cvs = [self.problem_description.symbol_list['OBJ']]
            for cv in self.problem_description.symbol_list['CON']:
                cvs.append(cv)
            for cv in self.problem_description.symbol_list['CV']:
                cvs.append(cv)
        self.model_simulator = PyomoSimulator(PyomoModelWithModifiers(model,modifier_type, mvs, cvs))
        self.model_optimizer = PyomoOptimizer(PyomoModelWithModifiers(model,modifier_type, mvs, cvs))
        self.pe_estimator = PyomoGradientParamEstimator(PyomoModelWithModifiers(model,modifier_type, mvs, cvs),\
                                                        perturbation_method.number_of_data_points,
                                                        use_obj_and_con_info=True)
        self.perturbation_method = perturbation_method
        self.noise_generator = noise_generator
        self.solver_executable = solver_executable
        self.spec_function = spec_function

        self.set_weight(output_weight, parameter_weight, parameter_initial_guess,
                        fixed_parameter_values)


    def initialize_simulation(self, starting_point, initial_parameter_value):
        # build model
        self.plant_simulator.build(self.problem_description)
        self.model_simulator.build(self.problem_description)
        self.model_optimizer.build(self.problem_description)
        self.pe_estimator.build(self.problem_description)

        # set solver
        # TODO: solver parameter tuning
        default_options={'max_iter':100}
        solver1 = SolverFactory('ipopt', executable=self.solver_executable)
        self.plant_simulator.set_solver(solver1, tee=False, default_options=default_options)
        solver2 = SolverFactory('ipopt', executable=self.solver_executable)
        self.model_simulator.set_solver(solver2, tee=False, default_options=default_options)
        solver3 = SolverFactory('ipopt', executable=self.solver_executable)
        self.model_optimizer.set_solver(solver3, tee=False, default_options=default_options)
        solver4 = SolverFactory('ipopt', executable=self.solver_executable)
        pe_options={'max_iter':100,
                    }
        self.pe_estimator.set_solver(solver4, tee=False, default_options=pe_options)

        self.iter_count = 0
        self.input_history_data[0] = {}
        self.set_current_point(starting_point)

        self.iter_count = 1

        self.model_history_data[0] = {}
        self.set_initial_model_parameters(initial_parameter_value)

        # set basepoint
        self.model_simulator.set_base_point(self.current_point)
        self.model_optimizer.set_base_point(self.current_point)
        self.pe_estimator.set_base_point(self.current_point)


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

        self.current_point_model_simulation(trial_points[0])

        # set basepoint
        self.model_simulator.set_base_point(self.current_point)
        self.model_optimizer.set_base_point(self.current_point)
        self.pe_estimator.set_base_point(self.current_point)

        # parameter estimation
        self.current_parameter_value = self.pe_estimator.estimate_parameter(plant_data,
                                             fixed_param_values=self.fixed_parameter_values,
                                             use_homo=self.options["homotopy_simulation"],
                                             pre_simulation=self.options["pre-simulation_before_pe"])

        # set prior theta
        if self.options["prior_theta_strategy"] == "Adapted":
            self.pe_estimator.set_parameter_guess(self.current_parameter_value)

        self.store_model_adaptation_data()

        optimized_input, solve_status = self.optimize_for_u()

        mv_bounds = self.problem_description.bounds
        filtered_input = self.filter_mv(optimized_input, mv_bounds)

        # set specification and store input data
        self.set_current_point(filtered_input)

        # iter count
        self.iter_count += 1
