#!/usr/bin/env python
#-*- coding:utf-8 -*-
from rtolib.core.algo.common import MA_type_Algorithm
from rtolib.core.pyomo_model import *
from rtolib.core.solve import PyomoSimulator,PyomoOptimizer

class ModifierAdaptation(MA_type_Algorithm):

    def register_option(self):
        self.available_options["filtering_factor"] = ("[0,1]",0.5)
        self.available_options["homotopy_simulation"] = ("bool",True)
        self.available_options["homotopy_optimization"] = ("bool",True)

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
        default_options={'max_iter':100}
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




