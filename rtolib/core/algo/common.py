#!/usr/bin/env python
#-*- coding:utf-8 -*-
from enum import Enum
from pyomo.environ import value

__name__ = ['AlgoIterReturnStatus', 'Algorithm']

class AlgoIterReturnStatus(Enum):
    OK = 1
    OPTIMIZATION_FAILED = 2
    OTHER = 0


class Algorithm():
    def __init__(self):
        self.model_history_data = {}
        self.plant_history_data = {}
        self.input_history_data = {}
        self.available_options = {}
        self.register_option()
        self.iter_count = -1
        self.spec_function=None
        self.plant_simulator=None
        self.model_simulator=None
        self.model_optimizer=None
        self.problem_description=None
        self.noise_generator=None
        self.perturbation_method=None


    def register_option(self):
        '''
        an example:
        self.available_options["filtering_factor"] = ("[0,1]",0.5) # accepted value, default value
        :return:
        '''
        self.available_options["filtering_factor"] = ("[0,1]", 0.5)
        self.available_options["homotopy_simulation"] = ("bool", True)
        self.available_options["homotopy_optimization"] = ("bool", True)
        self.available_options["noise_adding_fashion"] = ("integrated/added", "added")

    def set_algorithm_option(self, options=None):
        if options is None:
            options = {}
        self.options = {}
        for k in options.keys():
            if k not in self.available_options.keys():
                raise NameError("no such option %s"%k)
            self.options[k] = options[k]
        for k in self.available_options.keys():
            if k not in options.keys():
                self.options[k] = self.available_options[k][1]

    def set_problem(self, *args, **kwargs):
        return NotImplementedError

    def initialize_simulation(self, starting_point, initial_parameter_value):
        raise NotImplementedError()

    def one_step_simulation(self, *options):
        '''
        one step simulation gets the next trial point
        :param options:
        :return: AlgoIterReturnStatus
        '''
        raise NotImplementedError()


    def set_current_point(self, current_point):
        # set mv
        self.current_point = current_point
        for k, v in self.current_point.items():
            self.input_history_data[self.iter_count][k] = v

        # set specification
        if self.spec_function is not None:
            self.current_spec = self.spec_function.__call__(self.iter_count+1)
            for k, v in self.current_spec.items():
                self.input_history_data[self.iter_count][k] = v

    def set_initial_model_parameters(self, initial_parameter_value):
        # set parameter
        self.current_parameter_value = {}
        for p in self.model_simulator.pyomo_model.parameters.keys():
            if p in initial_parameter_value.keys():
                self.current_parameter_value[p] = initial_parameter_value[p]
            else:
                self.current_parameter_value[p] = self.model_simulator.pyomo_model.default_value[p]
        for k, v in self.current_parameter_value.items():
            self.model_history_data[0][k] = v

    def available_measurements(self):
        for k in self.problem_description.symbol_list['CV']:
            yield k
        yield self.problem_description.symbol_list['OBJ']
        for k in self.problem_description.symbol_list['CON']:
            yield k

    def filter_mv(self, optimized_input, mv_bounds):
        filtered_input = {}
        for k in self.problem_description.symbol_list['MV']:
            filtered_input[k] = self.current_point[k] + \
                                (optimized_input[k] - self.current_point[k]) * \
                                self.options["filtering_factor"]
            if mv_bounds[k][0] is not None and filtered_input[k] <= mv_bounds[k][0]:
                print("MV %s reaches its lower bound." % k)
                filtered_input[k] = mv_bounds[k][0]
            if mv_bounds[k][1] is not None and filtered_input[k] >= mv_bounds[k][1]:
                print("MV %s reaches its upper bound." % k)
                filtered_input[k] = mv_bounds[k][1]
        return filtered_input

    def adapt_to_bound_mv(self, optimized_input, mv_bounds):
        bounded_input = {}
        # for k in self.problem_description.symbol_list['MV']:
        for k in optimized_input.keys():
            bounded_input[k] = optimized_input[k]
            if mv_bounds[k][0] is not None and optimized_input[k] <= mv_bounds[k][0]:
                print("MV %s reaches its lower bound." % k)
                bounded_input[k] = mv_bounds[k][0]
            if mv_bounds[k][0] is not None and optimized_input[k] >= mv_bounds[k][1]:
                print("MV %s reaches its upper bound." % k)
                bounded_input[k] = mv_bounds[k][1]
        return bounded_input

    def optimize_for_u(self):
        input_values = {}
        if self.spec_function is not None:
            for k, v in self.current_spec.items():
                input_values[k] = v
        optimized_input, solve_status = self.model_optimizer.optimize(input_values,
                                                                      param_values=self.current_parameter_value,
                                                                      use_homo=self.options["homotopy_optimization"])
        # TODO:deal with solve status, fallback strategy
        return optimized_input, solve_status

    def get_trial_point(self):
        trial_points = self.perturbation_method.get_trial_points(self.current_point)
        if self.spec_function is not None:
            for i in range(len(trial_points)):
                for k,v in self.current_spec.items():
                    trial_points[i][k]=v
        return trial_points

    def get_plant_simulation_result(self, trial_points):
        plant_output_data = [{} for i in range(len(trial_points))]
        # self.plant_simulator.model.display(r"F:\Research\RTOdemo\debug\WO\result\display.txt")
        for i, p in enumerate(trial_points):
            if self.options["noise_adding_fashion"] == "integrated":
                for k,v in self.plant_simulator.pyomo_model.output_noise.items():
                    noise = self.noise_generator.get_noise(self.iter_count, i, k)
                    var=self.plant_simulator.pyomo_model.output_noise[k].__call__(self.plant_simulator.model)
                    var.fix(noise)
            outputs, solve_status = self.plant_simulator.simulate(p, param_values=None,
                                                                  use_homo=self.options["homotopy_simulation"])
            # TODO:deal with solve status
            if self.options["noise_adding_fashion"] == "integrated":
                noised_outputs = {}
                for op in self.plant_simulator.pyomo_model.noised_outputs.keys():
                    var = self.plant_simulator.pyomo_model.noised_outputs[op].__call__(self.plant_simulator.model)
                    noised_outputs[op] = value(var)

            # print(noised_outputs)
            # print(outputs)
            if i == 0:
                for k, v in outputs.items():
                    self.plant_history_data[self.iter_count][k] = v
            for k in self.available_measurements():
                if self.options["noise_adding_fashion"] == "integrated":
                    plant_output_data[i][k] = noised_outputs[k]
                elif self.options["noise_adding_fashion"] == "added":
                    plant_output_data[i][k] = outputs[k]

        # add noise to the plant data
        if self.options["noise_adding_fashion"] == "added":
            for trial_point_no in range(len(trial_points)):
                for k in plant_output_data[0].keys():
                    noise = self.noise_generator.get_noise(self.iter_count, trial_point_no, k)
                    plant_output_data[trial_point_no][k] += noise

        return plant_output_data

    def store_model_adaptation_data(self):
        for k, v in self.current_parameter_value.items():
            self.model_history_data[self.iter_count][k] = v

    def current_point_model_simulation(self, current_point):
        # get model simulation result with and without modifiers
        outputs, solve_status = self.model_simulator.simulate(current_point,
                                                              param_values=self.current_parameter_value,
                                                              use_homo=self.options["homotopy_simulation"])
        # TODO:deal with solve status
        for k, v in outputs.items():
            self.model_history_data[self.iter_count][k] = v

    def get_model_simulation_result(self, trial_points):
        model_output_data = [{} for i in range(len(trial_points))]
        for i, p in enumerate(trial_points):
            outputs, solve_status = self.model_simulator.simulate(p, param_values=self.current_parameter_value,
                                                                  use_homo=self.options["homotopy_simulation"])
            # TODO:deal with solve status
            if i == 0:
                for k, v in outputs.items():
                    self.model_history_data[self.iter_count][k] = v

            for k in self.available_measurements():
                if k + '_unmodified' in outputs.keys():
                    model_output_data[i][k] = outputs[k + '_unmodified']
        return model_output_data
            

class PE_type_Algorithm(Algorithm):
    def __init__(self):
        super().__init__()
        self.pe_estimator=None

    def register_option(self):
        super().register_option()
        self.available_options["pre-simulation_before_pe"] = ("bool", True)


    def set_initial_model_parameters(self, pyomo_model_parameters):
        super().set_initial_model_parameters(pyomo_model_parameters)

        self.pe_estimator.set_weight(self.output_weight, self.parameter_weight)
        self.pe_estimator.set_parameter_guess(self.parameter_initial_guess)
        # self.pe_estimator.model.display("pe_model.txt")
        # exit()

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
    def __init__(self):
        super().__init__()
        self.skipped_modifiers=[]

    def update_modifiers(self, plant_output_data, model_output_data):
        modifiers = self.perturbation_method.calculate_modifiers(plant_output_data, model_output_data)
        for k, v in modifiers.items():
            if k[1] is None:
                if k[0]+"_eps" not in self.skipped_modifiers:
                    self.current_parameter_value[k[0] + "_eps"] = v
            else:
                if k[1] + "_" + k[0] + "_lam" not in self.skipped_modifiers:
                    self.current_parameter_value[k[1] + "_" + k[0] + "_lam"] = v



