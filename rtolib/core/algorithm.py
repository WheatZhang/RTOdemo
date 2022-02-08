from enum import Enum

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

    def register_option(self):
        '''
        an example:
        self.available_options["filtering_factor"] = ("[0,1]",0.5) # accepted value, default value
        :return:
        '''
        raise NotImplementedError()

    def set_algorithm_option(self, options={}):
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

    def initialize_simulation(self):
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
            outputs, solve_status = self.plant_simulator.simulate(p, param_values=None,
                                                                  use_homo=self.options["homotopy_simulation"])
            # TODO:deal with solve status
            if i == 0:
                for k, v in outputs.items():
                    self.plant_history_data[self.iter_count][k] = v
            for k in self.problem_description.symbol_list['CV']:
                plant_output_data[i][k] = outputs[k]

        # add noise to the plant data
        for trial_point_no in range(len(trial_points)):
            for k in self.problem_description.symbol_list['CV']:
                noise = self.noise_generator.get_noise(self.iter_count, trial_point_no, k)
                plant_output_data[trial_point_no][k] += noise

        return plant_output_data

    def store_model_adaptation_data(self):
        for k, v in self.current_parameter_value.items():
            self.model_history_data[self.iter_count][k] = v

    def current_point_model_simulation(self, current_point):
        # get model simulation result with and without modifiers
        self.model_simulator.model.display("result/display.txt")
        outputs, solve_status = self.model_simulator.simulate(current_point,
                                                              param_values=self.current_parameter_value,
                                                              use_homo=self.options["homotopy_simulation"])
        # TODO:deal with solve status
        for k, v in outputs.items():
            self.model_history_data[self.iter_count][k] = v