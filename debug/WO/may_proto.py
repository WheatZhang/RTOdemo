#!/usr/bin/env python
#-*- coding:utf-8 -*-
import copy
from rtolib.core import NoiseGenerator, \
    SimpleFiniteDiffPerturbation,ModifierType
from rtolib.core.algo.ma import ModifierAdaptation
from rtolib.model.wo_reactor import RTO_Plant_WO_reactor, RTO_Mismatched_WO_reactor, \
    default_WOR_description
from rtolib.util.misc import save_iteration_data_in_dict

def generate_noise_file():
    noise_level = {
            'profit': 0.1,
            'XFr_A': 1e-4,
            'XFr_B': 1e-4,
            'XFr_E': 1e-4,
            'XFr_P': 1e-4,
            'XFr_G': 1e-4,
    }
    max_iter = 200
    noise_generator = NoiseGenerator(noise_level)
    for i in range(max_iter):
        for j in range(len(noise_level)):
            for k in noise_level.keys():
                noise_generator.get_noise(i, j, k)
    noise_generator.save_noise("noise/noise1.txt")


def generate_zero_noise_file():
    noise_level = {
            'profit': 0,
            'XFr_A': 0,
            'XFr_B': 0,
            'XFr_E': 0,
            'XFr_P': 0,
            'XFr_G': 0,
    }
    max_iter = 200
    noise_generator = NoiseGenerator(noise_level)
    for i in range(max_iter):
        for j in range(len(noise_level)):
            for k in noise_level.keys():
                noise_generator.get_noise(i, j, k)
    noise_generator.save_noise("noise/uncon-wo-noise0.txt")


def do_test():
    problem_description = copy.deepcopy(default_WOR_description)

    rto_algorithm = ModifierAdaptation()
    options = {
        "homotopy_simulation": True,
        "homotopy_optimization": True,
        "filtering_factor": 0.5,
    }
    rto_algorithm.set_algorithm_option(options)

    pertubation_stepsize = {
        "Fb": 0.2,
        "Tr": 2,
    }
    ffd_perturb = SimpleFiniteDiffPerturbation(pertubation_stepsize, problem_description)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise("noise/uncon-wo-noise0.txt")

    solver_executable=r"../../external/bin/ipopt.exe"

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=RTO_Plant_WO_reactor(),
        model=RTO_Mismatched_WO_reactor(),
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        solver_executable=solver_executable,
        spec_function=None,
        modifier_type=ModifierType.OUTPUT,
        )

    homotopy_var=['Fb','Tr']
    rto_algorithm.set_homotopy_var(homotopy_var)

    starting_point={
        "Fb": 4,
        "Tr": 75,
    }
    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    max_iter = 20
    for i in range(max_iter):
        rto_algorithm.one_step_simulation()
        print(rto_algorithm.model_history_data)
        print(rto_algorithm.plant_history_data)
        print(rto_algorithm.input_history_data)

    # save data
    filename_header = "result/"
    save_iteration_data_in_dict(rto_algorithm.model_history_data, filename_header + "model_data.txt")
    save_iteration_data_in_dict(rto_algorithm.plant_history_data, filename_header + "plant_data.txt")
    save_iteration_data_in_dict(rto_algorithm.input_history_data, filename_header + "input_data.txt")



if __name__ == "__main__":
    # generate_zero_noise_file()
    # generate_noise_file()
    do_test()