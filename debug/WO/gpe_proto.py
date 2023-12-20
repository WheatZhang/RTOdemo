#!/usr/bin/env python
#-*- coding:utf-8 -*-
import copy
from rtolib.core import NoiseGenerator, \
    SimpleFiniteDiffPerturbation,ModifierType
from rtolib.core.algo.gpe import GeneralizedParameterEstimation
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
    problem_description.symbol_list['CV']=[]

    pertubation_stepsize = {
        "Fb": 0.2,
        "Tr": 2,
    }
    ffd_perturb = SimpleFiniteDiffPerturbation(pertubation_stepsize, problem_description)

    rto_algorithm = GeneralizedParameterEstimation()
    options={
        "homotopy_optimization": True,
        "pre-simulation_before_pe": False,
        "filtering_factor": 0.5,
    }
    rto_algorithm.set_algorithm_option(options)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise("noise/uncon-wo-noise0.txt")

    solver_executable=r"F:\Research\RTOdemo\external\bin\ipopt.exe"

    profit_noise_level = 0.6
    output_noise_level = 1e-4
    ka_uncertainty = 0.01  # after scaling
    kb_uncertainty = 0.1  # after scaling
    n = 10

    parameter_weight = {
        "Ka1": 1 / ka_uncertainty / ka_uncertainty/2.189e8/2.189e8,
        "Ka2": 1 / ka_uncertainty / ka_uncertainty/4.310e13/4.310e13,
        "Kb1": 1 / kb_uncertainty / kb_uncertainty/8077.6/8077.6,
        "Kb2": 1 / kb_uncertainty / kb_uncertainty/12438/12438,
        "profit_eps":n * (1 / profit_noise_level / profit_noise_level),
        "Fb_profit_lam":n * (1 / profit_noise_level / profit_noise_level / 2 * (pertubation_stepsize['Fb'] ** 2)),
        "Tr_profit_lam":n * (1 / profit_noise_level / profit_noise_level / 2 * (pertubation_stepsize['Tr'] ** 2)),
    }

    output_weight = {
        'profit': 1 / profit_noise_level / profit_noise_level,
    }

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=RTO_Plant_WO_reactor(),
        model=RTO_Mismatched_WO_reactor(),
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        solver_executable=solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        parameter_weight=parameter_weight,
        output_weight=output_weight,
        )

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