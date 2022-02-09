#!/usr/bin/env python
#-*- coding:utf-8 -*-
from rtolib.model.wo_reactor import RTO_Plant_WO_reactor, RTO_Mismatched_WO_reactor,\
    default_WOR_description
from rtolib.core import NoiseGenerator, ModifierType,\
                    SimpleFiniteDiffPerturbation
from rtolib.core.algo import ModifierAdaptation, GeneralizedParameterEstimation, ISOPE_Algorithm,\
        ITSParameterEstimation
import copy
from rtolib.util.misc import save_iteration_data_in_dict


def generate_noise_file(profit_noise_level, composition_noise_level, noise_filename):
    noise_level = {
            'profit': profit_noise_level,
            'XFr_A': composition_noise_level,
            'XFr_B': composition_noise_level,
            'XFr_E': composition_noise_level,
            'XFr_P': composition_noise_level,
            'XFr_G': composition_noise_level,
    }
    max_iter = 200
    noise_generator = NoiseGenerator(noise_level)
    for i in range(max_iter):
        for j in range(len(noise_level)):
            for k in noise_level.keys():
                noise_generator.get_noise(i, j, k)
    noise_generator.save_noise(noise_filename)


def do_test_MA(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, solver_executable, print_iter_data, max_iter,\
               result_filename_header):
    '''

    :param perturbation_stepsize: dict
    :param starting_point: dict
    :param filtering_factor: double within [0,1]
    :param noise_filename: string
    :param solver_executable: string
    :param print_iter_data: bool
    :param max_iter: int, e.g. 20
    :param result_filename_header: string
    :return:
    '''
    problem_description = copy.deepcopy(default_WOR_description)
    problem_description.symbol_list['CV']=('profit',)

    rto_algorithm = ModifierAdaptation()
    options = {
        "homotopy_simulation": True,
        "homotopy_optimization": True,
        "filtering_factor": filtering_factor,
    }
    rto_algorithm.set_algorithm_option(options)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=RTO_Plant_WO_reactor(),
        model=RTO_Mismatched_WO_reactor(),
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        solver_executable=solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        )

    homotopy_var=['Fb','Tr']
    rto_algorithm.set_homotopy_var(homotopy_var)

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    for i in range(max_iter):
        rto_algorithm.one_step_simulation()
        if print_iter_data:
            print(rto_algorithm.model_history_data)
            print(rto_algorithm.plant_history_data)
            print(rto_algorithm.input_history_data)

    # save data
    save_iteration_data_in_dict(rto_algorithm.model_history_data, result_filename_header + "model_data.txt")
    save_iteration_data_in_dict(rto_algorithm.plant_history_data, result_filename_header + "plant_data.txt")
    save_iteration_data_in_dict(rto_algorithm.input_history_data, result_filename_header + "input_data.txt")


def do_test_MAy(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, solver_executable, print_iter_data, max_iter,\
               result_filename_header):
    problem_description = copy.deepcopy(default_WOR_description)

    rto_algorithm = ModifierAdaptation()
    options = {
        "homotopy_simulation": True,
        "homotopy_optimization": True,
        "filtering_factor": filtering_factor,
    }
    rto_algorithm.set_algorithm_option(options)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

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

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    for i in range(max_iter):
        rto_algorithm.one_step_simulation()
        if print_iter_data:
            print(rto_algorithm.model_history_data)
            print(rto_algorithm.plant_history_data)
            print(rto_algorithm.input_history_data)

    # save data
    save_iteration_data_in_dict(rto_algorithm.model_history_data, result_filename_header + "model_data.txt")
    save_iteration_data_in_dict(rto_algorithm.plant_history_data, result_filename_header + "plant_data.txt")
    save_iteration_data_in_dict(rto_algorithm.input_history_data, result_filename_header + "input_data.txt")


def do_test_PE(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, solver_executable, print_iter_data, max_iter,\
               result_filename_header, output_noise_level):
    problem_description = copy.deepcopy(default_WOR_description)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    rto_algorithm = ITSParameterEstimation()
    options={
        "homotopy_optimization": True,
        "pre-simulation_before_pe": True,
        "filtering_factor": filtering_factor,
    }
    rto_algorithm.set_algorithm_option(options)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    parameter_weight = {
        "Ka1": 0,
        "Ka2": 0,
        "Kb1": 0,
        "Kb2": 0,
    }

    output_weight = {
        "XFr_A": 1/output_noise_level/output_noise_level,
        "XFr_B": 1/output_noise_level/output_noise_level,
        "XFr_E": 1/output_noise_level/output_noise_level,
        "XFr_P": 1/output_noise_level/output_noise_level,
        "XFr_G": 1/output_noise_level/output_noise_level,
    }

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=RTO_Plant_WO_reactor(),
        model=RTO_Mismatched_WO_reactor(),
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        solver_executable=solver_executable,
        spec_function=None,
        parameter_weight=parameter_weight,
        output_weight=output_weight,
        )

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    for i in range(max_iter):
        rto_algorithm.one_step_simulation()
        if print_iter_data:
            print(rto_algorithm.model_history_data)
            print(rto_algorithm.plant_history_data)
            print(rto_algorithm.input_history_data)

    # save data
    save_iteration_data_in_dict(rto_algorithm.model_history_data, result_filename_header + "model_data.txt")
    save_iteration_data_in_dict(rto_algorithm.plant_history_data, result_filename_header + "plant_data.txt")
    save_iteration_data_in_dict(rto_algorithm.input_history_data, result_filename_header + "input_data.txt")


def do_test_GPE(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, solver_executable, print_iter_data, max_iter,\
               result_filename_header,profit_noise_level,output_noise_level,\
                ka_relative_uncertainty,kb_relative_uncertainty,factor_n):
    problem_description = copy.deepcopy(default_WOR_description)
    problem_description.symbol_list['CV']=[]

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    rto_algorithm = GeneralizedParameterEstimation()
    options={
        "homotopy_optimization": True,
        "pre-simulation_before_pe": False,
        "filtering_factor": filtering_factor,
    }
    rto_algorithm.set_algorithm_option(options)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    

    parameter_weight = {
        "Ka1": 1 / ka_relative_uncertainty / ka_relative_uncertainty/2.189e8/2.189e8,
        "Ka2": 1 / ka_relative_uncertainty / ka_relative_uncertainty/4.310e13/4.310e13,
        "Kb1": 1 / kb_relative_uncertainty / kb_relative_uncertainty/8077.6/8077.6,
        "Kb2": 1 / kb_relative_uncertainty / kb_relative_uncertainty/12438/12438,
        "profit_eps":factor_n * (1 / profit_noise_level / profit_noise_level),
        "Fb_profit_lam":factor_n * (1 / profit_noise_level / profit_noise_level / 2 * (perturbation_stepsize['Fb'] ** 2)),
        "Tr_profit_lam":factor_n * (1 / profit_noise_level / profit_noise_level / 2 * (perturbation_stepsize['Tr'] ** 2)),
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

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    for i in range(max_iter):
        rto_algorithm.one_step_simulation()
        if print_iter_data:
            print(rto_algorithm.model_history_data)
            print(rto_algorithm.plant_history_data)
            print(rto_algorithm.input_history_data)

    # save data
    save_iteration_data_in_dict(rto_algorithm.model_history_data, result_filename_header + "model_data.txt")
    save_iteration_data_in_dict(rto_algorithm.plant_history_data, result_filename_header + "plant_data.txt")
    save_iteration_data_in_dict(rto_algorithm.input_history_data, result_filename_header + "input_data.txt")

    
def do_test_ISOPE(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, solver_executable, print_iter_data, max_iter,\
               result_filename_header, output_noise_level):
    problem_description = copy.deepcopy(default_WOR_description)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    rto_algorithm = ISOPE_Algorithm()
    options={
        "homotopy_optimization": True,
        "pre-simulation_before_pe": True,
        "filtering_factor": filtering_factor,
    }
    rto_algorithm.set_algorithm_option(options)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    parameter_weight = {
        "Ka1": 0,
        "Ka2": 0,
        "Kb1": 0,
        "Kb2": 0,
    }

    output_weight = {
        "XFr_A": 1/output_noise_level/output_noise_level,
        "XFr_B": 1/output_noise_level/output_noise_level,
        "XFr_E": 1/output_noise_level/output_noise_level,
        "XFr_P": 1/output_noise_level/output_noise_level,
        "XFr_G": 1/output_noise_level/output_noise_level,
    }

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=RTO_Plant_WO_reactor(),
        model=RTO_Mismatched_WO_reactor(),
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        solver_executable=solver_executable,
        spec_function=None,
        parameter_weight=parameter_weight,
        output_weight=output_weight,
        fixed_parameter_values=None,
        )

    homotopy_var = ['Fb', 'Tr']
    rto_algorithm.set_homotopy_var(homotopy_var)

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    for i in range(max_iter):
        rto_algorithm.one_step_simulation()
        if print_iter_data:
            print(rto_algorithm.model_history_data)
            print(rto_algorithm.plant_history_data)
            print(rto_algorithm.input_history_data)

    # save data
    save_iteration_data_in_dict(rto_algorithm.model_history_data, result_filename_header + "model_data.txt")
    save_iteration_data_in_dict(rto_algorithm.plant_history_data, result_filename_header + "plant_data.txt")
    save_iteration_data_in_dict(rto_algorithm.input_history_data, result_filename_header + "input_data.txt")
