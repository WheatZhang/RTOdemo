#!/usr/bin/env python
#-*- coding:utf-8 -*-
from rtolib.model.wo_reactor_con import RTO_Plant_WO_reactor, RTO_Mismatched_WO_reactor,\
    RTO_Mismatched_WO_reactor_QC, default_WOR_description
from rtolib.core import NoiseGenerator, ModifierType,\
                    SimpleFiniteDiffPerturbation
from rtolib.core.algo import ModifierAdaptationTR, ModifierAdaptationPenaltyTR,\
    ModifierAdaptationCompoStepTR, ModifierAdaptation
import copy
from rtolib.util.misc import save_iteration_data_in_dict
import os
import pandas
import matplotlib.pyplot as plt
import numpy


global_parameter={
        "eta1": 0.01,
        "eta2": 0.9,
        "gamma1": 0.5,
        "gamma2": 1,
        "gamma3": 2,
        "MA_filtering": 0.2, #0.2
}

def generate_noise_file():
    noise_level = {
        'cost': 0,
        'con': 0,
        'XFr_A': 0,
        'XFr_B': 0,
        'XFr_E': 0,
        'XFr_P': 0,
        'XFr_G': 0,
    }
    max_iter = 200
    noise_generator = NoiseGenerator(noise_level)
    noise_generator.generate_noise(max_iter, 5)
    noise_generator.save_noise("noise/noise_0_wo.txt")

def algo3_TR_MA(perturbation_stepsize, starting_point, initial_trust_radius, max_trust_radius,\
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

    rto_algorithm = ModifierAdaptationTR()
    options = {
        "homotopy_simulation": False,
        "homotopy_optimization": False,
        "noise_adding_fashion": "integrated",
        "eta1": global_parameter["eta1"],
        "eta2": global_parameter["eta2"],
        "gamma1": global_parameter["gamma1"],
        "gamma2": global_parameter["gamma2"],
        "gamma3": global_parameter["gamma3"],
        "max_iter": max_iter,
        "max_trust_radius": max_trust_radius,
        "initial_trust_radius": initial_trust_radius,
    }
    rto_algorithm.set_algorithm_option(options)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=RTO_Plant_WO_reactor(),
        model=RTO_Mismatched_WO_reactor_QC(),
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        solver_executable=solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        )

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    while rto_algorithm.iter_count <= max_iter:
        rto_algorithm.one_step_simulation()
        if print_iter_data:
            print(rto_algorithm.model_history_data)
            print(rto_algorithm.plant_history_data)
            print(rto_algorithm.input_history_data)

    # save data
    save_iteration_data_in_dict(rto_algorithm.model_history_data, result_filename_header + "model_data.txt")
    save_iteration_data_in_dict(rto_algorithm.plant_history_data, result_filename_header + "plant_data.txt")
    save_iteration_data_in_dict(rto_algorithm.input_history_data, result_filename_header + "input_data.txt")


def algo2_TR_MA(perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,\
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

    rto_algorithm = ModifierAdaptationPenaltyTR()
    options = {
        "homotopy_simulation": False,
        "homotopy_optimization": False,
        "noise_adding_fashion": "integrated",
        "eta1": global_parameter["eta1"],
        "eta2": global_parameter["eta2"],
        "gamma1": global_parameter["gamma1"],
        "gamma2": global_parameter["gamma2"],
        "gamma3": global_parameter["gamma3"],
        "max_iter": max_iter,
        "sigma": sigma,
        "max_trust_radius":max_trust_radius,
        "initial_trust_radius":initial_trust_radius,
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

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    while rto_algorithm.iter_count <= max_iter:
        rto_algorithm.one_step_simulation()
        if print_iter_data:
            print(rto_algorithm.model_history_data)
            print(rto_algorithm.plant_history_data)
            print(rto_algorithm.input_history_data)

    # save data
    save_iteration_data_in_dict(rto_algorithm.model_history_data, result_filename_header + "model_data.txt")
    save_iteration_data_in_dict(rto_algorithm.plant_history_data, result_filename_header + "plant_data.txt")
    save_iteration_data_in_dict(rto_algorithm.input_history_data, result_filename_header + "input_data.txt")

def algo1_TR_MA(perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,xi_N,\
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

    rto_algorithm = ModifierAdaptationCompoStepTR()
    options = {
        "homotopy_simulation": False,
        "homotopy_optimization": False,
        "noise_adding_fashion": "integrated",
        "xi_N": xi_N,
        "eta1": global_parameter["eta1"],
        "eta2": global_parameter["eta2"],
        "gamma1": global_parameter["gamma1"],
        "gamma2": global_parameter["gamma2"],
        "gamma3": global_parameter["gamma3"],
        "max_iter": max_iter,
        "sigma": sigma,
        "max_trust_radius":max_trust_radius,
        "initial_trust_radius":initial_trust_radius,
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

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    while rto_algorithm.iter_count <= max_iter:
        rto_algorithm.one_step_simulation()
        if print_iter_data:
            print(rto_algorithm.model_history_data)
            print(rto_algorithm.plant_history_data)
            print(rto_algorithm.input_history_data)

    # save data
    save_iteration_data_in_dict(rto_algorithm.model_history_data, result_filename_header + "model_data.txt")
    save_iteration_data_in_dict(rto_algorithm.plant_history_data, result_filename_header + "plant_data.txt")
    save_iteration_data_in_dict(rto_algorithm.input_history_data, result_filename_header + "input_data.txt")


def original_MA(perturbation_stepsize, starting_point, \
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

    rto_algorithm = ModifierAdaptation()
    options = {
        "homotopy_simulation": False,
        "homotopy_optimization": False,
        "filtering_factor": global_parameter["MA_filtering"],
        "noise_adding_fashion": "integrated",
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

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    while rto_algorithm.iter_count <= max_iter:
        rto_algorithm.one_step_simulation()
        if print_iter_data:
            print(rto_algorithm.model_history_data)
            print(rto_algorithm.plant_history_data)
            print(rto_algorithm.input_history_data)

    # save data
    save_iteration_data_in_dict(rto_algorithm.model_history_data, result_filename_header + "model_data.txt")
    save_iteration_data_in_dict(rto_algorithm.plant_history_data, result_filename_header + "plant_data.txt")
    save_iteration_data_in_dict(rto_algorithm.input_history_data, result_filename_header + "input_data.txt")

def original_MA_with_QC(perturbation_stepsize, starting_point, \
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

    rto_algorithm = ModifierAdaptation()
    options = {
        "homotopy_simulation": False,
        "homotopy_optimization": False,
        "filtering_factor": global_parameter["MA_filtering"],
        "noise_adding_fashion": "integrated",
    }
    rto_algorithm.set_algorithm_option(options)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=RTO_Plant_WO_reactor(),
        model=RTO_Mismatched_WO_reactor_QC(),
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        solver_executable=solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        )

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    while rto_algorithm.iter_count <= max_iter:
        rto_algorithm.one_step_simulation()
        if print_iter_data:
            print(rto_algorithm.model_history_data)
            print(rto_algorithm.plant_history_data)
            print(rto_algorithm.input_history_data)

    # save data
    save_iteration_data_in_dict(rto_algorithm.model_history_data, result_filename_header + "model_data.txt")
    save_iteration_data_in_dict(rto_algorithm.plant_history_data, result_filename_header + "plant_data.txt")
    save_iteration_data_in_dict(rto_algorithm.input_history_data, result_filename_header + "input_data.txt")


def do_test():
    # ------------------------------------
    noise_filename = "noise/noise_0_wo.txt"
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    result_filename_folder="data/wo_example/"
    if not os.path.exists(result_filename_folder):
        os.makedirs(result_filename_folder)

    # ------------------------------------
    starting_point = {
        "Fa": 3.6,
        "Fb": 10,
        "Tr": 85,
    }
    # starting_point = {
    #     "Fa": 3.8,
    #     "Fb": 9.3,
    #     "Tr": 91,
    # }

    # ------------------------------------
    perturbation_stepsize = {
        "Fa": 0.0001,
        "Fb": 0.0001,
        "Tr": 0.0001,
    }
    # ------------------------------------
    print_iter_data = False
    max_iter = 80
    initial_trust_radius = 0.1
    sigma = 100
    xi_N = 0.5
    max_trust_radius=1

    # # # ------------------------------------
    print("\nTesting TR_MA")
    result_filename_header = result_filename_folder + "TR_MA_"
    algo3_TR_MA(perturbation_stepsize, starting_point, initial_trust_radius,max_trust_radius,\
               noise_filename, solver_executable, print_iter_data, max_iter, \
               result_filename_header)

    # # ------------------------------------
    # print("\nTesting Penalty_TR_MA")
    # result_filename_header = result_filename_folder + "Penalty_TR_MA_"
    # algo2_TR_MA(perturbation_stepsize, starting_point, sigma, initial_trust_radius,max_trust_radius,\
    #             noise_filename, solver_executable, print_iter_data, max_iter, \
    #             result_filename_header)
    #
    # # ------------------------------------
    # print("\nTesting CompoStep_TR_MA")
    # result_filename_header = result_filename_folder + "CompoStep_TR_MA_"
    # algo1_TR_MA(perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,\
    #             xi_N,\
    #             noise_filename, solver_executable, print_iter_data, max_iter, \
    #             result_filename_header)

    # ------------------------------------
    print("\nTesting Original_MA")
    result_filename_header = result_filename_folder + "Original_MA_"
    original_MA(perturbation_stepsize, starting_point, \
                noise_filename, solver_executable, print_iter_data, max_iter, \
                result_filename_header)

def batch_test(batch_name, initial_trust_radius, sigma, test_algo3=True):
    # ------------------------------------
    noise_filename = "noise/noise_0_wo.txt"
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    result_filename_folder="data/wo_example/batch/"
    if not os.path.exists(result_filename_folder):
        os.makedirs(result_filename_folder)

    # ------------------------------------
    starting_point = {
        "Fa": 3.6,
        "Fb": 10,
        "Tr": 85,
    }

    # ------------------------------------
    perturbation_stepsize = {
        "Fa": 0.0001,
        "Fb": 0.0001,
        "Tr": 0.0001,
    }
    # ------------------------------------
    print_iter_data = False
    max_iter = 80
    max_trust_radius=1
    xi_N=0.5

    # # ------------------------------------
    if test_algo3:
        print("\nTesting TR_MA")
        result_filename_header = result_filename_folder + batch_name+"_"+"TR_MA_"
        algo3_TR_MA(perturbation_stepsize, starting_point, initial_trust_radius,max_trust_radius,\
                   noise_filename, solver_executable, print_iter_data, max_iter, \
                   result_filename_header)

    # ------------------------------------
    print("\nTesting Penalty_TR_MA")
    result_filename_header = result_filename_folder + batch_name+"_"+ "Penalty_TR_MA_"
    algo2_TR_MA(perturbation_stepsize, starting_point, sigma, initial_trust_radius,max_trust_radius,\
                noise_filename, solver_executable, print_iter_data, max_iter, \
                result_filename_header)

    # ------------------------------------
    print("\nTesting CompoStep_TR_MA")
    result_filename_header = result_filename_folder + batch_name+"_"+ "CompoStep_TR_MA_"
    algo1_TR_MA(perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,\
                xi_N,\
                noise_filename, solver_executable, print_iter_data, max_iter, \
                result_filename_header)


def plot_profile():
    max_iter=80
    global_marker_size = 2
    linewidth=1
    compo_step_plant_data = pandas.read_csv("data/wo_example/CompoStep_TR_MA_plant_data.txt", \
                                      index_col=0, header=0, sep='\t')
    penalty_plant_data = pandas.read_csv("data/wo_example/Penalty_TR_MA_plant_data.txt", \
                                   index_col=0, header=0, sep='\t')
    inf_averse_plant_data = pandas.read_csv("data/wo_example/TR_MA_plant_data.txt", \
                                         index_col=0, header=0, sep='\t')
    original_ma_plant_data = pandas.read_csv("data/wo_example/Original_MA_plant_data.txt", \
                                            index_col=0, header=0, sep='\t')

    compo_step_model_data = pandas.read_csv("data/wo_example/CompoStep_TR_MA_model_data.txt", \
                                      index_col=0, header=0, sep='\t')
    penalty_model_data = pandas.read_csv("data/wo_example/Penalty_TR_MA_model_data.txt", \
                                   index_col=0, header=0, sep='\t')
    inf_averse_model_data = pandas.read_csv("data/wo_example/TR_MA_model_data.txt", \
                                            index_col=0, header=0, sep='\t')
    original_ma_model_data = pandas.read_csv("data/wo_example/Original_MA_model_data.txt", \
                                             index_col=0, header=0, sep='\t')

    compo_step_input_data = pandas.read_csv("data/wo_example/CompoStep_TR_MA_input_data.txt", \
                                            index_col=0, header=0, sep='\t')
    penalty_input_data = pandas.read_csv("data/wo_example/Penalty_TR_MA_input_data.txt", \
                                         index_col=0, header=0, sep='\t')
    inf_averse_input_data = pandas.read_csv("data/wo_example/TR_MA_input_data.txt", \
                                            index_col=0, header=0, sep='\t')
    original_ma_input_data = pandas.read_csv("data/wo_example/Original_MA_input_data.txt", \
                                             index_col=0, header=0, sep='\t')

    fig = plt.figure(figsize=(6,12))
    plt.subplot(711)
    optimal = -210.33 * numpy.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1,max_iter+1), compo_step_plant_data.loc[1:max_iter, 'cost'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), penalty_plant_data.loc[1:max_iter, 'cost'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), inf_averse_plant_data.loc[1:max_iter, 'cost'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), original_ma_plant_data.loc[1:max_iter, 'cost'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel("plant cost")
    plt.subplot(712)
    optimal = 0 * numpy.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_plant_data.loc[1:max_iter, 'con'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), penalty_plant_data.loc[1:max_iter, 'con'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), inf_averse_plant_data.loc[1:max_iter, 'con'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), original_ma_plant_data.loc[1:max_iter, 'con'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel("plant constraints")
    plt.subplot(713)
    optimal = 3.887 * numpy.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_input_data.loc[1:max_iter, 'Fa'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), penalty_input_data.loc[1:max_iter, 'Fa'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), inf_averse_input_data.loc[1:max_iter, 'Fa'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), original_ma_input_data.loc[1:max_iter, 'Fa'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"Fa")
    plt.subplot(714)
    optimal = 9.369 * numpy.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_input_data.loc[1:max_iter, 'Fb'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), penalty_input_data.loc[1:max_iter, 'Fb'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), inf_averse_input_data.loc[1:max_iter, 'Fb'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), original_ma_input_data.loc[1:max_iter, 'Fb'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"Fb")
    plt.subplot(715)
    optimal = 91.2 * numpy.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_input_data.loc[1:max_iter, 'Tr'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), penalty_input_data.loc[1:max_iter, 'Tr'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), inf_averse_input_data.loc[1:max_iter, 'Tr'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), original_ma_input_data.loc[1:max_iter, 'Tr'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"Tr")
    plt.subplot(716)
    plt.plot(range(1, max_iter + 1), compo_step_model_data.loc[1:max_iter, 'rho'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), penalty_model_data.loc[1:max_iter, 'rho'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), inf_averse_model_data.loc[1:max_iter, 'rho'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"$\rho$")
    plt.subplot(717)
    plt.plot(range(1, max_iter + 1), compo_step_model_data.loc[1:max_iter, 'tr'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), penalty_model_data.loc[1:max_iter, 'tr'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), inf_averse_model_data.loc[1:max_iter, 'tr'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel("radius")
    plt.xlabel("#iteration")

    plt.savefig("pic/wo-example/profile", dpi=600)
    plt.close()

def batch_plot_profile(batch_i, batch_j):
    max_iter=80
    global_marker_size = 2
    linewidth=1
    batch_prefix = "%d_%d_"%(batch_i, batch_j)
    batch_prefix_for_algo3 = "%d_0_"%(batch_i)
    compo_step_plant_data = pandas.read_csv("data/wo_example/batch/"+batch_prefix+"CompoStep_TR_MA_plant_data.txt", \
                                      index_col=0, header=0, sep='\t')
    penalty_plant_data = pandas.read_csv("data/wo_example/batch/"+batch_prefix+"Penalty_TR_MA_plant_data.txt", \
                                   index_col=0, header=0, sep='\t')
    inf_averse_plant_data = pandas.read_csv("data/wo_example/batch/"+batch_prefix_for_algo3+"TR_MA_plant_data.txt", \
                                         index_col=0, header=0, sep='\t')
    original_ma_plant_data = pandas.read_csv("data/wo_example/Original_MA_plant_data.txt", \
                                            index_col=0, header=0, sep='\t')

    compo_step_model_data = pandas.read_csv("data/wo_example/batch/"+batch_prefix+"CompoStep_TR_MA_model_data.txt", \
                                      index_col=0, header=0, sep='\t')
    penalty_model_data = pandas.read_csv("data/wo_example/batch/"+batch_prefix+"Penalty_TR_MA_model_data.txt", \
                                   index_col=0, header=0, sep='\t')
    inf_averse_model_data = pandas.read_csv("data/wo_example/batch/"+batch_prefix_for_algo3+"TR_MA_model_data.txt", \
                                            index_col=0, header=0, sep='\t')
    original_ma_model_data = pandas.read_csv("data/wo_example/Original_MA_model_data.txt", \
                                             index_col=0, header=0, sep='\t')

    compo_step_input_data = pandas.read_csv("data/wo_example/batch/"+batch_prefix+"CompoStep_TR_MA_input_data.txt", \
                                            index_col=0, header=0, sep='\t')
    penalty_input_data = pandas.read_csv("data/wo_example/batch/"+batch_prefix+"Penalty_TR_MA_input_data.txt", \
                                         index_col=0, header=0, sep='\t')
    inf_averse_input_data = pandas.read_csv("data/wo_example/batch/"+batch_prefix_for_algo3+"TR_MA_input_data.txt", \
                                            index_col=0, header=0, sep='\t')
    original_ma_input_data = pandas.read_csv("data/wo_example/Original_MA_input_data.txt", \
                                             index_col=0, header=0, sep='\t')

    fig = plt.figure(figsize=(6,12))
    plt.subplot(711)
    optimal = -210.33 * numpy.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1,max_iter+1), compo_step_plant_data.loc[1:max_iter, 'cost'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), penalty_plant_data.loc[1:max_iter, 'cost'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), inf_averse_plant_data.loc[1:max_iter, 'cost'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), original_ma_plant_data.loc[1:max_iter, 'cost'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel("plant cost")
    plt.subplot(712)
    optimal = 0 * numpy.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_plant_data.loc[1:max_iter, 'con'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), penalty_plant_data.loc[1:max_iter, 'con'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), inf_averse_plant_data.loc[1:max_iter, 'con'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), original_ma_plant_data.loc[1:max_iter, 'con'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel("plant constraints")
    plt.subplot(713)
    optimal = 3.887 * numpy.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_input_data.loc[1:max_iter, 'Fa'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), penalty_input_data.loc[1:max_iter, 'Fa'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), inf_averse_input_data.loc[1:max_iter, 'Fa'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), original_ma_input_data.loc[1:max_iter, 'Fa'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"Fa")
    plt.subplot(714)
    optimal = 9.369 * numpy.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_input_data.loc[1:max_iter, 'Fb'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), penalty_input_data.loc[1:max_iter, 'Fb'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), inf_averse_input_data.loc[1:max_iter, 'Fb'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), original_ma_input_data.loc[1:max_iter, 'Fb'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"Fb")
    plt.subplot(715)
    optimal = 91.2 * numpy.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_input_data.loc[1:max_iter, 'Tr'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), penalty_input_data.loc[1:max_iter, 'Tr'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), inf_averse_input_data.loc[1:max_iter, 'Tr'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), original_ma_input_data.loc[1:max_iter, 'Tr'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"Tr")
    plt.subplot(716)
    plt.plot(range(1, max_iter + 1), compo_step_model_data.loc[1:max_iter, 'rho'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), penalty_model_data.loc[1:max_iter, 'rho'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), inf_averse_model_data.loc[1:max_iter, 'rho'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"$\rho$")
    plt.subplot(717)
    plt.plot(range(1, max_iter + 1), compo_step_model_data.loc[1:max_iter, 'tr'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), penalty_model_data.loc[1:max_iter, 'tr'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), inf_averse_model_data.loc[1:max_iter, 'tr'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel("radius")
    plt.xlabel("#iteration")

    plt.savefig("pic/wo-example/batch/"+batch_prefix+"profile", dpi=600)
    plt.close()

def do_batch_test():
    initial_trust_radius_options = [0.001, 0.005 ,0.01, 0.05, 0.1, 0.5]
    sigma_options = [0, 1, 2, 5, 1e1, 2e1, 5e1, 1e2, 2e2, 5e2, 1e3, 2e3, 5e3]
    for i, initial_trust_radius in enumerate(initial_trust_radius_options):
        for j, sigma in enumerate(sigma_options):
            if j == 0:
                batch_test("%d_%d"%(i,j), initial_trust_radius=initial_trust_radius, sigma=sigma,\
                           test_algo3=True)
            else:
                batch_test("%d_%d"%(i,j), initial_trust_radius=initial_trust_radius, sigma=sigma,\
                           test_algo3=False)
            batch_plot_profile(i, j)


def rearrange_pic():
    import shutil

    initial_trust_radius_options = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5]
    sigma_options = [0, 1, 2, 5, 1e1, 2e1, 5e1, 1e2, 2e2, 5e2, 1e3, 2e3, 5e3]
    for i, initial_trust_radius in enumerate(initial_trust_radius_options):
        for j, sigma in enumerate(sigma_options):
            batch_prefix = "%d_%d_" % (i, j)
            shutil.copyfile("pic/wo-example/batch/"+batch_prefix+"profile.png", \
                            "pic/wo-example/selected_pic/"+"%dS_"%j+batch_prefix+"profile.png")


if __name__ == "__main__":
    # generate_noise_file()
    do_test()
    plot_profile()

    # do_batch_test()

    # rearrange_pic()