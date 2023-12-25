from rtolib.model.parallel_wo_uncon import RTO_Plant_Parallel_WO, RTO_Mismatched_BlackBox_Parallel_WO,\
    RTO_Mismatched_Parallel_WO, default_WOR_description
from rtolib.core import NoiseGenerator, ModifierType,\
                    SimpleFiniteDiffPerturbation
from rtolib.core.algo import ModifierAdaptationCompoStepTR
import copy
from rtolib.util.misc import save_iteration_data_in_dict
import os
import pandas
import matplotlib.pyplot as plt
import numpy as np


global_parameter={
        "eta1": 0.01,
        "eta2": 0.9,
        "gamma1": 0.5,
        "gamma2": 1,
        "gamma3": 2,
        "feasibility_tol": 1e-3,
        "stationarity_tol": 1e-3,
}

def generate_noise_file():
    noise_level = {
        'cost': 0,
    }
    max_iter = 200
    noise_generator = NoiseGenerator(noise_level)
    noise_generator.generate_noise(max_iter, 10)
    noise_generator.save_noise("noise/noise_0_try_pso_uncon.txt")

def compo_step_TR_pyomo_model(perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,xi_N,\
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
        "noise_adding_fashion": "added",
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
        'feasibility_tol': global_parameter["feasibility_tol"],
        "stationarity_tol": global_parameter["stationarity_tol"],
        'adaptive_sigma': True,
    }
    rto_algorithm.set_algorithm_option(options)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=RTO_Plant_Parallel_WO(),
        model=RTO_Mismatched_Parallel_WO(),
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        solver_executable=solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        black_box_model=None,
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


def compo_step_TR_blackbox_model(perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,xi_N,\
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
        "noise_adding_fashion": "added",
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
        'feasibility_tol': global_parameter["feasibility_tol"],
        "stationarity_tol": global_parameter["stationarity_tol"],
        'adaptive_sigma': True,
    }
    rto_algorithm.set_algorithm_option(options)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=RTO_Plant_Parallel_WO(),
        model=RTO_Mismatched_Parallel_WO(),
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        solver_executable=solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        black_box_model=RTO_Mismatched_BlackBox_Parallel_WO(),
        )

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})
    pso_options = {
        'population': 200,
        'max_iter':10,
    }
    rto_algorithm.model_optimizer.set_solver_options('PSO', pso_options)

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

def compo_step_TR_blackbox_model2(perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,xi_N,\
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
        "noise_adding_fashion": "added",
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
        'feasibility_tol': global_parameter["feasibility_tol"],
        "stationarity_tol": global_parameter["stationarity_tol"],
        'adaptive_sigma': True,
    }
    rto_algorithm.set_algorithm_option(options)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=RTO_Plant_Parallel_WO(),
        model=RTO_Mismatched_Parallel_WO(),
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        solver_executable=solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        black_box_model=RTO_Mismatched_BlackBox_Parallel_WO(),
        )

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})
    pso_options = {
        'population': 10,
        'max_iter':5,
    }
    rto_algorithm.model_optimizer.set_solver_options('PSO', pso_options)

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
    noise_filename = "noise/noise_0_try_pso_uncon.txt"
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    result_filename_folder="data/try_pso_uncon/"
    if not os.path.exists(result_filename_folder):
        os.makedirs(result_filename_folder)

    # ------------------------------------
    starting_point = {
        "Fb1": 4,
        "Fb2": 4,
        "Tr1": 80,
        "Tr2": 80,
        "Tr3": 80,
    }

    # ------------------------------------
    perturbation_stepsize = {
        "Fb1": 0.001,
        "Fb2": 0.001,
        "Tr1": 0.01,
        "Tr2": 0.01,
        "Tr3": 0.01,
    }
    # ------------------------------------
    print_iter_data = False
    max_iter = 20
    initial_trust_radius = 0.2
    sigma = 100
    xi_N = 0.5
    kappa_b=0.2
    max_trust_radius=1

    # ------------------------------------
    # print("\nTesting CompoStep_TR_MA with pyomo model")
    # result_filename_header = result_filename_folder + "CompoStep_TR_Pyomo_"
    # compo_step_TR_pyomo_model(perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,\
    #             xi_N,\
    #             noise_filename, solver_executable, print_iter_data, max_iter, \
    #             result_filename_header)
    # ------------------------------------
    # print("\nTesting CompoStep_TR_MA with black-box model")
    # result_filename_header = result_filename_folder + "CompoStep_TR_Blackbox_"
    # compo_step_TR_blackbox_model(perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,\
    #             xi_N,\
    #             noise_filename, solver_executable, print_iter_data, max_iter, \
    #             result_filename_header)
    # ------------------------------------
    print("\nTesting CompoStep_TR_MA with black-box model2")
    result_filename_header = result_filename_folder + "CompoStep_TR_Blackbox2_"
    compo_step_TR_blackbox_model2(perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius, \
                                 xi_N, \
                                 noise_filename, solver_executable, print_iter_data, max_iter, \
                                 result_filename_header)


def plot_profile():
    max_iter=20
    global_marker_size = 2
    linewidth=1
    compo_step_plant_data = pandas.read_csv("data/try_pso_uncon/CompoStep_TR_Pyomo_plant_data.txt", \
                                      index_col=0, header=0, sep='\t')
    compo_step_plant_data2 = pandas.read_csv("data/try_pso_uncon/CompoStep_TR_Blackbox_plant_data.txt", \
                                            index_col=0, header=0, sep='\t')
    compo_step_plant_data3 = pandas.read_csv("data/try_pso_uncon/CompoStep_TR_Blackbox2_plant_data.txt", \
                                             index_col=0, header=0, sep='\t')

    compo_step_model_data = pandas.read_csv("data/try_pso_uncon/CompoStep_TR_Pyomo_model_data.txt", \
                                      index_col=0, header=0, sep='\t')
    compo_step_model_data2 = pandas.read_csv("data/try_pso_uncon/CompoStep_TR_Blackbox_model_data.txt", \
                                            index_col=0, header=0, sep='\t')
    compo_step_model_data3 = pandas.read_csv("data/try_pso_uncon/CompoStep_TR_Blackbox2_model_data.txt", \
                                             index_col=0, header=0, sep='\t')

    compo_step_input_data = pandas.read_csv("data/try_pso_uncon/CompoStep_TR_Pyomo_input_data.txt", \
                                            index_col=0, header=0, sep='\t')
    compo_step_input_data2 = pandas.read_csv("data/try_pso_uncon/CompoStep_TR_Blackbox_input_data.txt", \
                                            index_col=0, header=0, sep='\t')
    compo_step_input_data3 = pandas.read_csv("data/try_pso_uncon/CompoStep_TR_Blackbox2_input_data.txt", \
                                             index_col=0, header=0, sep='\t')

    fig = plt.figure(figsize=(6,9))
    plt.subplot(511)
    # optimal = -210.33 * np.ones(max_iter + 1)
    # plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
    #          linestyle='--')
    plt.plot(range(1,max_iter+1), compo_step_plant_data.loc[1:max_iter, 'cost'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), compo_step_plant_data2.loc[1:max_iter, 'cost'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_plant_data3.loc[1:max_iter, 'cost'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel("plant cost")
    plt.subplot(512)
    # optimal = 0 * np.ones(max_iter + 1)
    # plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
    #          linestyle='--')
    # plt.plot(range(1, max_iter + 1), compo_step_plant_data.loc[1:max_iter, 'con'], \
    #          marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    # plt.plot(range(1, max_iter + 1), compo_step_plant_data2.loc[1:max_iter, 'con'], \
    #          marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel("plant constraints")
    plt.subplot(513)
    # optimal = 3.887 * np.ones(max_iter + 1)
    # plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
    #          linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_input_data.loc[1:max_iter, 'Fb1'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), compo_step_input_data.loc[1:max_iter, 'Fb2'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth * 2)
    plt.plot(range(1, max_iter + 1), compo_step_input_data2.loc[1:max_iter, 'Fb1'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_input_data2.loc[1:max_iter, 'Fb2'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_input_data3.loc[1:max_iter, 'Fb1'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_input_data3.loc[1:max_iter, 'Fb2'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"Fb")
    plt.subplot(514)
    # optimal = 9.369 * np.ones(max_iter + 1)
    # plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
    #          linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_input_data.loc[1:max_iter, 'Tr1'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth * 2)
    plt.plot(range(1, max_iter + 1), compo_step_input_data.loc[1:max_iter, 'Tr2'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth * 2)
    plt.plot(range(1, max_iter + 1), compo_step_input_data.loc[1:max_iter, 'Tr3'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth * 2)
    plt.plot(range(1, max_iter + 1), compo_step_input_data2.loc[1:max_iter, 'Tr1'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_input_data2.loc[1:max_iter, 'Tr2'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_input_data2.loc[1:max_iter, 'Tr3'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_input_data3.loc[1:max_iter, 'Tr1'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_input_data3.loc[1:max_iter, 'Tr2'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_input_data3.loc[1:max_iter, 'Tr3'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"Tr")
    plt.subplot(515)
    plt.plot(range(1, max_iter + 1), compo_step_model_data.loc[1:max_iter, 'tr'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), compo_step_model_data2.loc[1:max_iter, 'tr'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_model_data3.loc[1:max_iter, 'tr'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"tr")
    plt.xlabel("#iteration")

    plt.savefig("pic/try_pso/profile_uncon", dpi=600)
    plt.close()

if __name__ == "__main__":
    # generate_noise_file()
    do_test()
    plot_profile()