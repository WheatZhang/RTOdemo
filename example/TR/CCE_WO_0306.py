import os
import matplotlib.pyplot as plt
import pandas
import numpy
from rtolib.model.wo_reactor_con import RTO_Plant_WO_reactor, RTO_Mismatched_WO_reactor,\
    RTO_Mismatched_WO_reactor_QC, default_WOR_description
from rtolib.core import NoiseGenerator, ModifierType,\
                    SimpleFiniteDiffPerturbation
from rtolib.core.algo import ModifierAdaptationCompoStepTR, ModifierAdaptation
from rtolib.util.misc import save_iteration_data_in_dict
import copy

global_parameter={
        "eta1": 0.01,
        "eta2": 0.9,
        "gamma1": 0.5,
        "gamma2": 1,
        "gamma3": 2,
}

def compo_step_TR_MA(perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,xi_N,\
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
        "filtering_factor": 0.5,
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
        "filtering_factor": 1,
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
    result_filename_folder="data/CCE_wo_0306/"
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
    max_iter = 40
    initial_trust_radius = 0.1#0.1
    sigma = 100
    xi_N = 0.5
    max_trust_radius=0.5 #1


    print("\nTesting CompoStep_TR_MA")
    result_filename_header = result_filename_folder + "CompoStep_TR_MA_"
    compo_step_TR_MA(perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,\
                xi_N,\
                noise_filename, solver_executable, print_iter_data, max_iter, \
                result_filename_header)

    # ------------------------------------
    # print("\nTesting Original_MA")
    # result_filename_header = result_filename_folder + "Original_MA_"
    # original_MA(perturbation_stepsize, starting_point, \
    #             noise_filename, solver_executable, print_iter_data, max_iter, \
    #             result_filename_header)
    #
    # # ------------------------------------
    # print("\nTesting Original_MAQC")
    # result_filename_header = result_filename_folder + "Original_MAQC_"
    # original_MA_with_QC(perturbation_stepsize, starting_point, \
    #             noise_filename, solver_executable, print_iter_data, max_iter, \
    #             result_filename_header)

def plot_profile():
    max_iter=30
    global_marker_size = 2
    linewidth=1
    compo_step_plant_data = pandas.read_csv("data/CCE_wo_0306/CompoStep_TR_MA_plant_data.txt", \
                                      index_col=0, header=0, sep='\t')
    QC_ma_plant_data = pandas.read_csv("data/CCE_wo_0306/Original_MAQC_plant_data.txt", \
                                         index_col=0, header=0, sep='\t')
    original_ma_plant_data = pandas.read_csv("data/CCE_wo_0306/Original_MA_plant_data.txt", \
                                            index_col=0, header=0, sep='\t')

    compo_step_model_data = pandas.read_csv("data/CCE_wo_0306/CompoStep_TR_MA_model_data.txt", \
                                      index_col=0, header=0, sep='\t')
    QC_ma_model_data = pandas.read_csv("data/CCE_wo_0306/Original_MAQC_model_data.txt", \
                                            index_col=0, header=0, sep='\t')
    original_ma_model_data = pandas.read_csv("data/CCE_wo_0306/Original_MA_model_data.txt", \
                                             index_col=0, header=0, sep='\t')

    compo_step_input_data = pandas.read_csv("data/CCE_wo_0306/CompoStep_TR_MA_input_data.txt", \
                                            index_col=0, header=0, sep='\t')
    QC_ma_input_data = pandas.read_csv("data/CCE_wo_0306/Original_MAQC_input_data.txt", \
                                            index_col=0, header=0, sep='\t')
    original_ma_input_data = pandas.read_csv("data/CCE_wo_0306/Original_MA_input_data.txt", \
                                             index_col=0, header=0, sep='\t')

    fig = plt.figure(figsize=(6,12))
    plt.subplot(511)
    optimal = -210.33 * numpy.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1,max_iter+1), compo_step_plant_data.loc[1:max_iter, 'cost'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2,\
             label="composite-step TR")
    plt.plot(range(1, max_iter + 1), original_ma_plant_data.loc[1:max_iter, 'cost'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth, \
             label="MA")
    plt.plot(range(1, max_iter + 1), QC_ma_plant_data.loc[1:max_iter, 'cost'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth,\
             label="MA with convex model")
    plt.ylabel("plant cost")
    plt.legend()
    plt.subplot(512)
    optimal = 0 * numpy.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_plant_data.loc[1:max_iter, 'con'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), QC_ma_plant_data.loc[1:max_iter, 'con'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), original_ma_plant_data.loc[1:max_iter, 'con'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel("plant constraints")
    plt.subplot(513)
    optimal = 3.887 * numpy.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_input_data.loc[1:max_iter, 'Fa'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), QC_ma_input_data.loc[1:max_iter, 'Fa'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), original_ma_input_data.loc[1:max_iter, 'Fa'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"Fa")
    plt.subplot(514)
    optimal = 9.369 * numpy.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_input_data.loc[1:max_iter, 'Fb'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), QC_ma_input_data.loc[1:max_iter, 'Fb'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), original_ma_input_data.loc[1:max_iter, 'Fb'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"Fb")
    plt.subplot(515)
    optimal = 91.2 * numpy.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_input_data.loc[1:max_iter, 'Tr'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), QC_ma_input_data.loc[1:max_iter, 'Tr'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), original_ma_input_data.loc[1:max_iter, 'Tr'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"Tr")

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("pic/CCE_wo_0306/wo_profile", dpi=600)
    plt.close()

if __name__ == "__main__":
    do_test()
    plot_profile()