from rtolib.model.wo_reactor import RTO_Plant_WO_reactor, RTO_Mismatched_WO_reactor,\
    default_WOR_description, ZeroModel_WO_reactor
from rtolib.core import NoiseGenerator, ModifierType,\
                    SimpleFiniteDiffPerturbation
from rtolib.core.algo import ModifierAdaptationCompoStepTR
import copy
from rtolib.util.misc import save_iteration_data_in_dict
import matplotlib.pyplot as plt
import pandas
import os
import numpy

global_parameter={
        "eta1": 0.01,
        "eta2": 0.9,
        "gamma1": 0.5,
        "gamma2": 1,
        "gamma3": 2,
}

global_linewidth = 1
global_point_size= 4
pic_constant = 0.39370
font_factor = numpy.sqrt(1/pic_constant)
font_legend = {'family': 'Times New Roman',
         'weight': 'normal',
         'size': 17/font_factor/1.2,
         }
font_title = {'family': 'Times New Roman',
              'size': 17/font_factor,
              'weight': 'normal'
              }
font_label = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 17/font_factor,
             }
global_tick_size=17/font_factor
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'


def generate_noise_file():
    noise_filename = "noise/uncon_wo_noise0.txt"
    cost_noise_level=0
    composition_noise_level=0
    noise_level = {
            'cost': cost_noise_level,
            'XFr_A': composition_noise_level,
            'XFr_B': composition_noise_level,
            'XFr_E': composition_noise_level,
            'XFr_P': composition_noise_level,
            'XFr_G': composition_noise_level,
    }
    max_iter = 200
    noise_generator = NoiseGenerator(noise_level)
    noise_generator.generate_noise(max_iter, 5)
    noise_generator.save_noise(noise_filename)


def do_test_MA(model_name, perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,xi_N,\
               noise_filename, solver_executable, print_iter_data, max_iter,\
               result_filename_header):

    problem_description = copy.deepcopy(default_WOR_description)
    problem_description.symbol_list['CV']=('cost',)

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
        "adaptive_sigma": True,
        "max_trust_radius": max_trust_radius,
        "initial_trust_radius": initial_trust_radius,
    }
    rto_algorithm.set_algorithm_option(options)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    if model_name == "ODM":
        model = RTO_Mismatched_WO_reactor()
    elif model_name == "ZERO":
        model = ZeroModel_WO_reactor()

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=RTO_Plant_WO_reactor(),
        model=model,
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        solver_executable=solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
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


def get_simulation_result():
    # ------------------------------------
    noise_filename = "noise/uncon_wo_noise0.txt"
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    result_filename_folder = "data/zero_model_unconwo/"

    # ------------------------------------
    starting_point = {
        "Fb": 4,
        "Tr": 75,
    }
    # ------------------------------------
    perturbation_stepsize = {
        "Fb": 0.2,
        "Tr": 2,
    }

    max_trust_radius = 4
    xi_N = 0.5
    sigma = 1
    initial_trust_radius = 2
    # ------------------------------------
    print_iter_data = False
    max_iter = 20

    # ------------------------------------
    print("\nTesting MA-ODM")
    result_filename_header = result_filename_folder + "MA_ODM_"
    do_test_MA("ODM",perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius, xi_N, \
               noise_filename, solver_executable, print_iter_data, max_iter, \
               result_filename_header)

    print("\nTesting MA-ZEROMODEL")
    result_filename_header = result_filename_folder + "MA_ZERO_"
    do_test_MA("ZERO", perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius, xi_N, \
               noise_filename, solver_executable, print_iter_data, max_iter, \
               result_filename_header)

def draw_profile():
    max_iter = 20
    global_marker_size = 2
    linewidth = 1
    fig = plt.figure(figsize=(21 * pic_constant, 10 * pic_constant))

    optimal = -190.98 * numpy.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth/3, label='Optimal', color='black',
             linestyle='--')
    plt.ylim([-100, -200])
    plt.grid(ls='--', axis='y')

    color = {
        "MA_ODM": "black",
        "MA_ZERO": "red",
    }
    for no,prefix in enumerate(["MA_ODM", "MA_ZERO"]):
        plant_data = pandas.read_csv("data/zero_model_unconwo/" + prefix + "_plant_data.txt", \
                                          index_col=0, header=0, sep='\t')
        plt.plot(range(1, max_iter + 1), plant_data.loc[1:max_iter, 'cost'], \
                 marker='o', c=color[prefix], markersize=global_marker_size, linewidth=linewidth, label="Composite-step")

    plt.xticks([0, 5, 10, 15, 20])
    plt.tick_params(labelsize=global_tick_size)

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)

    plt.savefig("pic/zero_model/uncon_wo", dpi=600)
    plt.close()


if __name__== "__main__":
    # generate_noise_file()
    get_simulation_result()
    draw_profile()