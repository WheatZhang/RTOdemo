import numpy as np
import pickle
from rtolib.core.dc_cpwl_model import QuadraticBoostedDCCPWLFunction
from rtolib.core.algo.dccpwl_ma import DCCPWL_ModifierAdaptation
from rtolib.model.wo_reactor import RTO_Plant_WO_reactor, RTO_Mismatched_WO_reactor,\
    default_WOR_description
from rtolib.core import NoiseGenerator, ModifierType,\
                    SimpleFiniteDiffPerturbation
from rtolib.core.algo import ModifierAdaptation
import os
import copy
from rtolib.util.misc import save_iteration_data_in_dict
import matplotlib.pyplot as plt
import pandas
from matplotlib.pyplot import MultipleLocator

class UnconstrainedWO_obj_CPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "model/Unconstrained_WO_2d_processed_output_0.bin"
        self.load_from_cpwl_model_file(cpwl_file)
        self.input_variables = ['Tr', 'Fb']
        self.parameters = []

class UnconstrainedWO_obj_QCPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "qmodel/Unconstrained_WO_2d_cpwl_processed_output_0.bin"
        quadr_file = "qmodel/Unconstrained_WO_2d_quadr_processed_output_0.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['Tr', 'Fb']
        self.parameters = []

def do_test_CPWL_MA(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, nlp_solver_executable, qcqp_solver_executable,\
                     print_iter_data, max_iter,\
               result_filename_header):
    problem_description = copy.deepcopy(default_WOR_description)
    problem_description.symbol_list['CV']=('cost',)

    rto_algorithm = DCCPWL_ModifierAdaptation()
    options = {
        "homotopy_simulation": True,
        "homotopy_optimization": True,
        "filtering_factor": filtering_factor,
        "noise_adding_fashion": "integrated",
    }
    rto_algorithm.set_algorithm_option(options)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=RTO_Plant_WO_reactor(),
        model_dc_cpwl_functions={"cost":UnconstrainedWO_obj_CPWL()},
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        nlp_solver_executable=nlp_solver_executable,
        qcqp_solver_executable=qcqp_solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        parameter_set=[],
        )

    homotopy_var=['Fb','Tr']
    rto_algorithm.plant_simulator.homotopy_var = homotopy_var

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

def do_test_QCPWL_MA(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, nlp_solver_executable, qcqp_solver_executable,\
                     print_iter_data, max_iter,\
               result_filename_header):
    problem_description = copy.deepcopy(default_WOR_description)
    problem_description.symbol_list['CV']=('cost',)

    rto_algorithm = DCCPWL_ModifierAdaptation()
    options = {
        "homotopy_simulation": True,
        "homotopy_optimization": True,
        "filtering_factor": filtering_factor,
        "noise_adding_fashion": "integrated",
    }
    rto_algorithm.set_algorithm_option(options)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=RTO_Plant_WO_reactor(),
        model_dc_cpwl_functions={"cost":UnconstrainedWO_obj_QCPWL()},
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        nlp_solver_executable=nlp_solver_executable,
        qcqp_solver_executable=qcqp_solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        parameter_set=[],
        )

    homotopy_var=['Fb','Tr']
    rto_algorithm.plant_simulator.homotopy_var = homotopy_var

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
    problem_description.symbol_list['CV']=('cost',)

    rto_algorithm = ModifierAdaptation()
    options = {
        "homotopy_simulation": True,
        "homotopy_optimization": True,
        "filtering_factor": filtering_factor,
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

def algo_test_main():
    # ------------------------------------
    noise_filename = "noise/uncon-wo-noise0.txt"
    nlp_solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    qcqp_solver_executable = None
    result_filename_folder="data/uncon-wo/"
    if not os.path.exists(result_filename_folder):
        os.makedirs(result_filename_folder)

    # ------------------------------------
    filtering_factor = 0.5
    starting_point = {
        "Fb": 4,
        "Tr": 75,
    }

    # ------------------------------------
    perturbation_stepsize = {
        "Fb": 0.01,#0.2,
        "Tr": 0.1,#2,
    }
    # ------------------------------------
    print_iter_data = False
    max_iter = 20

    # ------------------------------------
    print("\nTesting MA")
    result_filename_header = result_filename_folder + "MA_"
    do_test_MA(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, nlp_solver_executable, print_iter_data, max_iter, \
               result_filename_header)

    print("\nTesting CPWL-MA")
    result_filename_header = result_filename_folder + "CPWL_MA_"
    do_test_CPWL_MA(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, nlp_solver_executable, qcqp_solver_executable,\
               print_iter_data, max_iter, \
               result_filename_header)

    print("\nTesting QCPWL-MA")
    result_filename_header = result_filename_folder + "QCPWL_MA_"
    do_test_QCPWL_MA(perturbation_stepsize, starting_point, filtering_factor, \
                     noise_filename, nlp_solver_executable, qcqp_solver_executable, \
                     print_iter_data, max_iter, \
                     result_filename_header)

def draw_overall_pic(folder, pic_filename):
    # pic_filename="pic/overall_pic0-0.png"
    # =======================================
    #            Load Data
    # =======================================
    # folder="data/0-0/"
    MA_plant_data = pandas.read_csv(folder + "MA_plant_data.txt", sep='\t', index_col=0, header=0)
    MA_input_data = pandas.read_csv(folder + "MA_input_data.txt", sep='\t', index_col=0, header=0)
    CPWL_MA_plant_data = pandas.read_csv(folder + "CPWL_MA_plant_data.txt", sep='\t', index_col=0, header=0)
    CPWL_MA_input_data = pandas.read_csv(folder + "CPWL_MA_input_data.txt", sep='\t', index_col=0, header=0)
    QCPWL_MA_plant_data = pandas.read_csv(folder + "QCPWL_MA_plant_data.txt", sep='\t', index_col=0, header=0)
    QCPWL_MA_input_data = pandas.read_csv(folder + "QCPWL_MA_input_data.txt", sep='\t', index_col=0, header=0)
    
    max_iter=20
    # #----------Optimal-------------
    ori_plant_cost = 190.98 * np.ones(max_iter+1)
    ori_FB = 4.7874 * np.ones(max_iter+1)
    ori_TR = 89.704 * np.ones(max_iter+1)
    # =======================================
    #            Draw Pictures
    # =======================================
    global_font_size = 20
    global_tick_size = 20
    global_linewidth = 1
    global_point_size = 4
    pic_constant = 0.39370
    font_factor = np.sqrt(1 / pic_constant)
    global_font_size = global_font_size / font_factor
    global_tick_size = global_tick_size / font_factor

    fig = plt.figure(figsize=(13 * 0.39370, 26 * 0.39370))
    font1 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': global_tick_size,
             }
    font2 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': global_font_size,
             }
    font3 = {'family': 'Times New Roman',
             'weight': 'bold',
             'size': global_font_size,
             }
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    plt.subplot(3, 1, 1)
    Optimal, = plt.plot(range(max_iter+1), ori_plant_cost, linewidth=global_linewidth, label='Optimal', color='gray',
                        linestyle='--')
    MA, = plt.plot(range(max_iter), MA_plant_data.loc[1:(max_iter+1),"cost"].multiply(-1), linewidth=global_linewidth, label='MA', linestyle='-',
                   marker='o', markersize=global_point_size)
    CPWL_MA, = plt.plot(range(max_iter), CPWL_MA_plant_data.loc[1:(max_iter+1),"cost"].multiply(-1), linewidth=global_linewidth, label='CPWL MA', linestyle='-',
                    marker='o', markersize=global_point_size)
    QCPWL_MA, = plt.plot(range(max_iter), QCPWL_MA_plant_data.loc[1:(max_iter+1),"cost"].multiply(-1), linewidth=global_linewidth, label='Quadratic+CPWL MA',
                      linestyle='-', marker='o', markersize=global_point_size)

    plt.legend(handles=[MA, CPWL_MA, QCPWL_MA, Optimal], prop=font1)
    plt.legend(loc='lower right',
               fancybox=False,
               edgecolor='k',
               fontsize=global_point_size,
               shadow=False,
               facecolor='w',
               framealpha=1.0,
               prop={'family': 'Times New Roman', 'size': global_tick_size})
    plt.ylabel('(a) Plant cost', font2)
    plt.xticks([0, 3, 6, 9, 12, 15, 18])
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.tick_params(labelsize=global_tick_size)
    ax = plt.gca()
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    plt.tick_params(top=['on', 'true'], right=['on', 'true'], which='both')
    plt.xlim(0, 20)
    plt.ylim(145, 195)

    plt.subplot(3, 1, 2)
    OptimalFB, = plt.plot(range(max_iter + 1), ori_FB, linewidth=global_linewidth, label='Optimal', color='gray',
                          linestyle='--')
    MA, = plt.plot(range(max_iter), MA_input_data.loc[range(max_iter),"Fb"], linewidth=global_linewidth, label='MA', linestyle='-',
                   marker='o', markersize=global_point_size)
    CPWL_MA, = plt.plot(range(max_iter), CPWL_MA_input_data.loc[range(max_iter),"Fb"], linewidth=global_linewidth, label='CPWL MA', linestyle='-',
                    marker='o', markersize=global_point_size)
    QCPWL_MA, = plt.plot(range(max_iter), QCPWL_MA_input_data.loc[range(max_iter),"Fb"], linewidth=global_linewidth, label='Quadratic+CPWL MA',
                      linestyle='-', marker='o', markersize=global_point_size)

    plt.ylabel('(b) Flowrates of B, $F_B$ (kg/s)', font2)
    plt.xticks([0, 3, 6, 9, 12, 15, 18])
    plt.yticks([4, 4.5, 5, 5.5])
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.tick_params(labelsize=global_tick_size)
    ax = plt.gca()
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    plt.tick_params(top=['on', 'true'], right=['on', 'true'], which='both')
    plt.xlim(0, 20)
    plt.ylim(3.9, 5.6)

    plt.subplot(3, 1, 3)
    OptimalTR, = plt.plot(range(max_iter + 1), ori_TR, linewidth=global_linewidth, label='Optimal', color='gray',
                          linestyle='--')
    MA, = plt.plot(range(max_iter), MA_input_data.loc[range(max_iter),"Tr"], linewidth=global_linewidth, label='MA', linestyle='-',
                   marker='o', markersize=global_point_size)
    CPWL_MA, = plt.plot(range(max_iter), CPWL_MA_input_data.loc[range(max_iter),"Tr"], linewidth=global_linewidth, label='CPWL MA', linestyle='-',
                    marker='o', markersize=global_point_size)
    QCPWL_MA, = plt.plot(range(max_iter), QCPWL_MA_input_data.loc[range(max_iter),"Tr"], linewidth=global_linewidth, label='Quadratic+CPWL MA',
                      linestyle='-', marker='o', markersize=global_point_size)

    # plt.legend(handles=[OptimalTR, MA, CPWL_MA, PA, QCPWL_MA, ISOQCPWL_MA], prop=font1)
    plt.xlabel('RTO Iterations', font2)
    plt.ylabel('(c) CSTR TemQCPWL_MArature, $T_R$ (℃)', font2)
    plt.xticks([0, 3, 6, 9, 12, 15, 18])
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    ax = plt.gca()
    plt.tick_params(labelsize=global_tick_size)
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    plt.tick_params(top=['on', 'true'], right=['on', 'true'], which='both')
    plt.xlim(0, 20)
    plt.ylim(74, 101)
    y_major_locator = MultipleLocator(5)
    ax = plt.gca()
    ax.yaxis.set_major_locator(y_major_locator)

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig(pic_filename, dpi=600)
    plt.close()


if __name__ == "__main__":
    # algo_test_main()
    draw_overall_pic("data/uncon-wo/","pic/uncon-wo")