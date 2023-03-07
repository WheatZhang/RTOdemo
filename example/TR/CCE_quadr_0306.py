import os
import matplotlib.pyplot as plt
from example.TR.draw_lib import plot_contour
import pandas
import numpy
from rtolib.model.quadratic_con import quadratic_con_problem_plant, quadratic_con_problem_model,\
    quadratic_con_problem_description, quadr_obj_con_wrong_curv_model,\
    quadr_con_wrong_curv_model,quadr_obj_wrong_curv_model,quadr_obj_wrong_curv_model2,\
    quadr_con_wrong_curv_model2
from rtolib.core import NoiseGenerator, ModifierType,\
                    SimpleFiniteDiffPerturbation
from rtolib.core.algo import ModifierAdaptationPenaltyTR, \
    ModifierAdaptationCompoStepTR, ModifierAdaptationMaxTR
import copy
from rtolib.util.misc import save_iteration_data_in_dict


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
         'size': 17/font_factor,
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

def original_MA_with_max_tr(model_name, perturbation_stepsize, starting_point, max_trust_radius,\
                            filtering_factor,\
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
    problem_description = copy.deepcopy(quadratic_con_problem_description)

    rto_algorithm = ModifierAdaptationMaxTR()
    options = {
        "homotopy_simulation": False,
        "homotopy_optimization": False,
        "noise_adding_fashion": "added",
        "filtering_factor": filtering_factor,
        "max_trust_radius": max_trust_radius,
    }
    rto_algorithm.set_algorithm_option(options)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    if model_name == "norminal":
        model = quadratic_con_problem_model()
    elif model_name == "obj_wrong_curv":
        model = quadr_obj_wrong_curv_model()
    elif model_name == "obj_partial_wrong_curv":
        model = quadr_obj_wrong_curv_model2()
    elif model_name == "con_wrong_curv":
        model = quadr_con_wrong_curv_model()
    elif model_name == "con_wrong_curv2":
        model = quadr_con_wrong_curv_model2()
    elif model_name == "obj_con_wrong_curv":
        model = quadr_obj_con_wrong_curv_model()
    else:
        raise ValueError("Wrong model name.")

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=quadratic_con_problem_plant(),
        model=model,
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


def algo2_TR_MA(model_name, perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,\
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
    problem_description = copy.deepcopy(quadratic_con_problem_description)

    rto_algorithm = ModifierAdaptationPenaltyTR()
    options = {
        "homotopy_simulation": False,
        "homotopy_optimization": False,
        "noise_adding_fashion": "added",
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

    if model_name == "norminal":
        model = quadratic_con_problem_model()
    elif model_name == "obj_wrong_curv":
        model = quadr_obj_wrong_curv_model()
    elif model_name == "obj_partial_wrong_curv":
        model = quadr_obj_wrong_curv_model2()
    elif model_name == "con_wrong_curv":
        model = quadr_con_wrong_curv_model()
    elif model_name == "obj_con_wrong_curv":
        model = quadr_obj_con_wrong_curv_model()
    else:
        raise ValueError("Wrong model name.")

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=quadratic_con_problem_plant(),
        model=model,
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

def algo1_TR_MA(model_name, perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,xi_N,\
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
    problem_description = copy.deepcopy(quadratic_con_problem_description)

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
    }
    rto_algorithm.set_algorithm_option(options)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    if model_name == "norminal":
        model = quadratic_con_problem_model()
    elif model_name == "obj_wrong_curv":
        model = quadr_obj_wrong_curv_model()
    elif model_name == "obj_partial_wrong_curv":
        model = quadr_obj_wrong_curv_model2()
    elif model_name == "con_wrong_curv":
        model = quadr_con_wrong_curv_model()
    elif model_name == "con_wrong_curv2":
        model = quadr_con_wrong_curv_model2()
    elif model_name == "obj_con_wrong_curv":
        model = quadr_obj_con_wrong_curv_model()
    else:
        raise ValueError("Wrong model name.")

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=quadratic_con_problem_plant(),
        model=model,
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


def do_test_global_convergence():
    u1_0 = 2
    u2_0 = -2
    r_0 = 1
    sigma = 10
    xi_N = 0.5
    r_max = 2
    filtering_factor=0.1
    # ------------------------------------
    noise_filename = "noise/noise_0.txt"
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    result_filename_folder = "data/CCE_quadr_0306/global_conv/"
    if not os.path.exists(result_filename_folder):
        os.makedirs(result_filename_folder)

    # ------------------------------------
    starting_point = {
        "u1": u1_0,
        "u2": u2_0,
    }

    # ------------------------------------
    perturbation_stepsize = {
        "u1": 0.01,
        "u2": 0.01,
    }
    # ------------------------------------
    print_iter_data = False
    max_iter = 20
    initial_trust_radius = r_0
    sigma = sigma
    xi_N = xi_N


    # ------------------------------------
    for model_name in ["norminal", "obj_partial_wrong_curv", "con_wrong_curv2"]:
        print("\nTesting Original_MA_with_Rmax")
        result_filename_header = result_filename_folder + model_name+"_MA_"
        original_MA_with_max_tr(model_name, perturbation_stepsize, starting_point, r_max, filtering_factor, \
                    noise_filename, solver_executable, print_iter_data, max_iter, \
                    result_filename_header)

        # ------------------------------------
        print("\nTesting CompoStep_TR_MA")
        result_filename_header = result_filename_folder + model_name+"_CompoStep_"
        algo1_TR_MA(model_name, perturbation_stepsize, starting_point, sigma, initial_trust_radius,r_max, xi_N, \
                    noise_filename, solver_executable, print_iter_data, max_iter, \
                    result_filename_header)

def do_test_sigma():
    f_u1_0 = 2
    f_u2_0 = -2
    r_0 = 1
    xi_N = 0.5
    r_max = 2
    filtering_factor=0.1
    # ------------------------------------
    noise_filename = "noise/noise_0.txt"
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    result_filename_folder = "data/CCE_quadr_0306/penalty_coeff/"
    if not os.path.exists(result_filename_folder):
        os.makedirs(result_filename_folder)

    # ------------------------------------
    feasible_starting_point = {
        "u1": f_u1_0,
        "u2": f_u2_0,
    }
    # ------------------------------------
    perturbation_stepsize = {
        "u1": 0.01,
        "u2": 0.01,
    }
    # ------------------------------------
    print_iter_data = False
    max_iter = 20
    initial_trust_radius = r_0
    xi_N = xi_N
    model_name = "norminal"

    # ------------------------------------
    name=["LP","MP","SP"]
    for no,sigma in enumerate([1,0.1,0.01]):
        print("\nTesting Penalty_MA")
        result_filename_header = result_filename_folder + name[no]+"Penalty_"
        algo2_TR_MA(model_name, perturbation_stepsize, feasible_starting_point, sigma, initial_trust_radius, r_max, \
                    noise_filename, solver_executable, print_iter_data, max_iter, \
                    result_filename_header)

        print("\nTesting CompoStep_TR_MA")
        result_filename_header = result_filename_folder + name[no]+"CompoStep_"
        algo1_TR_MA(model_name, perturbation_stepsize, feasible_starting_point, sigma, initial_trust_radius,r_max, xi_N, \
                    noise_filename, solver_executable, print_iter_data, max_iter, \
                    result_filename_header)

def do_test_xi():
    inf_u1_0 = -1
    inf_u2_0 = 0
    r_0 = 1
    sigma = 1
    r_max = 2
    # ------------------------------------
    noise_filename = "noise/noise_0.txt"
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    result_filename_folder = "data/CCE_quadr_0306/xi/"
    if not os.path.exists(result_filename_folder):
        os.makedirs(result_filename_folder)

    # ------------------------------------
    infeasible_starting_point = {
        "u1": inf_u1_0,
        "u2": inf_u2_0,
    }

    # ------------------------------------
    perturbation_stepsize = {
        "u1": 0.01,
        "u2": 0.01,
    }
    # ------------------------------------
    print_iter_data = False
    max_iter = 20
    initial_trust_radius = r_0
    model_name = "norminal"

    # ------------------------------------
    print("\nTesting Penalty_MA")
    result_filename_header = result_filename_folder + name[no] + "Penalty_"
    algo2_TR_MA(model_name, perturbation_stepsize, infeasible_starting_point, sigma, initial_trust_radius, r_max, \
                noise_filename, solver_executable, print_iter_data, max_iter, \
                result_filename_header)

    name=["3","6","9"]
    for no,xi_N in enumerate([0.3,0.6,0.9]):
        print("\nTesting CompoStep_TR_MA")
        result_filename_header = result_filename_folder + name[no]+"CompoStep_"
        algo1_TR_MA(model_name, perturbation_stepsize, feasible_starting_point, sigma, initial_trust_radius,r_max, xi_N, \
                    noise_filename, solver_executable, print_iter_data, max_iter, \
                    result_filename_header)

def do_test_trust_radius():
    f_u1_0 = 2
    f_u2_0 = -2
    r_0 = 1
    xi_N = 0.5
    sigma = 1
    r_max = 2
    # ------------------------------------
    noise_filename = "noise/noise_0.txt"
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    result_filename_folder = "data/CCE_quadr_0306/trust_radius/"
    if not os.path.exists(result_filename_folder):
        os.makedirs(result_filename_folder)

    # ------------------------------------
    infeasible_starting_point = {
        "u1": f_u1_0,
        "u2": f_u2_0,
    }

    # ------------------------------------
    perturbation_stepsize = {
        "u1": 0.01,
        "u2": 0.01,
    }
    # ------------------------------------
    print_iter_data = False
    max_iter = 20
    initial_trust_radius = r_0
    model_name = "obj_partial_wrong_curv"

    name=["1","2","5","10"]
    for no,coeff in enumerate([1,2,5,10]):
        print("\nTesting CompoStep_TR_MA")
        result_filename_header = result_filename_folder + name[no]+"CompoStep_"
        algo1_TR_MA(model_name, perturbation_stepsize, feasible_starting_point,\
                    sigma, initial_trust_radius/coeff,r_max/coeff, xi_N, \
                    noise_filename, solver_executable, print_iter_data, max_iter, \
                    result_filename_header)

def generate_all_data():
    do_test_global_convergence()
    do_test_sigma()
    do_test_xi()
    do_test_trust_radius()


def draw_profile_global_convergence():
    max_iter = 20
    global_marker_size = 2
    linewidth = 1
    fig = plt.figure(figsize=(21 * pic_constant, 10 * pic_constant))
    for ax_no in [231, 232, 233]:
        plt.subplot(ax_no)
        optimal = 1.453659e-01 * numpy.ones(max_iter + 1)
        plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='blue',
                 linestyle='--')
        plt.ylim([-1, 5])
        plt.grid(ls='--', axis='y')
    for ax_no in [234, 235, 236]:
        plt.subplot(ax_no)
        optimal = 0 * numpy.ones(max_iter + 1)
        plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='blue',
                 linestyle='--')
        plt.ylim([-2, 2])
        plt.grid(ls='--', axis='y')
    for no,prefix in enumerate(["norminal", "obj_partial_wrong_curv", "con_wrong_curv2"]):
        compo_step_data = pandas.read_csv("data/CCE_quadr_0306/global_conv/" + prefix + "_CompoStep_plant_data.txt", \
                                          index_col=0, header=0, sep='\t')
        origin_ma_data = pandas.read_csv("data/CCE_quadr_0306/global_conv/" + prefix + "_MA_plant_data.txt", \
                                          index_col=0, header=0, sep='\t')
        plt.subplot(2, 3, no + 1)
        plt.plot(range(1, max_iter + 1), compo_step_data.loc[1:max_iter, 'cost'], \
                 marker='o', c='black', markersize=global_marker_size, linewidth=linewidth, label="Composite-step")
        plt.plot(range(1, max_iter + 1), origin_ma_data.loc[1:max_iter, 'cost'], \
                 marker='o', c='red', markersize=global_marker_size / 3, linewidth=linewidth / 3, label="Modifier adaptation")
        if no == 0:
            plt.ylabel("cost", font_label)
            plt.title(r"Model with correct curvatures", fontdict=font_title)
        elif no == 1:
            plt.title(r"Model with wrong objective curvature", fontdict=font_title)
        elif no == 2:
            plt.title(r"Model with wrong constraint curvature", fontdict=font_title)
            plt.legend(loc='upper right', prop=font_legend)
        plt.xticks([0, 5, 10, 15, 20])
        plt.yticks([-1, 1, 3, 5])
        plt.tick_params(labelsize=global_tick_size)

        plt.subplot(2, 3, no + 4)
        plt.plot(range(1, max_iter + 1), compo_step_data.loc[1:max_iter, 'con'], \
                 marker='o', c='black', markersize=global_marker_size, linewidth=linewidth)
        plt.plot(range(1, max_iter + 1), origin_ma_data.loc[1:max_iter, 'con'], \
                 marker='o', c='red', markersize=global_marker_size / 3, linewidth=linewidth / 3)
        if no == 0:
            plt.ylabel("constraint", font_label)
        plt.xlabel("RTO iteration", font_label)
        plt.xticks([0, 5, 10, 15, 20])
        plt.tick_params(labelsize=global_tick_size)

    for i in range(1, 7):
        plt.subplot(2, 3, i)
        ax = plt.gca()
        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontname('Times New Roman') for label in labels]

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)

    plt.savefig("pic/CCE_quadr_0306/global_convengence_profile", dpi=600)
    plt.close()

def draw_profile_sigma():
    max_iter = 20
    global_marker_size = 2
    linewidth = 1
    fig = plt.figure(figsize=(21 * pic_constant, 10 * pic_constant))
    for ax_no in [231,232,233]:
        plt.subplot(ax_no)
        optimal = 1.453659e-01 * numpy.ones(max_iter + 1)
        plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='blue',
                 linestyle='--')
        plt.ylim([-1,5])
        plt.grid(ls='--',axis='y')
    for ax_no in [234,235,236]:
        plt.subplot(ax_no)
        optimal = 0 * numpy.ones(max_iter + 1)
        plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='blue',
                 linestyle='--')
        plt.ylim([-2, 2])
        plt.grid(ls='--',axis='y')
    for no,prefix in enumerate(["LPF","MPF","SPF"]):
        compo_step_data = pandas.read_csv("data/CCE_quadr_0306/penalty_coeff/"+ prefix+"_CompoStep_plant_data.txt", \
                                          index_col=0, header=0, sep='\t')
        penalty_tr_data = pandas.read_csv("data/CCE_quadr_0306/penalty_coeff/"+ prefix+"_Penalty_plant_data.txt", \
                                       index_col=0, header=0, sep='\t')
        plt.subplot(2,3,no+1)
        plt.plot(range(1, max_iter + 1), compo_step_data.loc[1:max_iter, 'cost'], \
                 marker='o', c='black', markersize=global_marker_size, linewidth=linewidth,label="Composite-step")
        plt.plot(range(1, max_iter + 1), penalty_tr_data.loc[1:max_iter, 'cost'], \
                 marker='o', c='red', markersize=global_marker_size/3, linewidth=linewidth/3,label="Penalty")
        if no ==0:
            plt.ylabel("cost",font_label)
            plt.title(r"$\sigma=1$", fontdict=font_title)
        elif no == 1:
            plt.title(r"$\sigma=0.1$", fontdict=font_title)
        elif no == 2:
            plt.title(r"$\sigma=0.01$", fontdict=font_title)
            plt.legend(loc='upper right', prop=font_legend)
        plt.xticks([0, 5, 10, 15, 20])
        plt.yticks([-1,1,3,5])
        plt.tick_params(labelsize=global_tick_size)

        plt.subplot(2,3, no + 4)
        plt.plot(range(1, max_iter + 1), compo_step_data.loc[1:max_iter, 'con'], \
                 marker='o', c='black', markersize=global_marker_size, linewidth=linewidth)
        plt.plot(range(1, max_iter + 1), penalty_tr_data.loc[1:max_iter, 'con'], \
                 marker='o', c='red', markersize=global_marker_size/3, linewidth=linewidth/3)
        if no == 0:
            plt.ylabel("constraint",font_label)
        plt.xlabel("RTO iteration",font_label)
        plt.xticks([0, 5, 10, 15, 20])
        plt.tick_params(labelsize=global_tick_size)

    for i in range(1,7):
        plt.subplot(2,3,i)
        ax = plt.gca()
        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontname('Times New Roman') for label in labels]

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("pic/CCE_quadr_0306/penalty_comparison", dpi=600)
    plt.close()

if __name__ == "__main__":
    # do_test_global_convergence()
    # draw_profile_global_convergence()
    # do_test_penalty()
    draw_profile_penalty()