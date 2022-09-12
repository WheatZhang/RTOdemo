#!/usr/bin/env python
#-*- coding:utf-8 -*-
from rtolib.model.quadratic_con import quadratic_con_problem_plant, quadratic_con_problem_model,\
    quadratic_con_problem_description, quadr_obj_con_wrong_curv_model,\
    quadr_con_wrong_curv_model,quadr_obj_wrong_curv_model,quadr_obj_wrong_curv_model2
from rtolib.core import NoiseGenerator, ModifierType,\
                    SimpleFiniteDiffPerturbation
from rtolib.core.algo import ModifierAdaptationTR, ModifierAdaptationPenaltyTR, ModifierAdaptationCompoStepTR
import copy
from rtolib.util.misc import save_iteration_data_in_dict
import os


global_parameter={
        "eta1": 0.01,
        "eta2": 0.9,
        "gamma1": 0.5,
        "gamma2": 1,
        "gamma3": 2,
}

def generate_noise_file():
    noise_level = {
            'cost': 0,
            'con': 0,
    }
    max_iter = 200
    noise_generator = NoiseGenerator(noise_level)
    noise_generator.generate_noise(max_iter, 5)
    noise_generator.save_noise("noise/noise_0.txt")

def algo3_TR_MA(model_name, perturbation_stepsize, starting_point, initial_trust_radius, max_trust_radius,\
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

    rto_algorithm = ModifierAdaptationTR()
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
        "max_trust_radius": max_trust_radius,
        "initial_trust_radius": initial_trust_radius,
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


def do_batch_test_algo12(model_name, batch_label, param, u1_0, u2_0, r_0, r_max, sigma, xi_N):
    # ------------------------------------
    noise_filename = "noise/noise_0.txt"
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    result_filename_folder = "data/batch12/"+batch_label+"%d/"%param
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
    max_iter = 25
    initial_trust_radius = r_0
    sigma = sigma
    xi_N = xi_N

    # ------------------------------------
    print("\nTesting Penalty_TR_MA")
    result_filename_header = result_filename_folder + "Penalty_TR_MA_"
    algo2_TR_MA(model_name, perturbation_stepsize, starting_point, sigma, initial_trust_radius, r_max,\
                noise_filename, solver_executable, print_iter_data, max_iter, \
                result_filename_header)

    # ------------------------------------
    print("\nTesting CompoStep_TR_MA")
    result_filename_header = result_filename_folder + "CompoStep_TR_MA_"
    algo1_TR_MA(model_name, perturbation_stepsize, starting_point, sigma, initial_trust_radius,r_max, xi_N, \
                noise_filename, solver_executable, print_iter_data, max_iter, \
                result_filename_header)


def do_batch_test_algo13(model_name, batch_label, param, u1_0, u2_0, r_0, r_max, sigma, xi_N):
    # ------------------------------------
    noise_filename = "noise/noise_0.txt"
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    result_filename_folder = "data/batch13/"+batch_label+"%d/"%param
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
    max_iter = 25
    initial_trust_radius = r_0
    sigma = sigma
    xi_N = xi_N

    # ------------------------------------
    print("\nTesting InfAverse_TR_MA")
    result_filename_header = result_filename_folder + "InfAverse_TR_MA_"
    algo3_TR_MA(model_name, perturbation_stepsize, starting_point, initial_trust_radius, r_max,\
                noise_filename, solver_executable, print_iter_data, max_iter, \
                result_filename_header)

    # ------------------------------------
    print("\nTesting CompoStep_TR_MA")
    result_filename_header = result_filename_folder + "CompoStep_TR_MA_"
    algo1_TR_MA(model_name, perturbation_stepsize, starting_point, sigma, initial_trust_radius,r_max, xi_N, \
                noise_filename, solver_executable, print_iter_data, max_iter, \
                result_filename_header)


def do_batch_test_algo12_main(model_name):
    default_u1_0 = 1.1
    default_u2_0 = 0.1
    default_r_0 = 1
    default_sigma = 10
    default_xi_N = 0.5
    default_r_max = 4

    if model_name == "norminal":
        batch_prefix = "N_"
    elif model_name == "obj_wrong_curv":
        batch_prefix = "O1_"
    elif model_name == "obj_partial_wrong_curv":
        batch_prefix = "O2_"
    elif model_name == "con_wrong_curv":
        batch_prefix = "C_"
    elif model_name == "obj_con_wrong_curv":
        batch_prefix = "OC_"
    else:
        raise ValueError("Wrong model name.")

    do_batch_test_algo12(model_name, batch_prefix+"U", 0, 1.1, 0.1, default_r_0,default_r_max, default_sigma, default_xi_N)
    do_batch_test_algo12(model_name, batch_prefix+"U", 1, -1.1, -0.1, default_r_0,default_r_max, default_sigma, default_xi_N)
    do_batch_test_algo12(model_name, batch_prefix+"U", 2, 1.5, -1.5, default_r_0,default_r_max, default_sigma, default_xi_N)
    do_batch_test_algo12(model_name, batch_prefix+"U", 3, 1.5, 1.5, default_r_0,default_r_max, default_sigma, default_xi_N)
    do_batch_test_algo12(model_name, batch_prefix+"U", 4, -1.5, -1.5, default_r_0,default_r_max, default_sigma, default_xi_N)
    do_batch_test_algo12(model_name, batch_prefix+"U", 5, 3.5, -2.5, default_r_0,default_r_max, default_sigma, default_xi_N)
    do_batch_test_algo12(model_name, batch_prefix+"U", 6, 3.5, -1, default_r_0,default_r_max, default_sigma, default_xi_N)
    do_batch_test_algo12(model_name, batch_prefix+"U", 7, 3.5, 0.5, default_r_0,default_r_max, default_sigma, default_xi_N)

    do_batch_test_algo12(model_name, batch_prefix+"FP", 0, default_u1_0, default_u2_0, default_r_0,default_r_max, default_sigma, default_xi_N)
    do_batch_test_algo12(model_name, batch_prefix+"FP", 1, default_u1_0, default_u2_0, default_r_0,default_r_max, 5, default_xi_N)
    do_batch_test_algo12(model_name, batch_prefix+"FP", 2, default_u1_0, default_u2_0, default_r_0,default_r_max, 1, default_xi_N)
    do_batch_test_algo12(model_name, batch_prefix+"FP", 3, default_u1_0, default_u2_0, default_r_0,default_r_max, 0.5, default_xi_N)
    do_batch_test_algo12(model_name, batch_prefix+"FP", 4, default_u1_0, default_u2_0, default_r_0,default_r_max, 0.1, default_xi_N)

    do_batch_test_algo12(model_name, batch_prefix+"IFP", 0, -1.1, -0.1, default_r_0,default_r_max, default_sigma, default_xi_N)
    do_batch_test_algo12(model_name, batch_prefix+"IFP", 1, -1.1, -0.1, default_r_0,default_r_max, 5, default_xi_N)
    do_batch_test_algo12(model_name, batch_prefix+"IFP", 2, -1.1, -0.1, default_r_0,default_r_max, 1, default_xi_N)
    do_batch_test_algo12(model_name, batch_prefix+"IFP", 3, -1.1, -0.1, default_r_0,default_r_max, 0.5, default_xi_N)
    do_batch_test_algo12(model_name, batch_prefix+"IFP", 4, -1.1, -0.1, default_r_0,default_r_max, 0.1, default_xi_N)

    do_batch_test_algo12(model_name, batch_prefix+"R", 0, default_u1_0, default_u2_0, default_r_0, \
                         default_r_max, default_sigma, default_xi_N)
    do_batch_test_algo12(model_name, batch_prefix+"R", 1, default_u1_0, default_u2_0, default_r_0/2, \
                         default_r_max/2, default_sigma, default_xi_N)
    do_batch_test_algo12(model_name, batch_prefix+"R", 2, default_u1_0, default_u2_0, default_r_0 / 5, \
                         default_r_max / 5, default_sigma, default_xi_N)
    do_batch_test_algo12(model_name, batch_prefix+"R", 3, default_u1_0, default_u2_0, default_r_0 / 10, \
                         default_r_max / 10, default_sigma, default_xi_N)
    do_batch_test_algo12(model_name, batch_prefix+"R", 4, default_u1_0, default_u2_0, default_r_0 / 20, \
                         default_r_max / 20, default_sigma, default_xi_N)

    do_batch_test_algo12(model_name, batch_prefix+"XINM", 0, default_u1_0, default_u2_0, default_r_0, \
                         default_r_max, default_sigma, default_xi_N)
    do_batch_test_algo12(model_name, batch_prefix+"XINM", 1, default_u1_0, default_u2_0, default_r_0, \
                         default_r_max, default_sigma, 0.1)
    do_batch_test_algo12(model_name, batch_prefix+"XINM", 2, default_u1_0, default_u2_0, default_r_0, \
                         default_r_max, default_sigma, 0.9)

    do_batch_test_algo12(model_name, batch_prefix+"XISR", 0, default_u1_0, default_u2_0, default_r_0/10, \
                         default_r_max/10, default_sigma, default_xi_N)
    do_batch_test_algo12(model_name, batch_prefix+"XISR", 1, default_u1_0, default_u2_0, default_r_0/10, \
                         default_r_max/10, default_sigma, 0.1)
    do_batch_test_algo12(model_name, batch_prefix+"XISR", 2, default_u1_0, default_u2_0, default_r_0/10, \
                         default_r_max/10, default_sigma, 0.9)


def do_batch_test_algo13_main(model_name):
    default_u1_0 = 1.1
    default_u2_0 = 0.1
    default_r_0 = 1
    default_sigma = 10
    default_xi_N = 0.5
    default_r_max = 4

    if model_name == "norminal":
        batch_prefix = "N_"
    elif model_name == "obj_wrong_curv":
        batch_prefix = "O1_"
    elif model_name == "obj_partial_wrong_curv":
        batch_prefix = "O2_"
    elif model_name == "con_wrong_curv":
        batch_prefix = "C_"
    elif model_name == "obj_con_wrong_curv":
        batch_prefix = "OC_"
    else:
        raise ValueError("Wrong model name.")

    do_batch_test_algo13(model_name, batch_prefix+"U", 0, 3.5, -2.5, default_r_0,default_r_max, default_sigma, default_xi_N)
    do_batch_test_algo13(model_name, batch_prefix+"U", 1, 3.5, -1, default_r_0,default_r_max, default_sigma, default_xi_N)
    do_batch_test_algo13(model_name, batch_prefix+"U", 2, 3.5, 0.5, default_r_0,default_r_max, default_sigma, default_xi_N)
    do_batch_test_algo13(model_name, batch_prefix+"U", 3, 2, -1, default_r_0,default_r_max, default_sigma, default_xi_N)
    do_batch_test_algo13(model_name, batch_prefix+"U", 4, 0.5, -1.5, default_r_0,default_r_max, default_sigma, default_xi_N)


def do_all_batches_for_all_model():
    model_names = ["norminal", "obj_wrong_curv", "obj_partial_wrong_curv", "con_wrong_curv", \
                   "obj_con_wrong_curv"]
    for model_name in model_names:
        do_batch_test_algo12_main(model_name)
        do_batch_test_algo13_main(model_name)

def do_test():
    # ------------------------------------
    noise_filename = "noise/noise_0.txt"
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    result_filename_folder="data/"
    if not os.path.exists(result_filename_folder):
        os.makedirs(result_filename_folder)

    # ------------------------------------
    starting_point = {
        "u1": -1.1,
        "u2": -0.1,
    }

    # ------------------------------------
    perturbation_stepsize = {
        "u1": 0.01,
        "u2": 0.01,
    }
    # ------------------------------------
    print_iter_data = False
    max_iter = 25
    initial_trust_radius = 1
    sigma = 10
    xi_N = 0.5
    max_trust_radius=4
    model_name = "norminal"

    # # ------------------------------------
    print("\nTesting TR_MA")
    result_filename_header = result_filename_folder + "TR_MA_"
    algo3_TR_MA(model_name, perturbation_stepsize, starting_point, initial_trust_radius,max_trust_radius,\
               noise_filename, solver_executable, print_iter_data, max_iter, \
               result_filename_header)

    # ------------------------------------
    print("\nTesting Penalty_TR_MA")
    result_filename_header = result_filename_folder + "Penalty_TR_MA_"
    algo2_TR_MA(model_name, perturbation_stepsize, starting_point, sigma, initial_trust_radius,max_trust_radius,\
                noise_filename, solver_executable, print_iter_data, max_iter, \
                result_filename_header)

    # ------------------------------------
    print("\nTesting CompoStep_TR_MA")
    result_filename_header = result_filename_folder + "CompoStep_TR_MA_"
    algo1_TR_MA(model_name, perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,\
                xi_N,\
                noise_filename, solver_executable, print_iter_data, max_iter, \
                result_filename_header)


if __name__ == "__main__":
    # generate_noise_file()
    do_test()
    # do_batch_test_algo12_main()
    # do_batch_test_algo13_main()