from rtolib.model.wo_reactor_con import RTO_Plant_WO_reactor, RTO_Mismatched_WO_reactor,\
    RTO_Mismatched_WO_reactor_QC, default_WOR_description, ZeroModel_WO_reactor,\
    QuadraticModel_WO_reactor
from rtolib.core import NoiseGenerator, ModifierType,\
                    SimpleFiniteDiffPerturbation
from rtolib.core.algo import ModifierAdaptationCompoStepTR, MACompoStepTRBackupModel,\
                ModifierAdaptation
import copy
from rtolib.util.misc import save_iteration_data_in_dict
import os
import pandas
import matplotlib.pyplot as plt
import numpy as np
from rtolib.util.cvx_quadr import fit_convex_quadratic_model, get_hypercube_sampling
from rtolib.core.solve import PyomoSimulator
from pyomo.environ import SolverFactory

global_parameter={
        "eta1": 0.01,
        "eta2": 0.9,
        "gamma1": 0.5,
        "gamma2": 1,
        "gamma3": 2,
        "feasibility_tol": 1e-4,
        "stationarity_tol": 1e-4,
}
pic_constant = 0.39370
font_factor = np.sqrt(1/pic_constant)
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
    noise_generator.save_noise("noise/noise_0_wo_single.txt")

def original_MA(model_name, perturbation_stepsize, starting_point, \
               noise_filename, solver_executable, print_iter_data, max_iter,\
               result_filename_header):
    if model_name == "rigorous_nonlinear":
        model = RTO_Mismatched_WO_reactor()
    elif model_name == "quadratic":
        model = QuadraticModel_WO_reactor()
    elif model_name == "linear":
        model = ZeroModel_WO_reactor()
    else:
        raise ValueError("unknown model name %s"%model_name)

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

def compo_step_TR_MA(model_name, perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,xi_N,\
               noise_filename, solver_executable, print_iter_data, max_iter,\
               result_filename_header):
    if model_name == "rigorous_nonlinear":
        model = RTO_Mismatched_WO_reactor()
    elif model_name == "quadratic":
        model = QuadraticModel_WO_reactor()
    elif model_name == "linear":
        model = ZeroModel_WO_reactor()
    else:
        raise ValueError("unknown model name %s"%model_name)

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
        plant=RTO_Plant_WO_reactor(),
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

def compo_step_TR_MA_Backup_Model(backup_model_name, perturbation_stepsize, starting_point, sigma, initial_trust_radius, \
                                  max_trust_radius,xi_N,kappa_b,\
               noise_filename, solver_executable, print_iter_data, max_iter,\
               result_filename_header, separate_tr_management, skip_backup=False):
    if backup_model_name == "rigorous_nonlinear":
        backup_model = RTO_Mismatched_WO_reactor()
    elif backup_model_name == "quadratic":
        backup_model = QuadraticModel_WO_reactor()
    elif backup_model_name == "linear":
        backup_model = ZeroModel_WO_reactor()
    else:
        raise ValueError("unknown model name %s"%model_name)
    problem_description = copy.deepcopy(default_WOR_description)

    rto_algorithm = MACompoStepTRBackupModel()
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
        "kappa_b": kappa_b,
        'feasibility_tol': global_parameter["feasibility_tol"],
        "stationarity_tol": global_parameter["stationarity_tol"],
        'adaptive_sigma': True,
        'separate_tr_management': separate_tr_management,
        "skip_backup": skip_backup,
    }
    rto_algorithm.set_algorithm_option(options)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=RTO_Plant_WO_reactor(),
        model=RTO_Mismatched_WO_reactor(),
        backup_model=backup_model,
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        solver_executable=solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        )

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    solver = rto_algorithm.model_optimizer.solver
    solver.options["max_iter"] = 5

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
    noise_filename = "noise/noise_0_wo_single.txt"
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    result_filename_folder="data/wo_single/"
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
        "Fa": 0.00001,
        "Fb": 0.00001,
        "Tr": 0.00001,
    }
    # ------------------------------------
    print_iter_data = False
    max_iter = 30
    initial_trust_radius = 0.1
    sigma = 10 #100
    xi_N = 0.5
    kappa_b=0.2
    max_trust_radius=0.4
    model_types = ['rigorous_nonlinear', 'quadratic', 'linear']
    suffix = {
        'rigorous_nonlinear':"ODM",
        'quadratic':"Q",
        'linear':"L",
    }

    # ------------------------------------
    # print("\nTesting Original_MA")
    # for model in model_types:
    #     result_filename_header = result_filename_folder + "MA_"+suffix[model]+"_"
    #     original_MA(model, perturbation_stepsize, starting_point, \
    #                 noise_filename, solver_executable, print_iter_data, max_iter, \
    #                 result_filename_header)
    # ------------------------------------
    # print("\nTesting CompoStep_TR")
    # for model in model_types:
    #     result_filename_header = result_filename_folder + "TR_"+suffix[model]+"_"
    #     compo_step_TR_MA(model, perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,\
    #             xi_N,\
    #             noise_filename, solver_executable, print_iter_data, max_iter, \
    #             result_filename_header)
    # ------------------------------------
    # print("\nTesting CompoStep_TR with Backup model")
    for model in ['quadratic', 'linear']:
        result_filename_header = result_filename_folder + "BackupSepDelta_"+suffix[model]+"_"
        compo_step_TR_MA_Backup_Model(model, perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius, \
                         xi_N, kappa_b,\
                         noise_filename, solver_executable, print_iter_data, max_iter, \
                         result_filename_header, separate_tr_management=True)

        result_filename_header = result_filename_folder + "BackupIndDelta_" + suffix[model] + "_"
        compo_step_TR_MA_Backup_Model(model, perturbation_stepsize, starting_point, sigma, initial_trust_radius,
                                      max_trust_radius, \
                                      xi_N, kappa_b, \
                                      noise_filename, solver_executable, print_iter_data, max_iter, \
                                      result_filename_header, separate_tr_management=False)
    # result_filename_header = result_filename_folder + "NoBackup"
    # compo_step_TR_MA_Backup_Model('quadratic', perturbation_stepsize, starting_point, sigma, initial_trust_radius,
    #                               max_trust_radius, \
    #                          xi_N, kappa_b,\
    #                      noise_filename, solver_executable, print_iter_data, max_iter, \
    #                      result_filename_header, separate_tr_management=False,\
    #                               skip_backup=True)

def plot_profile_MA():
    max_iter=30
    global_marker_size = 2
    linewidth=1
    ma_plant_data1 = pandas.read_csv("data/wo_single/MA_ODM_plant_data.txt", \
                                      index_col=0, header=0, sep='\t')
    ma_plant_data2 = pandas.read_csv("data/wo_single/MA_Q_plant_data.txt", \
                                            index_col=0, header=0, sep='\t')
    ma_plant_data3 = pandas.read_csv("data/wo_single/MA_L_plant_data.txt", \
                                   index_col=0, header=0, sep='\t')

    ma_model_data1 = pandas.read_csv("data/wo_single/MA_ODM_model_data.txt", \
                                      index_col=0, header=0, sep='\t')
    ma_model_data2 = pandas.read_csv("data/wo_single/MA_Q_model_data.txt", \
                                            index_col=0, header=0, sep='\t')
    ma_model_data3 = pandas.read_csv("data/wo_single/MA_L_model_data.txt", \
                                   index_col=0, header=0, sep='\t')

    ma_input_data1 = pandas.read_csv("data/wo_single/MA_ODM_input_data.txt", \
                                            index_col=0, header=0, sep='\t')
    ma_input_data2 = pandas.read_csv("data/wo_single/MA_Q_input_data.txt", \
                                            index_col=0, header=0, sep='\t')
    ma_input_data3 = pandas.read_csv("data/wo_single/MA_L_input_data.txt", \
                                         index_col=0, header=0, sep='\t')

    fig = plt.figure(figsize=(6,10))
    plt.subplot(511)
    optimal = -210.33 * np.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1,max_iter+1), ma_plant_data1.loc[1:max_iter, 'cost'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2, label='Detailed')
    plt.plot(range(1, max_iter + 1), ma_plant_data2.loc[1:max_iter, 'cost'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth, label='Quadratic')
    plt.plot(range(1, max_iter + 1), ma_plant_data3.loc[1:max_iter, 'cost'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth/2, label='Zero')
    plt.ylabel("plant cost")
    plt.legend(prop=font_legend)
    plt.subplot(512)
    optimal = 0 * np.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), ma_plant_data1.loc[1:max_iter, 'con'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), ma_plant_data2.loc[1:max_iter, 'con'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), ma_plant_data3.loc[1:max_iter, 'con'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth/2)
    plt.ylabel("plant constraints")
    plt.subplot(513)
    optimal = 3.887 * np.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), ma_input_data1.loc[1:max_iter, 'Fa'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), ma_input_data2.loc[1:max_iter, 'Fa'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), ma_input_data3.loc[1:max_iter, 'Fa'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth/2)
    plt.ylabel(r"$F_A$")
    plt.subplot(514)
    optimal = 9.369 * np.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), ma_input_data1.loc[1:max_iter, 'Fb'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), ma_input_data2.loc[1:max_iter, 'Fb'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), ma_input_data3.loc[1:max_iter, 'Fb'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth/2)
    plt.ylabel(r"$F_B$")
    plt.subplot(515)
    optimal = 91.2 * np.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), ma_input_data1.loc[1:max_iter, 'Tr'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), ma_input_data2.loc[1:max_iter, 'Tr'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), ma_input_data3.loc[1:max_iter, 'Tr'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth/2)
    plt.ylabel(r"$T_R$")
    plt.xlabel("#iteration")

    plt.savefig("pic/single_wo/profile_MA", dpi=600)
    plt.close()

def plot_profile_TR():
    max_iter=30
    global_marker_size = 2
    linewidth=1
    compo_step_plant_data1 = pandas.read_csv("data/wo_single/TR_ODM_plant_data.txt", \
                                      index_col=0, header=0, sep='\t')
    compo_step_plant_data2 = pandas.read_csv("data/wo_single/TR_Q_plant_data.txt", \
                                            index_col=0, header=0, sep='\t')
    compo_step_plant_data3 = pandas.read_csv("data/wo_single/TR_L_plant_data.txt", \
                                   index_col=0, header=0, sep='\t')
    compo_step_plant_data4 = pandas.read_csv("data/wo_single/NoBackupplant_data.txt", \
                                             index_col=0, header=0, sep='\t')

    compo_step_model_data1 = pandas.read_csv("data/wo_single/TR_ODM_model_data.txt", \
                                      index_col=0, header=0, sep='\t')
    compo_step_model_data2 = pandas.read_csv("data/wo_single/TR_Q_model_data.txt", \
                                            index_col=0, header=0, sep='\t')
    compo_step_model_data3 = pandas.read_csv("data/wo_single/TR_L_model_data.txt", \
                                   index_col=0, header=0, sep='\t')
    compo_step_model_data4 = pandas.read_csv("data/wo_single/NoBackupmodel_data.txt", \
                                             index_col=0, header=0, sep='\t')

    compo_step_input_data1 = pandas.read_csv("data/wo_single/TR_ODM_input_data.txt", \
                                            index_col=0, header=0, sep='\t')
    compo_step_input_data2 = pandas.read_csv("data/wo_single/TR_Q_input_data.txt", \
                                            index_col=0, header=0, sep='\t')
    compo_step_input_data3 = pandas.read_csv("data/wo_single/TR_L_input_data.txt", \
                                         index_col=0, header=0, sep='\t')
    compo_step_input_data4 = pandas.read_csv("data/wo_single/NoBackupinput_data.txt", \
                                             index_col=0, header=0, sep='\t')

    fig = plt.figure(figsize=(6,11))
    plt.subplot(611)
    optimal = -210.33 * np.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1,max_iter+1), compo_step_plant_data1.loc[1:max_iter, 'cost'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2, label='Detailed')
    plt.plot(range(1, max_iter + 1), compo_step_plant_data2.loc[1:max_iter, 'cost'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth, label='Quadratic')
    plt.plot(range(1, max_iter + 1), compo_step_plant_data3.loc[1:max_iter, 'cost'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth, label='Zero')
    plt.plot(range(1, max_iter + 1), compo_step_plant_data4.loc[1:max_iter, 'cost'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth, label='MaxIterExcceeded')
    plt.ylabel("plant cost")
    plt.legend(prop=font_legend)
    plt.subplot(612)
    optimal = 0 * np.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_plant_data1.loc[1:max_iter, 'con'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), compo_step_plant_data2.loc[1:max_iter, 'con'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_plant_data3.loc[1:max_iter, 'con'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_plant_data4.loc[1:max_iter, 'con'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel("plant constraints")
    plt.subplot(613)
    optimal = 3.887 * np.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_input_data1.loc[1:max_iter, 'Fa'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), compo_step_input_data2.loc[1:max_iter, 'Fa'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_input_data3.loc[1:max_iter, 'Fa'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_input_data4.loc[1:max_iter, 'Fa'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"$F_A$")
    plt.subplot(614)
    optimal = 9.369 * np.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_input_data1.loc[1:max_iter, 'Fb'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), compo_step_input_data2.loc[1:max_iter, 'Fb'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_input_data3.loc[1:max_iter, 'Fb'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_input_data4.loc[1:max_iter, 'Fb'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"$F_B$")
    plt.subplot(615)
    optimal = 91.2 * np.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_input_data1.loc[1:max_iter, 'Tr'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), compo_step_input_data2.loc[1:max_iter, 'Tr'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_input_data3.loc[1:max_iter, 'Tr'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_input_data4.loc[1:max_iter, 'Tr'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"$T_R$")
    plt.subplot(616)
    plt.plot(range(1, max_iter + 1), compo_step_model_data1.loc[1:max_iter, 'tr'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), compo_step_model_data2.loc[1:max_iter, 'tr'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_model_data3.loc[1:max_iter, 'tr'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_model_data4.loc[1:max_iter, 'tr'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"$||\Delta_k||$")
    plt.xlabel("#iteration")

    plt.savefig("pic/single_wo/profile_TR", dpi=600)
    plt.close()

def plot_profile_Backup():
    max_iter=30
    global_marker_size = 2
    linewidth=1
    compo_step_plant_data1 = pandas.read_csv("data/wo_single/BackupSepDelta_Q_plant_data.txt", \
                                      index_col=0, header=0, sep='\t')
    compo_step_plant_data2 = pandas.read_csv("data/wo_single/BackupSepDelta_L_plant_data.txt", \
                                             index_col=0, header=0, sep='\t')
    compo_step_plant_data3 = pandas.read_csv("data/wo_single/BackupIndDelta_Q_plant_data.txt", \
                                            index_col=0, header=0, sep='\t')
    compo_step_plant_data4 = pandas.read_csv("data/wo_single/BackupIndDelta_L_plant_data.txt", \
                                   index_col=0, header=0, sep='\t')

    compo_step_model_data1 = pandas.read_csv("data/wo_single/BackupSepDelta_Q_model_data.txt", \
                                             index_col=0, header=0, sep='\t')
    compo_step_model_data2 = pandas.read_csv("data/wo_single/BackupSepDelta_L_model_data.txt", \
                                             index_col=0, header=0, sep='\t')
    compo_step_model_data3 = pandas.read_csv("data/wo_single/BackupIndDelta_Q_model_data.txt", \
                                             index_col=0, header=0, sep='\t')
    compo_step_model_data4 = pandas.read_csv("data/wo_single/BackupIndDelta_L_model_data.txt", \
                                             index_col=0, header=0, sep='\t')

    compo_step_input_data1 = pandas.read_csv("data/wo_single/BackupSepDelta_Q_input_data.txt", \
                                             index_col=0, header=0, sep='\t')
    compo_step_input_data2 = pandas.read_csv("data/wo_single/BackupSepDelta_L_input_data.txt", \
                                             index_col=0, header=0, sep='\t')
    compo_step_input_data3 = pandas.read_csv("data/wo_single/BackupIndDelta_Q_input_data.txt", \
                                             index_col=0, header=0, sep='\t')
    compo_step_input_data4 = pandas.read_csv("data/wo_single/BackupIndDelta_L_input_data.txt", \
                                             index_col=0, header=0, sep='\t')

    fig = plt.figure(figsize=(6,13))
    plt.subplot(711)
    optimal = -210.33 * np.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1,max_iter+1), compo_step_plant_data1.loc[1:max_iter, 'cost'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2, label="Quadratic,individual $\Delta_k$")
    plt.plot(range(1, max_iter + 1), compo_step_plant_data2.loc[1:max_iter, 'cost'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth, label="Linear,individual $\Delta_k$")
    plt.plot(range(1, max_iter + 1), compo_step_plant_data3.loc[1:max_iter, 'cost'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth, label="Quadratic,shared $\Delta_k$")
    plt.plot(range(1, max_iter + 1), compo_step_plant_data4.loc[1:max_iter, 'cost'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth, label="Linear,shared $\Delta_k$")
    plt.ylabel("plant cost")
    plt.legend(prop=font_legend)
    plt.subplot(712)
    optimal = 0 * np.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_plant_data1.loc[1:max_iter, 'con'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), compo_step_plant_data2.loc[1:max_iter, 'con'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_plant_data3.loc[1:max_iter, 'con'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_plant_data4.loc[1:max_iter, 'con'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel("plant constraints")
    plt.subplot(713)
    optimal = 3.887 * np.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_input_data1.loc[1:max_iter, 'Fa'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), compo_step_input_data2.loc[1:max_iter, 'Fa'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_input_data3.loc[1:max_iter, 'Fa'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_input_data4.loc[1:max_iter, 'Fa'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"$F_A$")
    plt.subplot(714)
    optimal = 9.369 * np.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_input_data1.loc[1:max_iter, 'Fb'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), compo_step_input_data2.loc[1:max_iter, 'Fb'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_input_data3.loc[1:max_iter, 'Fb'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_input_data4.loc[1:max_iter, 'Fb'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"$F_B$")
    plt.subplot(715)
    optimal = 91.2 * np.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_input_data1.loc[1:max_iter, 'Tr'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), compo_step_input_data2.loc[1:max_iter, 'Tr'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_input_data3.loc[1:max_iter, 'Tr'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_input_data4.loc[1:max_iter, 'Tr'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"$T_R$")
    plt.subplot(716)
    plt.plot(range(1, max_iter + 1), compo_step_model_data1.loc[1:max_iter, 'tr'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), compo_step_model_data1.loc[1:max_iter, 'tr_b'], \
             marker='^', c='black', ls="--", markersize=global_marker_size, linewidth=linewidth * 2)
    plt.plot(range(1, max_iter + 1), compo_step_model_data2.loc[1:max_iter, 'tr'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_model_data2.loc[1:max_iter, 'tr_b'], \
             marker='^', c='blue', ls="--", markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_model_data3.loc[1:max_iter, 'tr'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_model_data4.loc[1:max_iter, 'tr'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"$||\Delta_k||$")
    plt.subplot(717)
    plt.plot(range(1, max_iter + 1), compo_step_model_data1.loc[1:max_iter, 'selected_model'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth * 2)
    plt.plot(range(1, max_iter + 1), compo_step_model_data2.loc[1:max_iter, 'selected_model'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_model_data3.loc[1:max_iter, 'selected_model'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), compo_step_model_data4.loc[1:max_iter, 'selected_model'], \
             marker='o', c='green', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel("selected model")
    plt.yticks([1,2], ["m","b"])
    plt.ylim(0.5,2.5)
    plt.xlabel("#iteration")

    plt.savefig("pic/single_wo/profile_Backup", dpi=600)
    plt.close()

if __name__ == "__main__":
    # generate_noise_file()
    do_test()
    # plot_profile_MA()
    # plot_profile_TR()
    plot_profile_Backup()