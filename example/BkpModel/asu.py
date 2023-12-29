from rtolib.model.asu17000 import default_ASU_description, ASU_Plant, ASU_Quadratic_Model,\
                ASU_Model
from rtolib.core import NoiseGenerator, ModifierType,\
                    SimpleFiniteDiffPerturbation
from rtolib.core.algo import ModifierAdaptationCompoStepTR, MACompoStepTRBackupModel
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
def generate_noise_file():
    noise_level = {
        "cost":0,
        "OxygenPrdtPurityCon":0,
        "OxygenPrdtImpurityCon":0,
        "NitrogenPrdtPurityCon":0,
        "NitrogenPrdtImpurityCon":0,
        "ArgonPrdtImpurityCon":0,
        "CrudeArgonImpurityCon":0,
        "LOX_Con":0,
        "GAN_Con":0,
        "GAR_Con":0,
        "MainCoolerTempDiffCon":0,
        "HeatLPCTempDiffCon":0,
        "HPA_GOX_LB_Con":0,
        "HPA_GOX_UB_Con":0,
        "Feed_LB_Con":0,
        "Feed_UB_Con":0,
        "HPA_LB_Con":0,
        "HPA_UB_Con":0,
        "TA_LB_Con":0,
        "TA_UB_Con":0,
        "OverallHeatCon":0,
    }
    max_iter = 200
    noise_generator = NoiseGenerator(noise_level)
    noise_generator.generate_noise(max_iter, 30)
    noise_generator.save_noise("noise/noise_0_asu.txt")

def compo_step_TR_quadratic_model(perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,xi_N,\
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
    problem_description = copy.deepcopy(default_ASU_description)

    rto_algorithm = ModifierAdaptationCompoStepTR()
    options = {
        "homotopy_simulation": True,
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
        plant=ASU_Plant(),
        model=ASU_Quadratic_Model(),
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        solver_executable=solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        black_box_model=None,
        )

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    solver = rto_algorithm.plant_simulator.solver
    solver.options['tol'] = 1e-8  # this could be smaller if necessary
    solver.options['bound_push'] = 1e-10  # this is also crucial, for both convergence and speed
    solver.options['max_iter'] = 1000
    solver.options['linear_solver'] = 'ma57'
    solver = rto_algorithm.model_simulator.solver
    solver.options['tol'] = 1e-8  # this could be smaller if necessary
    solver.options['bound_push'] = 1e-10  # this is also crucial, for both convergence and speed
    solver.options['max_iter'] = 1000
    solver.options['linear_solver'] = 'ma57'
    solver = rto_algorithm.model_optimizer.solver
    solver.options['tol'] = 1e-8  # this could be smaller if necessary
    solver.options['bound_push'] = 1e-10  # this is also crucial, for both convergence and speed
    solver.options['max_iter'] = 1000
    solver.options['linear_solver'] = 'ma57'

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

def compo_step_TR_MA_Backup_Model(perturbation_stepsize, starting_point, sigma, initial_trust_radius, \
                                  max_trust_radius,xi_N,kappa_b,\
               noise_filename, solver_executable, print_iter_data, max_iter,\
               result_filename_header):
    problem_description = copy.deepcopy(default_ASU_description)

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
        'separate_tr_management': True,
        "skip_backup": False,
        "use_premature_solution": False,
    }
    rto_algorithm.set_algorithm_option(options)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=ASU_Plant(),
        model=ASU_Model(),
        backup_model=ASU_Quadratic_Model(),
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        solver_executable=solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        )

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    solver = rto_algorithm.plant_simulator.solver
    solver.options['tol'] = 1e-8  # this could be smaller if necessary
    solver.options['bound_push'] = 1e-10  # this is also crucial, for both convergence and speed
    solver.options['max_iter'] = 1000
    solver.options['linear_solver'] = 'ma57'
    solver = rto_algorithm.model_simulator.solver
    solver.options['tol'] = 1e-8  # this could be smaller if necessary
    solver.options['bound_push'] = 1e-10  # this is also crucial, for both convergence and speed
    solver.options['max_iter'] = 1000
    solver.options['linear_solver'] = 'ma57'
    solver = rto_algorithm.model_optimizer.solver
    solver.options['tol'] = 1e-8  # this could be smaller if necessary
    solver.options['bound_push'] = 1e-10  # this is also crucial, for both convergence and speed
    solver.options['max_iter'] = 1000
    solver.options['linear_solver'] = 'ma57'
    solver = rto_algorithm.backup_model_simulator.solver
    solver.options['max_iter'] = 1000
    solver.options['linear_solver'] = 'ma57'
    solver = rto_algorithm.backup_model_optimizer.solver
    solver.options['max_iter'] = 1000
    solver.options['linear_solver'] = 'ma57'

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
    noise_filename = "noise/noise_0_asu.txt"
    solver_executable = r"F:\Research\ModularTask\wys_fit\WholeOpt\ipopt38.exe"
    result_filename_folder="data/asu/"
    if not os.path.exists(result_filename_folder):
        os.makedirs(result_filename_folder)

    # ------------------------------------
    starting_point = {
        'FeedMFlow':4100,
        'TAFeedFrac':0.19,
        'MAFeedFrac':0.5,
        '40_LIN':1450,
        '52_WN':1800,
        '47_ARC':600,
        'ASC_RefluxRatio':28,
        'LOX_Frac':0.09,
    }
    # ------------------------------------
    perturbation_stepsize = {
        'FeedMFlow':2,
        'TAFeedFrac':0.001,
        'MAFeedFrac':0.001,
        '40_LIN':2,
        '52_WN':2,
        '47_ARC':1,
        'ASC_RefluxRatio':0.2,
        'LOX_Frac':0.001,
    }
    # ------------------------------------
    print_iter_data = False
    max_iter = 20
    initial_trust_radius = 0.5
    sigma = 10 #100
    xi_N = 0.5
    max_trust_radius=2
    kappa_b = 0.2

    # ------------------------------------
    print("\nTesting CompoStep_TR")
    result_filename_header = result_filename_folder + "TR_quadratic_"
    compo_step_TR_quadratic_model(perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,\
            xi_N,\
            noise_filename, solver_executable, print_iter_data, max_iter, \
            result_filename_header)
    # ------------------------------------
    print("\nTesting Backuped CompoStep_TR")
    result_filename_header = result_filename_folder + "TR_odm_backuped_"
    compo_step_TR_MA_Backup_Model(perturbation_stepsize, starting_point, sigma, initial_trust_radius, \
                                  max_trust_radius, xi_N, kappa_b, \
                                  noise_filename, solver_executable, print_iter_data, max_iter, \
                                  result_filename_header)


if __name__ == "__main__":
    # generate_noise_file()
    do_test()