from rtolib.model.wo_reactor_con import RTO_Plant_WO_reactor, RTO_Mismatched_WO_reactor,\
    RTO_Mismatched_WO_reactor_QC, default_WOR_description, ZeroModel_WO_reactor,\
    QuadraticModel_WO_reactor
from rtolib.core import NoiseGenerator, ModifierType,\
                    SimpleFiniteDiffPerturbation
from rtolib.core.algo import ModifierAdaptationCompoStepTR, MACompoStepTRBackupModel
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

def fit_convex_quadratic_surrogate():
    model_simulator = PyomoSimulator(RTO_Mismatched_WO_reactor())
    model_simulator.build(default_WOR_description)
    default_options = {'max_iter': 100}
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    solver = SolverFactory('ipopt', executable=solver_executable)
    model_simulator.set_solver(solver, tee=False, default_options=default_options)

    mvs = default_WOR_description.symbol_list['MV']
    dimension = len(mvs)
    n_data = 100
    lb = []
    ub = []
    for mv in mvs:
        lb.append(default_WOR_description.bounds[mv][0])
        ub.append(default_WOR_description.bounds[mv][1])
    X = get_hypercube_sampling(dimension, n_data, lb, ub, seed=1)

    cvs=['cost', 'con']
    Y = np.zeros((n_data,len(cvs)), dtype=float)
    param_values = {}
    for k in model_simulator.pyomo_model.parameters.keys():
        param_values[k] = model_simulator.pyomo_model.default_value[k]
    for i in range(n_data):
        input_dict = {}
        for j,mv in enumerate(mvs):
            input_dict[mv] = X[i,j]
        outputs, solve_status = model_simulator.simulate(input_dict, param_values=param_values,
                                                         use_homo=False)
        for j, cv in enumerate(cvs):
            Y[i,j] = outputs[cv]

    mv_strings = ['m.Fa', 'm.Fb', 'm.Tr']
    for j, cv in enumerate(cvs):
        c, b, A, _ = fit_convex_quadratic_model(X, Y[:, j])
        print("fitting %s" % cv)
        print("c=")
        print(c)
        print("b=")
        print(b)
        print("A=")
        print(A)

        model_string = "%.6e"%c+"\\\n"
        for k1, mv1 in enumerate(mvs):
            model_string += "+%.6e" % b[k1] +"*"+ mv_strings[k1]+"\\\n"
        for k1, mv1 in enumerate(mvs):
            for k2, mv2 in enumerate(mvs):
                model_string+="+%.6e"%A[k1,k2]+"*"+mv_strings[k1]+"*"+mv_strings[k2]+"\\\n"
        print(model_string)


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

def compo_step_TR_MA_Lmodel(perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,xi_N,\
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
        model=ZeroModel_WO_reactor(),
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

def compo_step_TR_MA_Backup_Model(perturbation_stepsize, starting_point, sigma, initial_trust_radius, \
                                  max_trust_radius,xi_N,kappa_r,\
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
        "kappa_r": kappa_r,
        'feasibility_tol': global_parameter["feasibility_tol"],
        "stationarity_tol": global_parameter["stationarity_tol"],
        'adaptive_sigma': True,
        'separate_tr_management':True,
    }
    rto_algorithm.set_algorithm_option(options)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=RTO_Plant_WO_reactor(),
        model=RTO_Mismatched_WO_reactor(),
        backup_model=ZeroModel_WO_reactor(),# QuadraticModel_WO_reactor
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
    # starting_point = {
    #     "Fa": 3.8,
    #     "Fb": 9.3,
    #     "Tr": 91,
    # }

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
    sigma = 100
    xi_N = 0.5
    kappa_r=0.2
    max_trust_radius=0.4

    # ------------------------------------
    # print("\nTesting CompoStep_TR_MA")
    # result_filename_header = result_filename_folder + "CompoStep_TR_MA_"
    # compo_step_TR_MA(perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,\
    #             xi_N,\
    #             noise_filename, solver_executable, print_iter_data, max_iter, \
    #             result_filename_header)
    # ------------------------------------
    # print("\nTesting CompoStep_TR_MA with linear model")
    # result_filename_header = result_filename_folder + "CompoStep_TR_Linear_"
    # compo_step_TR_MA_Lmodel(perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,\
    #             xi_N,\
    #             noise_filename, solver_executable, print_iter_data, max_iter, \
    #             result_filename_header)
    # ------------------------------------
    print("\nTesting CompoStep_TR_MA with Backup model")
    result_filename_header = result_filename_folder + "CompoStep_TR_BkpModel_"
    compo_step_TR_MA_Backup_Model(perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius, \
                     xi_N, kappa_r,\
                     noise_filename, solver_executable, print_iter_data, max_iter, \
                     result_filename_header)


def plot_profile():
    max_iter=30
    global_marker_size = 2
    linewidth=1
    compo_step_plant_data = pandas.read_csv("data/wo_single/CompoStep_TR_MA_plant_data.txt", \
                                      index_col=0, header=0, sep='\t')
    compo_step_plant_data2 = pandas.read_csv("data/wo_single/CompoStep_TR_Linear_plant_data.txt", \
                                            index_col=0, header=0, sep='\t')
    backup_tr_plant_data = pandas.read_csv("data/wo_single/CompoStep_TR_BkpModel_plant_data.txt", \
                                   index_col=0, header=0, sep='\t')

    compo_step_model_data = pandas.read_csv("data/wo_single/CompoStep_TR_MA_model_data.txt", \
                                      index_col=0, header=0, sep='\t')
    compo_step_model_data2 = pandas.read_csv("data/wo_single/CompoStep_TR_Linear_model_data.txt", \
                                            index_col=0, header=0, sep='\t')
    backup_tr_model_data = pandas.read_csv("data/wo_single/CompoStep_TR_BkpModel_model_data.txt", \
                                   index_col=0, header=0, sep='\t')

    compo_step_input_data = pandas.read_csv("data/wo_single/CompoStep_TR_MA_input_data.txt", \
                                            index_col=0, header=0, sep='\t')
    compo_step_input_data2 = pandas.read_csv("data/wo_single/CompoStep_TR_Linear_input_data.txt", \
                                            index_col=0, header=0, sep='\t')
    backup_tr_input_data = pandas.read_csv("data/wo_single/CompoStep_TR_BkpModel_input_data.txt", \
                                         index_col=0, header=0, sep='\t')

    fig = plt.figure(figsize=(6,12))
    plt.subplot(711)
    optimal = -210.33 * np.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1,max_iter+1), compo_step_plant_data.loc[1:max_iter, 'cost'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), compo_step_plant_data2.loc[1:max_iter, 'cost'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), backup_tr_plant_data.loc[1:max_iter, 'cost'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel("plant cost")
    plt.subplot(712)
    optimal = 0 * np.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_plant_data.loc[1:max_iter, 'con'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), compo_step_plant_data2.loc[1:max_iter, 'con'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), backup_tr_plant_data.loc[1:max_iter, 'con'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel("plant constraints")
    plt.subplot(713)
    optimal = 3.887 * np.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_input_data.loc[1:max_iter, 'Fa'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), compo_step_input_data2.loc[1:max_iter, 'Fa'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), backup_tr_input_data.loc[1:max_iter, 'Fa'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"Fa")
    plt.subplot(714)
    optimal = 9.369 * np.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_input_data.loc[1:max_iter, 'Fb'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), compo_step_input_data2.loc[1:max_iter, 'Fb'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), backup_tr_input_data.loc[1:max_iter, 'Fb'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"Fb")
    plt.subplot(715)
    optimal = 91.2 * np.ones(max_iter + 1)
    plt.plot(range(max_iter + 1), optimal, linewidth=linewidth, label='Optimal', color='gray',
             linestyle='--')
    plt.plot(range(1, max_iter + 1), compo_step_input_data.loc[1:max_iter, 'Tr'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), compo_step_input_data2.loc[1:max_iter, 'Tr'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), backup_tr_input_data.loc[1:max_iter, 'Tr'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"Tr")
    plt.subplot(716)
    plt.plot(range(1, max_iter + 1), compo_step_model_data.loc[1:max_iter, 'tr'], \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), compo_step_model_data2.loc[1:max_iter, 'tr'], \
             marker='o', c='blue', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), backup_tr_model_data.loc[1:max_iter, 'tr'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.plot(range(1, max_iter + 1), backup_tr_model_data.loc[1:max_iter, 'tr_b'], \
             marker='o', c='orange', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel(r"tr")
    plt.subplot(717)
    plt.plot(range(1, max_iter + 1), np.ones(shape=(max_iter,)), \
             marker='o', c='black', markersize=global_marker_size, linewidth=linewidth*2)
    plt.plot(range(1, max_iter + 1), backup_tr_model_data.loc[1:max_iter, 'selected_model'], \
             marker='o', c='red', markersize=global_marker_size, linewidth=linewidth)
    plt.ylabel("selected_model")
    plt.xlabel("#iteration")

    plt.savefig("pic/single_wo/profile", dpi=600)
    plt.close()

if __name__ == "__main__":
    # generate_noise_file()
    do_test()
    plot_profile()
    # fit_convex_quadratic_surrogate()