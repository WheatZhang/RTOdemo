from rtolib.model.asu17000 import default_ASU_description, ASU_Plant, ASU_Quadratic_Model,\
                ASU_Model, ASU_QuadObjLinearCon_Model
from rtolib.core import NoiseGenerator, ModifierType,\
                    SimpleFiniteDiffPerturbation
from rtolib.core.algo import ModifierAdaptationCompoStepTR, MACompoStepTRBackupModel
import copy
from rtolib.util.misc import save_iteration_data_in_dict
import os
import pandas
import matplotlib.pyplot as plt
import numpy as np
from rtolib.core.solve import PyomoOptimizer, PyomoSimulator
from pyomo.environ import SolverFactory, value


global_parameter={
        "eta1": 0.01,
        "eta2": 0.9,
        "gamma1": 0.5,
        "gamma2": 1,
        "gamma3": 2,
        "feasibility_tol": 1e-2,
        "stationarity_tol": 1e-2,
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
        'adaptive_sigma': False,
    }
    rto_algorithm.set_algorithm_option(options)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=ASU_Plant(),
        model=ASU_QuadObjLinearCon_Model(), #ASU_Quadratic_Model(),
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
    # solver = rto_algorithm.model_optimizer.solver
    # solver.options['tol'] = 1e-8  # this could be smaller if necessary
    # solver.options['bound_push'] = 1e-10  # this is also crucial, for both convergence and speed
    # solver.options['max_iter'] = 1000
    # solver.options['linear_solver'] = 'ma57'
    rto_algorithm.model_optimizer.solver = SolverFactory('gurobi')

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


def do_plant_optimization():
    pyomo_optimizer = PyomoOptimizer(ASU_Plant())
    pyomo_simulator = PyomoSimulator(ASU_Plant())
    problem_description = copy.deepcopy(default_ASU_description)

    pyomo_optimizer.build(problem_description)
    pyomo_simulator.build(problem_description)
    solver_executable = r"F:\Research\ModularTask\wys_fit\WholeOpt\ipopt38.exe"

    solver = SolverFactory('ipopt', executable=solver_executable)
    default_options = {'max_iter': 1000,
                       "tol": 1e-8,
                       'bound_push':1e-10,
                       'linear_solver':'ma57'}
    pyomo_optimizer.set_solver(solver, tee=True, default_options=default_options)
    pyomo_simulator.set_solver(solver, tee=True, default_options=default_options)
    homotopy_var = problem_description.symbol_list['MV']
    pyomo_simulator.homotopy_var = homotopy_var

    optimized_input, solve_status = pyomo_optimizer.optimize(input_values={},
                                                                 param_values={},
                                                                 use_homo=True)
    outputs, solve_status = pyomo_simulator.simulate(optimized_input, param_values=None,
                                                         use_homo=True)
    print(optimized_input)
    print(outputs)
    '''{'FeedMFlow': 4030.7615429745706, 'TAFeedFrac': 0.18402973846900386, 
    'MAFeedFrac': 0.517016433016136, '40_LIN': 1381.0459258448782, 
    '52_WN': 1750.8448614448741, '47_ARC': 647.4296886102883, 
    'ASC_RefluxRatio': 27.96633374135291, 'LOX_Frac': 0.08354167030750419}
    {'cost': 7378.742019977917, 'OxygenPrdtPurityCon': -1.621611782298693e-06, 
    'OxygenPrdtImpurityCon': -7.999899600381685e-06, 
    'NitrogenPrdtPurityCon': 1.0108766712590977e-08, 
    'NitrogenPrdtImpurityCon': -2.3768411263457776e-06, 
    'ArgonPrdtImpurityCon': 0.00015310183054054896, 
    'CrudeArgonImpurityCon': -1.6744272376809002e-06, 
    'LOX_Con': 0.015713949209839484, 
    'GAN_Con': -0.01702088294518944, 
    'GAR_Con': -0.9135942036078895, 
    'MainCoolerTempDiffCon': -0.32969827003697105, 
    'HeatLPCTempDiffCon': -2.715844325550819, 'HPA_GOX_LB_Con': -401.66214375533605, 
    'HPA_GOX_UB_Con': -7.519419068557909e-05, 'Feed_LB_Con': -1353.817074751576, 
    'Feed_UB_Con': -2661.599525248423, 'HPA_LB_Con': -462.93658165439, 
    'HPA_UB_Con': -650.13841834561, 'TA_LB_Con': -0.011898809495505702,
     'TA_UB_Con': -1112.6581011905043, 'OverallHeatCon': -1726547.801505413}
     '''

def do_model_optimization():
    pyomo_optimizer = PyomoOptimizer(ASU_Model())
    pyomo_simulator = PyomoSimulator(ASU_Model())
    problem_description = copy.deepcopy(default_ASU_description)

    pyomo_optimizer.build(problem_description)
    pyomo_simulator.build(problem_description)
    solver_executable = r"F:\Research\ModularTask\wys_fit\WholeOpt\ipopt38.exe"

    solver = SolverFactory('ipopt', executable=solver_executable)
    default_options = {'max_iter': 1000,
                       "tol": 1e-8,
                       'bound_push':1e-10,
                       'linear_solver':'ma57'}
    pyomo_optimizer.set_solver(solver, tee=True, default_options=default_options)
    pyomo_simulator.set_solver(solver, tee=True, default_options=default_options)
    homotopy_var = problem_description.symbol_list['MV']
    pyomo_simulator.homotopy_var = homotopy_var

    optimized_input, solve_status = pyomo_optimizer.optimize(input_values={},
                                                                 param_values={},
                                                                 use_homo=True)
    outputs, solve_status = pyomo_simulator.simulate(optimized_input, param_values=None,
                                                         use_homo=True)
    print(optimized_input)
    print(outputs)
    '''
    {'FeedMFlow': 4086.00373608555, 'TAFeedFrac': 0.1815416824109121, 
    'MAFeedFrac': 0.5235463025900193, '40_LIN': 1402.9400273570284, 
    '52_WN': 1793.04123189593, '47_ARC': 633.5245177275037, 
    'ASC_RefluxRatio': 28.14249703352842, 'LOX_Frac': 0.09843645583068322}
    {'cost': 7231.003148581344, 'OxygenPrdtPurityCon': 
    -1.4013088117659223e-06, 'OxygenPrdtImpurityCon': -7.99983808385816e-06,
     'NitrogenPrdtPurityCon': 8.370325210727003e-09, 
     'NitrogenPrdtImpurityCon': -4.3400103902222425e-06, 
     'ArgonPrdtImpurityCon': 0.00015373475567276437, 
     'CrudeArgonImpurityCon': -1.2774193726565507e-06, 
     'LOX_Con': 0.014880154359616427, 'GAN_Con': -0.01640691216562118,
      'GAR_Con': -0.27813530041179035, 'MainCoolerTempDiffCon': -0.3297815476583139,
       'HeatLPCTempDiffCon': -2.677131176731507, 'HPA_GOX_LB_Con': -401.6618332012389, 
       'HPA_GOX_UB_Con': -0.0008270425503269507, 'Feed_LB_Con': -1409.059273510813, 
       'Feed_UB_Con': -2606.357326489186, 'HPA_LB_Con': -462.9371536888177, 
       'HPA_UB_Con': -650.1378463111823, 'TA_LB_Con': -0.01145476737417539, 
       'TA_UB_Con': -1112.6585452326258, 'OverallHeatCon': -1914750.528234668}
    '''


def do_quadratic_model_optimization():
    pyomo_optimizer = PyomoOptimizer(ASU_Quadratic_Model())
    pyomo_simulator = PyomoSimulator(ASU_Quadratic_Model())
    problem_description = copy.deepcopy(default_ASU_description)

    pyomo_optimizer.build(problem_description)
    pyomo_simulator.build(problem_description)
    solver_executable = r"F:\Research\ModularTask\wys_fit\WholeOpt\ipopt38.exe"

    solver = SolverFactory('ipopt', executable=solver_executable)
    default_options = {'max_iter': 1000,
                       "tol": 1e-8,
                       'bound_push':1e-10,
                       'linear_solver':'ma57'}
    pyomo_optimizer.set_solver(solver, tee=True, default_options=default_options)
    pyomo_simulator.set_solver(solver, tee=True, default_options=default_options)
    homotopy_var = problem_description.symbol_list['MV']
    pyomo_simulator.homotopy_var = homotopy_var

    optimized_input, solve_status = pyomo_optimizer.optimize(input_values={},
                                                                 param_values={},
                                                                 use_homo=False)
    outputs, solve_status = pyomo_simulator.simulate(optimized_input, param_values=None,
                                                         use_homo=False)
    print(optimized_input)
    print(outputs)

def do_test():
    # ------------------------------------
    noise_filename = "noise/noise_0_asu.txt"
    solver_executable = r"F:\Research\ModularTask\wys_fit\WholeOpt\ipopt38.exe"
    result_filename_folder="data/asu/"
    if not os.path.exists(result_filename_folder):
        os.makedirs(result_filename_folder)

    # ------------------------------------
    # optimal for LOX=17100
    starting_point = {'FeedMFlow': 4052.7912869166826,
                      'TAFeedFrac': 0.18302940863970732,
                      'MAFeedFrac': 0.5178927921861316,
                      '40_LIN': 1389.398643236176,
                      '52_WN': 1768.023215065772,
                      '47_ARC': 651.5039098061236,
                      'ASC_RefluxRatio': 27.983588697041924,
                      'LOX_Frac': 0.08337413512382041,
                      }
    # a random starting point
    # starting_point = {
    #     'FeedMFlow':4100,
    #     'TAFeedFrac':0.19,
    #     'MAFeedFrac':0.5,
    #     '40_LIN':1450,
    #     '52_WN':1800,
    #     '47_ARC':600,
    #     'ASC_RefluxRatio':28,
    #     'LOX_Frac':0.09,
    # }
    # ------------------------------------
    perturbation_stepsize = {
        'FeedMFlow':5,
        'TAFeedFrac':0.01,
        'MAFeedFrac':0.01,
        '40_LIN':2,
        '52_WN':2,
        '47_ARC':1,
        'ASC_RefluxRatio':0.5,
        'LOX_Frac':0.01,
    }
    # ------------------------------------
    print_iter_data = False
    max_iter = 20
    initial_trust_radius = 0.2
    sigma = 100 #100
    xi_N = 0.5
    max_trust_radius=0.8
    kappa_b = 0.2

    # ------------------------------------
    print("\nTesting CompoStep_TR")
    result_filename_header = result_filename_folder + "TR_quadratic_"
    compo_step_TR_quadratic_model(perturbation_stepsize, starting_point, sigma, initial_trust_radius, max_trust_radius,\
            xi_N,\
            noise_filename, solver_executable, print_iter_data, max_iter, \
            result_filename_header)
    # ------------------------------------
    # print("\nTesting Backuped CompoStep_TR")
    # result_filename_header = result_filename_folder + "TR_odm_backuped_"
    # compo_step_TR_MA_Backup_Model(perturbation_stepsize, starting_point, sigma, initial_trust_radius, \
    #                               max_trust_radius, xi_N, kappa_b, \
    #                               noise_filename, solver_executable, print_iter_data, max_iter, \
    #                               result_filename_header)


if __name__ == "__main__":
    # generate_noise_file()
    do_test()
    # do_plant_optimization()
    # do_model_optimization()
    # do_quadratic_model_optimization()