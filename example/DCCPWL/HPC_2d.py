import numpy as np
import pickle
import pandas
from rtolib.core.dc_cpwl_model import QuadraticBoostedDCCPWLFunction,\
QuadraticBoostedDCCPWLFunctionSubgrad
from rtolib.core.algo import ModifierAdaptation
from rtolib.core.algo.dccpwl_ma import DCCPWL_ModifierAdaptation
from rtolib.core.algo.subgrad import DCCPWL_ModifierAdaptationSubgrad,\
            DCCPWL_ModifierAdaptationTRPenalty
from rtolib.model.HPC import RTO_Plant_HPC, default_HPC_2d_description, RTO_Mismatched_HPC
from rtolib.core import NoiseGenerator, ModifierType,\
                    SimpleFiniteDiffPerturbation
from rtolib.core.solve import PyomoSimulator
from rtolib.core.dc_cpwl_model import QuadraticBoostedDCCPWL_RTOObject
import os
import copy
from rtolib.util.misc import save_iteration_data_in_dict
from rtolib.core.pyomo_model import SolverFactory
import matplotlib.pyplot as plt

# set specification
def specification_func(iter):
    ret = {}
    ret['GAN_Flow'] = 800
    ret['LN_Flow'] = 100
    return ret

global_parameter={
        "eta1": 0.01,
        "eta2": 0.9,
        "gamma1": 0.5,
        "gamma2": 1,
        "gamma3": 2,
        "sigma": 100,
        "max_trust_radius": 500,
        "initial_trust_radius": 100,
}

def plant_simulator_test():
    plant_simulator = PyomoSimulator(RTO_Mismatched_HPC())
    plant_simulator.build(default_HPC_2d_description)
    default_options = {'max_iter': 500,
                       "tol": 1e-10}
    nlp_solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    solver1 = SolverFactory('ipopt', executable=nlp_solver_executable)
    plant_simulator.set_solver(solver1, tee=True, default_options=default_options)
    homotopy_var = ['TA_Flow', 'MA_Flow', 'GAN_Flow']
    plant_simulator.homotopy_var = homotopy_var
    point={
        'TA_Flow': 500,
        'MA_Flow': 1600,
        # 'GAN_Flow': 800,
        #'LN_Flow': 100,
    }
    outputs, solve_status = plant_simulator.simulate(point, param_values=None,
                                            use_homo=True)
    print(outputs)
    
class HPC2_2d_validity_CPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "plant_data/HPC2_2d_validity.bin"
        self.load_from_cpwl_model_file(cpwl_file)
        self.input_variables = ['TA_Flow', 'MA_Flow']
        self.parameters = []

class HPC2_2d_cost_QCPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "qmodel/HPC_2d_cpwl_processed_output_0.bin"
        quadr_file = "qmodel/HPC_2d_quadr_processed_output_0.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow']
        self.parameters = []

class HPC2_2d_feed_con1_QCPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "qmodel/HPC_2d_cpwl_processed_output_1.bin"
        quadr_file = "qmodel/HPC_2d_quadr_processed_output_1.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow']
        self.parameters = []

class HPC2_2d_feed_con2_QCPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "qmodel/HPC_2d_cpwl_processed_output_2.bin"
        quadr_file = "qmodel/HPC_2d_quadr_processed_output_2.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow']
        self.parameters = []

class HPC2_2d_purity_con_QCPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "qmodel/HPC_2d_cpwl_processed_output_3.bin"
        quadr_file = "qmodel/HPC_2d_quadr_processed_output_3.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow']
        self.parameters = []

class HPC2_2d_drain_con_QCPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "qmodel/HPC_2d_cpwl_processed_output_4.bin"
        quadr_file = "qmodel/HPC_2d_quadr_processed_output_4.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow']
        self.parameters = []

class HPC2_2d_validity_CPWL_subgrad(QuadraticBoostedDCCPWLFunctionSubgrad):
    def __init__(self):
        cpwl_file = "plant_data/HPC2_2d_validity.bin"
        self.load_from_cpwl_model_file(cpwl_file)
        self.input_variables = ['TA_Flow', 'MA_Flow']
        self.parameters = []

class HPC2_2d_cost_QCPWL_subgrad(QuadraticBoostedDCCPWLFunctionSubgrad):
    def __init__(self):
        cpwl_file = "qmodel/HPC_2d_cpwl_processed_output_0.bin"
        quadr_file = "qmodel/HPC_2d_quadr_processed_output_0.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow']
        self.parameters = []

class HPC2_2d_feed_con1_QCPWL_subgrad(QuadraticBoostedDCCPWLFunctionSubgrad):
    def __init__(self):
        cpwl_file = "qmodel/HPC_2d_cpwl_processed_output_1.bin"
        quadr_file = "qmodel/HPC_2d_quadr_processed_output_1.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow']
        self.parameters = []

class HPC2_2d_feed_con2_QCPWL_subgrad(QuadraticBoostedDCCPWLFunctionSubgrad):
    def __init__(self):
        cpwl_file = "qmodel/HPC_2d_cpwl_processed_output_2.bin"
        quadr_file = "qmodel/HPC_2d_quadr_processed_output_2.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow']
        self.parameters = []

class HPC2_2d_purity_con_QCPWL_subgrad(QuadraticBoostedDCCPWLFunctionSubgrad):
    def __init__(self):
        cpwl_file = "qmodel/HPC_2d_cpwl_processed_output_3.bin"
        quadr_file = "qmodel/HPC_2d_quadr_processed_output_3.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow']
        self.parameters = []

class HPC2_2d_drain_con_QCPWL_subgrad(QuadraticBoostedDCCPWLFunctionSubgrad):
    def __init__(self):
        cpwl_file = "qmodel/HPC_2d_cpwl_processed_output_4.bin"
        quadr_file = "qmodel/HPC_2d_quadr_processed_output_4.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow']
        self.parameters = []

def generate_noise_file():
    noise_level = {
        "obj": 0,
        "feed_con1": 0,
        "feed_con2": 0,
        "purity_con": 0,
        "drain_con": 0,
    }
    max_iter = 200
    noise_generator = NoiseGenerator(noise_level)
    noise_generator.generate_noise(max_iter, 6)
    noise_generator.save_noise("noise/hpc-noise0.txt")

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
    problem_description = copy.deepcopy(default_HPC_2d_description)
    problem_description.symbol_list['CV'] = ('obj', "feed_con1", "feed_con2", \
                                             "purity_con", "drain_con")
    rto_algorithm = ModifierAdaptation()
    options = {
        "homotopy_simulation": True,
        "homotopy_optimization": True,
        "filtering_factor": filtering_factor,
        "noise_adding_fashion": "added",
    }
    rto_algorithm.set_algorithm_option(options)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=RTO_Plant_HPC(),
        model=RTO_Mismatched_HPC(),
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        solver_executable=solver_executable,
        spec_function=specification_func,
        modifier_type=ModifierType.RTO,
    )

    homotopy_var = ['TA_Flow', 'MA_Flow', 'GAN_Flow', 'LN_Flow']
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

def do_test_QCPWL_MA(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, nlp_solver_executable, qcqp_solver_executable,\
                     print_iter_data, max_iter,\
               result_filename_header):
    problem_description = copy.deepcopy(default_HPC_2d_description)
    problem_description.symbol_list['CV']=('obj',"feed_con1","feed_con2",\
                                           "purity_con","drain_con")

    rto_algorithm = DCCPWL_ModifierAdaptation()
    options = {
        "homotopy_simulation": True,
        "homotopy_optimization": True,
        "filtering_factor": filtering_factor,
        "noise_adding_fashion": "added",
    }
    rto_algorithm.set_algorithm_option(options)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=RTO_Mismatched_HPC(),
        model_dc_cpwl_functions={'obj':HPC2_2d_cost_QCPWL(),
                                "feed_con1":HPC2_2d_feed_con1_QCPWL(),
                                "feed_con2":HPC2_2d_feed_con2_QCPWL(),
                                "purity_con":HPC2_2d_purity_con_QCPWL(),
                                "drain_con":HPC2_2d_drain_con_QCPWL(),
                                 "validity_con": HPC2_2d_validity_CPWL_subgrad(),
                                 },
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        nlp_solver_executable=nlp_solver_executable,
        qcqp_solver_executable=qcqp_solver_executable,
        spec_function=specification_func,
        modifier_type=ModifierType.RTO,
        parameter_set=[],
        )

    homotopy_var=['TA_Flow', 'MA_Flow']
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

def do_test_QCPWL_subgrad(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, nlp_solver_executable, qcqp_solver_executable,\
                     print_iter_data, max_iter,\
               result_filename_header):
    problem_description = copy.deepcopy(default_HPC_2d_description)
    problem_description.symbol_list['CV']=('obj',"feed_con1","feed_con2",\
                                           "purity_con","drain_con")

    rto_algorithm = DCCPWL_ModifierAdaptationSubgrad()
    options = {
        "homotopy_simulation": True,
        "homotopy_optimization": True,
        "filtering_factor": filtering_factor,
        "noise_adding_fashion": "added",
    }
    rto_algorithm.set_algorithm_option(options)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=RTO_Plant_HPC(),
        model_dc_cpwl_functions={'obj':HPC2_2d_cost_QCPWL_subgrad(),
                                "feed_con1":HPC2_2d_feed_con1_QCPWL_subgrad(),
                                "feed_con2":HPC2_2d_feed_con2_QCPWL_subgrad(),
                                "purity_con":HPC2_2d_purity_con_QCPWL_subgrad(),
                                "drain_con":HPC2_2d_drain_con_QCPWL_subgrad(),
                                 "validity_con":HPC2_2d_validity_CPWL_subgrad(),
                                 },
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        nlp_solver_executable=nlp_solver_executable,
        qcqp_solver_executable=qcqp_solver_executable,
        spec_function=specification_func,
        modifier_type=ModifierType.RTO,
        parameter_set=[],
        model_mvs=['TA_Flow', 'MA_Flow'],
        )

    homotopy_var=['TA_Flow', 'MA_Flow', 'GAN_Flow', 'LN_Flow']
    rto_algorithm.plant_simulator.homotopy_var = homotopy_var

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    # draw 2d contour
    num_point = 40
    x_points = np.linspace(200, 400, num_point)
    y_points = np.linspace(1600, 2000, num_point)
    x_points, y_points = np.meshgrid(x_points, y_points)
    z_points = np.zeros((num_point, num_point))
    for i in range(num_point):
        for j in range(num_point):
            output, _ = rto_algorithm.DC_CPWL_RTO_model.simulate(
                {"TA_Flow": x_points[i, j], "MA_Flow": y_points[i, j]}, with_modifier=True)
            if (output["feed_con1"] <= 0) and (output["feed_con2"] <= 0) and \
                    (output["purity_con"] <= 0) and (output["drain_con"] <= 0):
                z_points[i, j] = output["obj"]
            else:
                z_points[i, j] = 3500
    plt.figure()
    plt.contourf(x_points, y_points, z_points, 20)
    plt.colorbar()
    plt.scatter(rto_algorithm.input_history_data[0]['TA_Flow'], \
                rto_algorithm.input_history_data[0]['MA_Flow'], s=10, c='black')
    plt.savefig("pic/HPC2_subgrad/HPC2_init.png")
    plt.close()

    for iter_no in range(max_iter):
        rto_algorithm.one_step_simulation()

        # draw 2d contour
        num_point = 40
        x_points = np.linspace(200, 400, num_point)
        y_points = np.linspace(1600, 2000, num_point)
        x_points, y_points = np.meshgrid(x_points, y_points)
        z_points = np.zeros((num_point, num_point))
        for i in range(num_point):
            for j in range(num_point):
                output, _ = rto_algorithm.DC_CPWL_RTO_model.simulate(
                    {"TA_Flow": x_points[i, j], "MA_Flow": y_points[i, j]}, with_modifier=True)
                if (output["feed_con1"] <= 0) and (output["feed_con2"] <= 0) and \
                        (output["purity_con"] <= 0) and (output["drain_con"] <= 0):
                    z_points[i, j] = output["obj"]
                else:
                    z_points[i, j] = 3500
        plt.figure()
        plt.contourf(x_points, y_points, z_points, 20)
        plt.colorbar()
        plt.scatter(rto_algorithm.input_history_data[iter_no+1]['TA_Flow'], \
                    rto_algorithm.input_history_data[iter_no+1]['MA_Flow'], s=10, c='black')
        plt.savefig("pic/HPC2_subgrad/HPC2_iter" + str(iter_no) + ".png")
        plt.close()

        if print_iter_data:
            print(rto_algorithm.model_history_data)
            print(rto_algorithm.plant_history_data)
            print(rto_algorithm.input_history_data)

    # save data
    save_iteration_data_in_dict(rto_algorithm.model_history_data, result_filename_header + "model_data.txt")
    save_iteration_data_in_dict(rto_algorithm.plant_history_data, result_filename_header + "plant_data.txt")
    save_iteration_data_in_dict(rto_algorithm.input_history_data, result_filename_header + "input_data.txt")


def do_test_QCPWL_subgrad_PenaltyTR(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, nlp_solver_executable, qcqp_solver_executable,\
                     print_iter_data, max_iter,\
               result_filename_header):
    problem_description = copy.deepcopy(default_HPC_2d_description)
    problem_description.symbol_list['CV']=('obj',"feed_con1","feed_con2",\
                                           "purity_con","drain_con")

    rto_algorithm = DCCPWL_ModifierAdaptationTRPenalty()
    options = {
        "homotopy_simulation": True,
        "homotopy_optimization": False,
        "noise_adding_fashion": "added",
        "eta1": global_parameter["eta1"],
        "eta2": global_parameter["eta2"],
        "gamma1": global_parameter["gamma1"],
        "gamma2": global_parameter["gamma2"],
        "gamma3": global_parameter["gamma3"],
        "max_iter": max_iter,
        "sigma": global_parameter["sigma"],
        "max_trust_radius": global_parameter["max_trust_radius"],
        "initial_trust_radius": global_parameter["initial_trust_radius"],
    }
    rto_algorithm.set_algorithm_option(options)

    ffd_perturb = SimpleFiniteDiffPerturbation(perturbation_stepsize, problem_description)

    noise_generator = NoiseGenerator()
    noise_generator.load_noise(noise_filename)

    rto_algorithm.set_problem(
        problem_description=problem_description,
        plant=RTO_Plant_HPC(),
        model_dc_cpwl_functions={'obj':HPC2_2d_cost_QCPWL_subgrad(),
                                "feed_con1":HPC2_2d_feed_con1_QCPWL_subgrad(),
                                "feed_con2":HPC2_2d_feed_con2_QCPWL_subgrad(),
                                "purity_con":HPC2_2d_purity_con_QCPWL_subgrad(),
                                "drain_con":HPC2_2d_drain_con_QCPWL_subgrad(),
                                 "validity_con":HPC2_2d_validity_CPWL_subgrad(),
                                 },
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        nlp_solver_executable=nlp_solver_executable,
        qcqp_solver_executable=qcqp_solver_executable,
        spec_function=specification_func,
        modifier_type=ModifierType.RTO,
        parameter_set=[],
        model_mvs=['TA_Flow', 'MA_Flow'],
        )

    homotopy_var=['TA_Flow', 'MA_Flow', 'GAN_Flow', 'LN_Flow']
    rto_algorithm.plant_simulator.homotopy_var = homotopy_var

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    # draw 2d contour
    num_point = 40
    x_points = np.linspace(200, 400, num_point)
    y_points = np.linspace(1600, 2000, num_point)
    x_points, y_points = np.meshgrid(x_points, y_points)
    z_points = np.zeros((num_point, num_point))
    for i in range(num_point):
        for j in range(num_point):
            output, _ = rto_algorithm.DC_CPWL_RTO_model.simulate(
                {"TA_Flow": x_points[i, j], "MA_Flow": y_points[i, j]}, with_modifier=True)
            if (output["feed_con1"] <= 0) and (output["feed_con2"] <= 0) and \
                    (output["purity_con"] <= 0) and (output["drain_con"] <= 0):
                z_points[i, j] = output["obj"]
            else:
                z_points[i, j] = 3500
    plt.figure()
    plt.contourf(x_points, y_points, z_points, 20)
    plt.colorbar()
    plt.scatter(rto_algorithm.input_history_data[0]['TA_Flow'], \
                rto_algorithm.input_history_data[0]['MA_Flow'], s=10, c='black')
    plt.savefig("pic/HPC2/HPC2_init.png")
    plt.close()

    for iter_no in range(max_iter):
        rto_algorithm.one_step_simulation()

        # draw 2d contour
        num_point = 40
        x_points = np.linspace(200, 400, num_point)
        y_points = np.linspace(1600, 2000, num_point)
        x_points, y_points = np.meshgrid(x_points, y_points)
        z_points = np.zeros((num_point, num_point))
        for i in range(num_point):
            for j in range(num_point):
                output, _ = rto_algorithm.DC_CPWL_RTO_model.simulate(
                    {"TA_Flow": x_points[i, j], "MA_Flow": y_points[i, j]}, with_modifier=True)
                if (output["feed_con1"] <= 0) and (output["feed_con2"] <= 0) and \
                        (output["purity_con"] <= 0) and (output["drain_con"] <= 0):
                    z_points[i, j] = output["obj"]
                else:
                    z_points[i, j] = 3500
        plt.figure()
        plt.contourf(x_points, y_points, z_points, 20)
        plt.colorbar()
        plt.scatter(rto_algorithm.input_history_data[iter_no+1]['TA_Flow'], \
                    rto_algorithm.input_history_data[iter_no+1]['MA_Flow'], s=10, c='black')
        plt.savefig("pic/HPC2/HPC2_iter" + str(iter_no) + ".png")
        plt.close()

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
    noise_filename = "noise/hpc-noise0.txt"
    nlp_solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    qcqp_solver_executable = None
    result_filename_folder="data/hpc2d/"
    if not os.path.exists(result_filename_folder):
        os.makedirs(result_filename_folder)

    # ------------------------------------
    filtering_factor = 1
    starting_point = {'TA_Flow': 300,
                      'MA_Flow': 1750,
                      }

    # ------------------------------------
    perturbation_stepsize = {'TA_Flow': 2,
                'MA_Flow': 2,
                }
    # ------------------------------------
    print_iter_data = False
    max_iter = 10

    # ------------------------------------
    # print("\nTesting MA")
    # result_filename_header = result_filename_folder + "MA_"
    # do_test_MA(perturbation_stepsize, starting_point, filtering_factor, \
    #            noise_filename, nlp_solver_executable, \
    #            print_iter_data, max_iter, \
    #            result_filename_header)
    # ------------------------------------
    # print("\nTesting QCPWL-MA")
    # result_filename_header = result_filename_folder + "QCPWL_MA_"
    # do_test_QCPWL_MA(perturbation_stepsize, starting_point, filtering_factor, \
    #                  noise_filename, nlp_solver_executable, qcqp_solver_executable, \
    #                  print_iter_data, max_iter, \
    #                  result_filename_header)
    # ------------------------------------
    print("\nTesting QCPWL-Subgrad")
    result_filename_header = result_filename_folder + "QCPWL_Subgrad_"
    do_test_QCPWL_subgrad(perturbation_stepsize, starting_point, filtering_factor, \
                          noise_filename, nlp_solver_executable, qcqp_solver_executable, \
                          print_iter_data, max_iter, \
                          result_filename_header)
    # ------------------------------------
    print("\nTesting QCPWL-Subgrad-TR")
    result_filename_header = result_filename_folder + "QCPWL_Subgrad_PenaltyTR_"
    do_test_QCPWL_subgrad_PenaltyTR(perturbation_stepsize, starting_point, filtering_factor, \
                     noise_filename, nlp_solver_executable, qcqp_solver_executable, \
                     print_iter_data, max_iter, \
                     result_filename_header)

def compare_MA_and_QCPWL(pic_name):
    pic_constant = 0.39370
    global_font_size = 15
    global_tick_size = 15
    global_linewidth = 1
    global_point_size = 4
    font_factor = np.sqrt(1 / pic_constant)
    fig = plt.figure(figsize=(13 * pic_constant, 30 * pic_constant))
    font2 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': global_font_size,
             }
    font3 = {'family': 'Times New Roman',
             'weight': 'bold',
             'size': global_font_size,
             }
    font_legend = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 16/font_factor,
             }
    font_title = {'family': 'Times New Roman',
                  'size': global_font_size,
                  'weight': 'bold'
                  }
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    algorighms = ['QCPWL_MA', 'MA','QCPWL_Subgrad']
    algo_name_dict = {
        'QCPWL_MA':'QCPWL-MA',
        'MA':"MA",
        'QCPWL_Subgrad':'QCPWL_Subgrad',
    }
    algo_color = {
        'QCPWL_MA':'red',
        'MA':"blue",
        'QCPWL_Subgrad':"black"
    }
    output_measurements = ['obj','Purity_GAN','Drain_Flow']
    inputs = ['TA_Flow', 'MA_Flow']
    mv_label_dict = {
        'TA_Flow': 'Turbine air flowrate (kmol/h)',
        'MA_Flow': 'Main air flowrate (kmol/h)',
    }
    y_label_dict = {
        'obj':'Cost',
        'Purity_GAN':'Purity of GAN',
        'Drain_Flow':'Drain flowrate (kmol/h)',
    }
    for op_index, op in enumerate(output_measurements):
        plt.subplot(5, 1, op_index + 1)
        for algo_index, algo in enumerate(algorighms):
            model_data = pandas.read_csv("data/hpc2d/" + algo + '_model_data' + ".txt", sep='\t', index_col=0, header=0)
            plant_data = pandas.read_csv("data/hpc2d/" + algo + '_plant_data' + ".txt", sep='\t', index_col=0, header=0)
            time = model_data.index
            plt.plot(time[:30], plant_data.loc[:30, op], linewidth=global_linewidth, label=algo,
                                       color=algo_color[algo],
                                       linestyle='-')
        plt.legend(prop=font_legend)
        plt.ylabel(y_label_dict[op], font2)
        plt.xlabel("RTO iteration", font2)
        plt.xticks([0,10,20,30])
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        plt.tick_params(labelsize=global_tick_size)
        ax = plt.gca()
        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontname('Times New Roman') for label in labels]
        plt.tick_params(top= ['on', 'true'], right= ['on', 'true'], which='both')
        plt.grid(ls='--',axis='y')

    for op_index, op in enumerate(inputs):
        plt.subplot(5, 1, op_index + 4)
        for algo_index, algo in enumerate(algorighms):
            model_data = pandas.read_csv("data/hpc2d/" + algo + '_model_data' + ".txt", sep='\t', index_col=0, header=0)
            input_data = pandas.read_csv("data/hpc2d/" + algo + '_input_data' + ".txt", sep='\t', index_col=0, header=0)
            time = model_data.index
            plt.plot(time[:30], input_data.loc[1:31, op], linewidth=global_linewidth, label=algo,
                                       color=algo_color[algo],
                                       linestyle='-')
        plt.legend(prop=font_legend)
        plt.ylabel(mv_label_dict[op], font2)
        plt.xlabel("RTO iteration", font2)
        plt.xticks([0,10,20,30])
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        plt.tick_params(labelsize=global_tick_size)
        ax = plt.gca()
        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontname('Times New Roman') for label in labels]
        plt.tick_params(top= ['on', 'true'], right= ['on', 'true'], which='both')
        plt.grid(ls='--',axis='y')

    plt.tight_layout(pad=0.5,w_pad=0.5, h_pad=0.5)
    plt.savefig(pic_name,dpi=600)
    plt.close()

if __name__ == "__main__":
    # plant_simulator_test()
    # generate_noise_file()
    algo_test_main()
    compare_MA_and_QCPWL("pic/HPC2/compare_MA_and_QCPWL.png")