import numpy as np
import pickle
import pandas
from rtolib.core.dc_cpwl_model import QuadraticBoostedDCCPWLFunction,\
QuadraticBoostedDCCPWLFunctionSubgrad
from rtolib.core.algo import ModifierAdaptation
from rtolib.core.algo.dccpwl_ma import DCCPWL_ModifierAdaptation
from rtolib.core.algo.subgrad import DCCPWL_ModifierAdaptationSubgrad,\
            DCCPWL_ModifierAdaptationTRPenalty, QCPWL_Subgrad_MINLP
from rtolib.model.HPC import RTO_Plant_HPC, default_HPC_2d_description, RTO_Mismatched_HPC
from rtolib.core import NoiseGenerator, ModifierType,\
                    SimpleFiniteDiffPerturbation
from rtolib.core.solve import PyomoSimulator, PyomoOptimizer
from rtolib.core.dc_cpwl_model import QuadraticBoostedDCCPWL_RTOObject
import os
import copy
from rtolib.util.misc import save_iteration_data_in_dict
from rtolib.core.pyomo_model import SolverFactory
import matplotlib.pyplot as plt
from rtolib.core.cpwl_minlp import QCPWLFunction_MINLP
import sys
import time

# set specification
def specification_func(iter):
    ret = {}
    if iter % 30 >= 20:
        ret['GAN_Flow'] = 750
        ret['LN_Flow'] = 100
    elif iter % 30 >= 10:
        ret['GAN_Flow'] = 800
        ret['LN_Flow'] = 100
    else:
        ret['GAN_Flow'] = 780
        ret['LN_Flow'] = 100
    # ret['GAN_Flow'] = 750
    # ret['LN_Flow'] = 100
    return ret

global_parameter={
        "eta1": 0.01,
        "eta2": 0.9,
        "gamma1": 0.5,
        "gamma2": 1,
        "gamma3": 2,
        "sigma": 10,
        "max_trust_radius": 500,
        "initial_trust_radius": 10,
}

def plant_simulator_test():
    plant_simulator = PyomoSimulator(RTO_Mismatched_HPC())
    plant_simulator.build(default_HPC_2d_description)
    default_options = {'max_iter': 500,
                       "tol": 1e-10}
    nlp_solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    solver1 = SolverFactory('ipopt', executable=nlp_solver_executable)
    plant_simulator.set_solver(solver1, tee=True, default_options=default_options)
    homotopy_var = ['TA_Flow', 'MA_Flow', 'GAN_Flow', 'LN_Flow']
    plant_simulator.homotopy_var = homotopy_var
    point={
        'TA_Flow': 500,
        'MA_Flow': 1600,
        'GAN_Flow': 800,
        'LN_Flow': 100,
    }
    outputs, solve_status = plant_simulator.simulate(point, param_values=None,
                                            use_homo=True)
    print(outputs)

def get_plant_optimum():
    plant_optimizer = PyomoOptimizer(RTO_Plant_HPC())
    plant_simulator = PyomoSimulator(RTO_Plant_HPC())
    problem_description = copy.deepcopy(default_HPC_2d_description)
    problem_description.symbol_list['CV'] = ('obj', "feed_con1", "feed_con2", \
                                             "purity_con", "drain_con")
    plant_optimizer.build(problem_description)
    plant_simulator.build(problem_description)
    nlp_solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    solver = SolverFactory('ipopt', executable=nlp_solver_executable)
    default_options = {'max_iter': 500,
                       "tol": 1e-10}
    plant_optimizer.set_solver(solver, tee=False, default_options=default_options)
    plant_simulator.set_solver(solver, tee=False, default_options=default_options)
    homotopy_var = ['TA_Flow', 'MA_Flow', 'GAN_Flow', 'LN_Flow']
    plant_simulator.homotopy_var = homotopy_var
    current_parameter_value=\
    {'obj_eps': 0, 'TA_Flow_obj_lam': 0,
        'MA_Flow_obj_lam': 0, 'feed_con1_eps': 0,
        'TA_Flow_feed_con1_lam': 0, 'MA_Flow_feed_con1_lam': 0,
        'feed_con2_eps': 0, 'TA_Flow_feed_con2_lam': 0,
        'MA_Flow_feed_con2_lam': 0, 'purity_con_eps': 0,
        'TA_Flow_purity_con_lam': 0, 'MA_Flow_purity_con_lam': 0,
        'drain_con_eps': 0, 'TA_Flow_drain_con_lam': 0,
        'MA_Flow_drain_con_lam': 0}

    for iter in [5,15,25]:
        print("iteration: ", iter)
        input_values = specification_func(iter)
        optimized_input, solve_status = plant_optimizer.optimize(input_values,
                                                                      param_values=current_parameter_value,
                                                                      use_homo=False)
        outputs, solve_status = plant_simulator.simulate(optimized_input, param_values=None,
                                                         use_homo=True)
        print(optimized_input)
        print(outputs)


class HPC3_3d_validity_CPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "hpc_model/HPC3_3d_validity.bin"
        self.load_from_cpwl_model_file(cpwl_file)
        self.input_variables = ['TA_Flow', 'MA_Flow', 'GAN_Flow']
        self.parameters = []

class HPC3_3d_cost_QCPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "hpc_model/HPC3_3d_cpwl_processed_output_0.bin"
        quadr_file = "hpc_model/HPC3_3d_quadr_processed_output_0.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow', 'GAN_Flow']
        self.parameters = []

class HPC3_3d_feed_con1_QCPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "hpc_model/HPC3_3d_cpwl_processed_output_1.bin"
        quadr_file = "hpc_model/HPC3_3d_quadr_processed_output_1.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow', 'GAN_Flow']
        self.parameters = []

class HPC3_3d_feed_con2_QCPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "hpc_model/HPC3_3d_cpwl_processed_output_2.bin"
        quadr_file = "hpc_model/HPC3_3d_quadr_processed_output_2.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow', 'GAN_Flow']
        self.parameters = []

class HPC3_3d_purity_con_QCPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "hpc_model/HPC3_3d_cpwl_processed_output_3.bin"
        quadr_file = "hpc_model/HPC3_3d_quadr_processed_output_3.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow', 'GAN_Flow']
        self.parameters = []

class HPC3_3d_drain_con_QCPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "hpc_model/HPC3_3d_cpwl_processed_output_4.bin"
        quadr_file = "hpc_model/HPC3_3d_quadr_processed_output_4.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow', 'GAN_Flow']
        self.parameters = []
        
class HPC3_3d_validity_CPWL_subgrad(QuadraticBoostedDCCPWLFunctionSubgrad):
    def __init__(self):
        cpwl_file = "hpc_model/HPC3_3d_validity.bin"
        self.load_from_cpwl_model_file(cpwl_file)
        self.input_variables = ['TA_Flow', 'MA_Flow', 'GAN_Flow']
        self.parameters = []

class HPC3_3d_cost_QCPWL_subgrad(QuadraticBoostedDCCPWLFunctionSubgrad):
    def __init__(self):
        cpwl_file = "hpc_model/HPC3_3d_cpwl_processed_output_0.bin"
        quadr_file = "hpc_model/HPC3_3d_quadr_processed_output_0.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow', 'GAN_Flow']
        self.parameters = []

class HPC3_3d_feed_con1_QCPWL_subgrad(QuadraticBoostedDCCPWLFunctionSubgrad):
    def __init__(self):
        cpwl_file = "hpc_model/HPC3_3d_cpwl_processed_output_1.bin"
        quadr_file = "hpc_model/HPC3_3d_quadr_processed_output_1.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow', 'GAN_Flow']
        self.parameters = []

class HPC3_3d_feed_con2_QCPWL_subgrad(QuadraticBoostedDCCPWLFunctionSubgrad):
    def __init__(self):
        cpwl_file = "hpc_model/HPC3_3d_cpwl_processed_output_2.bin"
        quadr_file = "hpc_model/HPC3_3d_quadr_processed_output_2.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow', 'GAN_Flow']
        self.parameters = []

class HPC3_3d_purity_con_QCPWL_subgrad(QuadraticBoostedDCCPWLFunctionSubgrad):
    def __init__(self):
        cpwl_file = "hpc_model/HPC3_3d_cpwl_processed_output_3.bin"
        quadr_file = "hpc_model/HPC3_3d_quadr_processed_output_3.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow', 'GAN_Flow']
        self.parameters = []

class HPC3_3d_drain_con_QCPWL_subgrad(QuadraticBoostedDCCPWLFunctionSubgrad):
    def __init__(self):
        cpwl_file = "hpc_model/HPC3_3d_cpwl_processed_output_4.bin"
        quadr_file = "hpc_model/HPC3_3d_quadr_processed_output_4.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow', 'GAN_Flow']
        self.parameters = []

class HPC3_3d_validity_MINLP(QCPWLFunction_MINLP):
    def __init__(self):
        cpwl_file = "hpc_model/HPC3_3d_validity.bin"
        self.load_from_cpwl_model_file(cpwl_file)
        self.input_variables = ['TA_Flow', 'MA_Flow', 'GAN_Flow']
        self.parameters = []
        self.bigM = 100000

class HPC3_3d_cost_MINLP(QCPWLFunction_MINLP):
    def __init__(self):
        cpwl_file = "hpc_model/HPC3_3d_cpwl_processed_output_0.bin"
        quadr_file = "hpc_model/HPC3_3d_quadr_processed_output_0.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow', 'GAN_Flow']
        self.parameters = []
        self.bigM = 100000

class HPC3_3d_feed_con1_MINLP(QCPWLFunction_MINLP):
    def __init__(self):
        cpwl_file = "hpc_model/HPC3_3d_cpwl_processed_output_1.bin"
        quadr_file = "hpc_model/HPC3_3d_quadr_processed_output_1.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow', 'GAN_Flow']
        self.parameters = []
        self.bigM = 100000

class HPC3_3d_feed_con2_MINLP(QCPWLFunction_MINLP):
    def __init__(self):
        cpwl_file = "hpc_model/HPC3_3d_cpwl_processed_output_2.bin"
        quadr_file = "hpc_model/HPC3_3d_quadr_processed_output_2.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow', 'GAN_Flow']
        self.parameters = []
        self.bigM = 100000

class HPC3_3d_purity_con_MINLP(QCPWLFunction_MINLP):
    def __init__(self):
        cpwl_file = "hpc_model/HPC3_3d_cpwl_processed_output_3.bin"
        quadr_file = "hpc_model/HPC3_3d_quadr_processed_output_3.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow', 'GAN_Flow']
        self.parameters = []
        self.bigM = 100

class HPC3_3d_drain_con_MINLP(QCPWLFunction_MINLP):
    def __init__(self):
        cpwl_file = "hpc_model/HPC3_3d_cpwl_processed_output_4.bin"
        quadr_file = "hpc_model/HPC3_3d_quadr_processed_output_4.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow', 'GAN_Flow']
        self.parameters = []
        self.bigM = 100000

def generate_noise_file():
    noise_level = {
        "obj": 2,
        "feed_con1": 0,
        "feed_con2": 0,
        "purity_con": 1e-5,
        "drain_con": 0.1,
        "validity_con":0,
    }
    # OK
    # noise_level = {
    #     "obj": 2,
    #     "feed_con1": 0,
    #     "feed_con2": 0,
    #     "purity_con": 1e-5,
    #     "drain_con": 0.1,
    #     "validity_con": 0,
    # }
    max_iter = 200
    noise_generator = NoiseGenerator(noise_level)
    noise_generator.generate_noise(max_iter, 6)
    noise_generator.save_noise("noise/hpc-noise-for-paper.txt")
    
def get_plant_contour_data():
    plant_simulator = PyomoSimulator(RTO_Plant_HPC())
    plant_simulator.build(default_HPC_2d_description)
    default_options = {'max_iter': 500,
                       "tol": 1e-10}
    nlp_solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    solver1 = SolverFactory('ipopt', executable=nlp_solver_executable)
    plant_simulator.set_solver(solver1, tee=False, default_options=default_options)
    homotopy_var = ['TA_Flow', 'MA_Flow', 'GAN_Flow', 'LN_Flow']
    plant_simulator.homotopy_var = homotopy_var

    HPC_validity = HPC3_3d_validity_CPWL()
    model_dc_cpwl_functions = {'obj': HPC3_3d_cost_QCPWL(),
                               "feed_con1": HPC3_3d_feed_con1_QCPWL(),
                               "feed_con2": HPC3_3d_feed_con2_QCPWL(),
                               "purity_con": HPC3_3d_purity_con_QCPWL(),
                               "drain_con": HPC3_3d_drain_con_QCPWL(),
                               }
    counter = 0
    for gan_flow in [750, 770, 800]:
        num_point = 20
        x_points = np.linspace(200, 400, num_point)
        y_points = np.linspace(1600, 2000, num_point)
        x_points, y_points = np.meshgrid(x_points, y_points)
        z_points = {}
        for output_name in model_dc_cpwl_functions.keys():
            z_points[output_name] = np.zeros((num_point, num_point))
        for i in range(num_point):
            for j in range(num_point):
                if counter % 10 == 0:
                    print(counter)
                counter += 1
                output = HPC_validity.simulate(
                    {"TA_Flow": x_points[i, j], "MA_Flow": y_points[i, j],
                     "GAN_Flow": gan_flow})
                if output < 0.5:
                    for output_name in model_dc_cpwl_functions.keys():
                        if output_name == "obj":
                            z_points[output_name][i, j] = 4500
                        else:
                            z_points[output_name][i, j] = 0
                else:
                    point = {
                        'TA_Flow': x_points[i, j],
                        'MA_Flow': y_points[i, j],
                        'GAN_Flow': gan_flow,
                        'LN_Flow': 100,
                    }
                    outputs, solve_status = plant_simulator.simulate(point, param_values=None,
                                                                     use_homo=True)
                    for output_name in model_dc_cpwl_functions.keys():
                        z_points[output_name][i, j] = outputs[output_name]
        for output_name in model_dc_cpwl_functions.keys():
            data_file = "plant_data/"+output_name+"_%d"%num_point+"_%d"%gan_flow+".bin"
            with open(data_file, "wb") as fp:
                pickle.dump(z_points[output_name], fp)

def draw_plant_model_mismatch():
    size = np.fromfile("hpc_model/HPC3_processed_size_3d.dat", dtype=int)
    X = np.fromfile("hpc_model/HPC3_processed_X_3d.dat", dtype=float).reshape((size[0], size[1]))
    Y = np.fromfile("hpc_model/HPC3_processed_Y_3d.dat", dtype=float).reshape((size[2], size[3]))
    output_name_and_index = {
        "obj": 0,
        "feed_con1": 1,
        "feed_con2": 2,
        "purity_con": 3,
        "drain_con": 4,
    }
    model_dc_cpwl_functions = {'obj': HPC3_3d_cost_QCPWL(),
                               "feed_con1": HPC3_3d_feed_con1_QCPWL(),
                               "feed_con2": HPC3_3d_feed_con2_QCPWL(),
                               "purity_con": HPC3_3d_purity_con_QCPWL(),
                               "drain_con": HPC3_3d_drain_con_QCPWL(),
                               }
    mvs = ('TA_Flow', 'MA_Flow', 'GAN_Flow')
    cvs = ('obj',"feed_con1","feed_con2", "purity_con","drain_con")
    DC_CPWL_RTO_model = QuadraticBoostedDCCPWL_RTOObject(model_dc_cpwl_functions, mvs, cvs, [])
    DC_CPWL_RTO_model.build()
    nlp_solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    solver = SolverFactory('ipopt', executable=nlp_solver_executable)
    DC_CPWL_RTO_model.set_solver(solver, tee=False, default_options={})

    for output_name in output_name_and_index.keys():
        plant_y = np.zeros((size[0], ))
        model_y = np.zeros((size[0], ))
        for i in range(size[0]):
            plant_y[i] = Y[i, output_name_and_index[output_name]]
            x = {
                "TA_Flow": X[i, 0],
                "MA_Flow": X[i, 1],
                "GAN_Flow": X[i, 2]
            }
            model_output, _ = DC_CPWL_RTO_model.simulate(x)
            model_y[i] = model_output[output_name]
        plt.figure()
        plt.plot(plant_y-model_y)
        plt.title(output_name)
        plt.savefig("pic/HPC3/plant_model_mismatch_"+output_name+".png")
        plt.close()

def draw_plant_model_mismatch2():
    model_dc_cpwl_functions = {'obj': HPC3_3d_cost_QCPWL(),
                               "feed_con1": HPC3_3d_feed_con1_QCPWL(),
                               "feed_con2": HPC3_3d_feed_con2_QCPWL(),
                               "purity_con": HPC3_3d_purity_con_QCPWL(),
                               "drain_con": HPC3_3d_drain_con_QCPWL(),
                               }
    for gan_flow in [750, 770, 800]:
        num_point = 20
        x_points = np.linspace(200, 400, num_point)
        y_points = np.linspace(1600, 2000, num_point)
        x_points, y_points = np.meshgrid(x_points, y_points)
        for output_name in model_dc_cpwl_functions.keys():
            data_file = "hpc_model/" + output_name + "_%d" % num_point + \
                "_%d" % gan_flow + ".bin"
            with open(data_file, "rb") as fp:
                plant_z = pickle.load(fp)
            model_z = np.zeros((num_point, num_point))
            for i in range(num_point):
                for j in range(num_point):
                    model_z[i, j] = model_dc_cpwl_functions[output_name].simulate(
                        {"TA_Flow": x_points[i, j], "MA_Flow": y_points[i, j], "GAN_Flow": gan_flow}
                    )
            plt.figure(figsize=(15,6))
            plt.subplot(131)
            plt.contourf(x_points, y_points, plant_z, 20)
            plt.colorbar()
            plt.title(output_name+"_plant")
            plt.subplot(132)
            plt.contourf(x_points, y_points, model_z, 20)
            plt.colorbar()
            plt.title(output_name + "_model")
            plt.subplot(133)
            plt.contourf(x_points, y_points, model_z-plant_z, 20)
            plt.colorbar()
            plt.title(output_name + "_error")
            plt.savefig("pic/HPC3/plant_model_mismatch2_"+output_name+"_ganflow%d"%gan_flow+".png")
            plt.close()

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
        "homotopy_optimization": False,
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
        spec_function=None,
        modifier_type=ModifierType.RTO,
        )

    homotopy_var=['TA_Flow', 'MA_Flow', 'GAN_Flow', 'LN_Flow']
    rto_algorithm.set_homotopy_var(homotopy_var)

    rto_algorithm.spec_function = specification_func

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
        plant=RTO_Plant_HPC(),
        model_dc_cpwl_functions={'obj':HPC3_3d_cost_QCPWL(),
                                "feed_con1":HPC3_3d_feed_con1_QCPWL(),
                                "feed_con2":HPC3_3d_feed_con2_QCPWL(),
                                "purity_con":HPC3_3d_purity_con_QCPWL(),
                                "drain_con":HPC3_3d_drain_con_QCPWL(),
                                 "validity_con":HPC3_3d_validity_CPWL(),
                                 },
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        nlp_solver_executable=nlp_solver_executable,
        qcqp_solver_executable=qcqp_solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        parameter_set=[],
        model_mvs=['TA_Flow', 'MA_Flow', 'GAN_Flow'],
        )

    homotopy_var=['TA_Flow', 'MA_Flow', 'GAN_Flow', 'LN_Flow']
    rto_algorithm.plant_simulator.homotopy_var = homotopy_var

    rto_algorithm.spec_function = specification_func

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    # draw 2d contour
    num_point = 40
    x_points = np.linspace(200, 400, num_point)
    y_points = np.linspace(1600, 2000, num_point)
    x_points, y_points = np.meshgrid(x_points, y_points)
    z_points = np.zeros((num_point, num_point))
    gan_flow = specification_func(0)['GAN_Flow']
    for i in range(num_point):
        for j in range(num_point):
            output, _ = rto_algorithm.DC_CPWL_RTO_model.simulate(
                {"TA_Flow": x_points[i, j], "MA_Flow": y_points[i, j],
                 "GAN_Flow": gan_flow}, with_modifier=True)
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
    plt.savefig("pic/HPC3/HPC3_init.png")
    plt.close()

    for iter_no in range(max_iter):
        rto_algorithm.one_step_simulation()

        # draw 2d contour
        num_point = 40
        x_points = np.linspace(200, 400, num_point)
        y_points = np.linspace(1600, 2000, num_point)
        x_points, y_points = np.meshgrid(x_points, y_points)
        z_points = np.zeros((num_point, num_point))
        gan_flow = specification_func(iter_no)['GAN_Flow']
        for i in range(num_point):
            for j in range(num_point):
                output, _ = rto_algorithm.DC_CPWL_RTO_model.simulate(
                    {"TA_Flow": x_points[i, j], "MA_Flow": y_points[i, j],
                     "GAN_Flow":gan_flow}, with_modifier=True)
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
        plt.savefig("pic/HPC3/HPC3_iter" + str(iter_no) + ".png")
        plt.close()

        if print_iter_data:
            print(rto_algorithm.model_history_data)
            print(rto_algorithm.plant_history_data)
            print(rto_algorithm.input_history_data)

    # save data
    save_iteration_data_in_dict(rto_algorithm.model_history_data, result_filename_header + "model_data.txt")
    save_iteration_data_in_dict(rto_algorithm.plant_history_data, result_filename_header + "plant_data.txt")
    save_iteration_data_in_dict(rto_algorithm.input_history_data, result_filename_header + "input_data.txt")


def do_test_QCPWL_Subgrad(perturbation_stepsize, starting_point, filtering_factor, \
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
        model_dc_cpwl_functions={'obj':HPC3_3d_cost_QCPWL_subgrad(),
                                "feed_con1":HPC3_3d_feed_con1_QCPWL_subgrad(),
                                "feed_con2":HPC3_3d_feed_con2_QCPWL_subgrad(),
                                "purity_con":HPC3_3d_purity_con_QCPWL_subgrad(),
                                "drain_con":HPC3_3d_drain_con_QCPWL_subgrad(),
                                 "validity_con":HPC3_3d_validity_CPWL_subgrad(),
                                 },
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        nlp_solver_executable=nlp_solver_executable,
        qcqp_solver_executable=qcqp_solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        parameter_set=[],
        model_mvs=['TA_Flow', 'MA_Flow', 'GAN_Flow'],
        )

    homotopy_var=['TA_Flow', 'MA_Flow', 'GAN_Flow', 'LN_Flow']
    rto_algorithm.plant_simulator.homotopy_var = homotopy_var

    rto_algorithm.spec_function = specification_func

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    # draw 2d contour
    num_point = 40
    x_points = np.linspace(200, 400, num_point)
    y_points = np.linspace(1600, 2000, num_point)
    x_points, y_points = np.meshgrid(x_points, y_points)
    z_points = np.zeros((num_point, num_point))
    gan_flow = specification_func(0)['GAN_Flow']
    for i in range(num_point):
        for j in range(num_point):
            output, _ = rto_algorithm.DC_CPWL_RTO_model.simulate(
                {"TA_Flow": x_points[i, j], "MA_Flow": y_points[i, j],
                 "GAN_Flow": gan_flow}, with_modifier=True)
            z_points[i, j] = output["obj"] / problem_description.scaling_factors["obj"] + \
                             global_parameter['sigma'] * ( \
                                         max(output["feed_con1"],0) / problem_description.scaling_factors["feed_con1"] + \
                                         max(output["feed_con2"],0) / problem_description.scaling_factors["feed_con2"] + \
                                         max(output["purity_con"],0) / problem_description.scaling_factors["purity_con"] + \
                                         max(output["drain_con"],0) / problem_description.scaling_factors["drain_con"]
                                 )
    plt.figure()
    plt.contourf(x_points, y_points, z_points, 20)
    plt.colorbar()
    plt.scatter(rto_algorithm.input_history_data[0]['TA_Flow'], \
                rto_algorithm.input_history_data[0]['MA_Flow'], s=10, c='black')
    plt.savefig("pic/HPC3_subgrad/HPC3_init.png")
    plt.close()

    for iter_no in range(max_iter):
        rto_algorithm.one_step_simulation()

        # draw 2d contour
        num_point = 40
        x_points = np.linspace(200, 400, num_point)
        y_points = np.linspace(1600, 2000, num_point)
        x_points, y_points = np.meshgrid(x_points, y_points)
        z_points = np.zeros((num_point, num_point))
        gan_flow = specification_func(iter_no)['GAN_Flow']
        for i in range(num_point):
            for j in range(num_point):
                output, _ = rto_algorithm.DC_CPWL_RTO_model.simulate(
                    {"TA_Flow": x_points[i, j], "MA_Flow": y_points[i, j],
                     "GAN_Flow":gan_flow}, with_modifier=True)
                z_points[i, j] = output["obj"] / problem_description.scaling_factors["obj"] + \
                                 global_parameter['sigma'] * ( \
                                             max(output["feed_con1"],0) / problem_description.scaling_factors["feed_con1"] + \
                                             max(output["feed_con2"],0) / problem_description.scaling_factors["feed_con2"] + \
                                             max(output["purity_con"],0) / problem_description.scaling_factors["purity_con"] + \
                                             max(output["drain_con"],0) / problem_description.scaling_factors["drain_con"]
                                     )
        plt.figure()
        plt.contourf(x_points, y_points, z_points, 20)
        plt.colorbar()
        plt.scatter(rto_algorithm.input_history_data[iter_no+1]['TA_Flow'], \
                    rto_algorithm.input_history_data[iter_no+1]['MA_Flow'], s=10, c='black')
        plt.savefig("pic/HPC3_subgrad/HPC3_iter" + str(iter_no) + ".png")
        plt.close()

        if print_iter_data:
            print(rto_algorithm.model_history_data)
            print(rto_algorithm.plant_history_data)
            print(rto_algorithm.input_history_data)

    # save data
    save_iteration_data_in_dict(rto_algorithm.model_history_data, result_filename_header + "model_data.txt")
    save_iteration_data_in_dict(rto_algorithm.plant_history_data, result_filename_header + "plant_data.txt")
    save_iteration_data_in_dict(rto_algorithm.input_history_data, result_filename_header + "input_data.txt")

def do_test_QCPWL_Subgrad_MINLP(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, nlp_solver_executable, minlp_solver_executable,\
                     print_iter_data, max_iter,\
               result_filename_header):
    problem_description = copy.deepcopy(default_HPC_2d_description)
    problem_description.symbol_list['CV']=('obj',"feed_con1","feed_con2",\
                                           "purity_con","drain_con")

    rto_algorithm = QCPWL_Subgrad_MINLP()
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
        model_dc_cpwl_functions={'obj':HPC3_3d_cost_MINLP(),
                                "feed_con1":HPC3_3d_feed_con1_MINLP(),
                                "feed_con2":HPC3_3d_feed_con2_MINLP(),
                                "purity_con":HPC3_3d_purity_con_MINLP(),
                                "drain_con":HPC3_3d_drain_con_MINLP(),
                                 "validity_con":HPC3_3d_validity_MINLP(),
                                 },
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        nlp_solver_executable=nlp_solver_executable,
        minlp_solver_executable=minlp_solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        parameter_set=[],
        model_mvs=['TA_Flow', 'MA_Flow', 'GAN_Flow'],
        )

    homotopy_var=['TA_Flow', 'MA_Flow', 'GAN_Flow', 'LN_Flow']
    rto_algorithm.plant_simulator.homotopy_var = homotopy_var

    rto_algorithm.spec_function = specification_func

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    # draw 2d contour
    num_point = 40
    x_points = np.linspace(200, 400, num_point)
    y_points = np.linspace(1600, 2000, num_point)
    x_points, y_points = np.meshgrid(x_points, y_points)
    z_points = np.zeros((num_point, num_point))
    gan_flow = specification_func(0)['GAN_Flow']
    for i in range(num_point):
        for j in range(num_point):
            output, _ = rto_algorithm.DC_CPWL_RTO_model.simulate(
                {"TA_Flow": x_points[i, j], "MA_Flow": y_points[i, j],
                 "GAN_Flow": gan_flow}, with_modifier=True)
            z_points[i, j] = output["obj"] / problem_description.scaling_factors["obj"] + \
                             global_parameter['sigma'] * ( \
                                         max(output["feed_con1"],0) / problem_description.scaling_factors["feed_con1"] + \
                                         max(output["feed_con2"],0) / problem_description.scaling_factors["feed_con2"] + \
                                         max(output["purity_con"],0) / problem_description.scaling_factors["purity_con"] + \
                                         max(output["drain_con"],0) / problem_description.scaling_factors["drain_con"]
                                 )
    plt.figure()
    plt.contourf(x_points, y_points, z_points, 20)
    plt.colorbar()
    plt.scatter(rto_algorithm.input_history_data[0]['TA_Flow'], \
                rto_algorithm.input_history_data[0]['MA_Flow'], s=10, c='black')
    plt.savefig("pic/HPC3_subgrad_MINLP/HPC3_init.png")
    plt.close()

    for iter_no in range(max_iter):
        rto_algorithm.one_step_simulation()

        # draw 2d contour
        num_point = 40
        x_points = np.linspace(200, 400, num_point)
        y_points = np.linspace(1600, 2000, num_point)
        x_points, y_points = np.meshgrid(x_points, y_points)
        z_points = np.zeros((num_point, num_point))
        gan_flow = specification_func(iter_no)['GAN_Flow']
        for i in range(num_point):
            for j in range(num_point):
                output, _ = rto_algorithm.DC_CPWL_RTO_model.simulate(
                    {"TA_Flow": x_points[i, j], "MA_Flow": y_points[i, j],
                     "GAN_Flow":gan_flow}, with_modifier=True)
                z_points[i, j] = output["obj"] / problem_description.scaling_factors["obj"] + \
                                 global_parameter['sigma'] * ( \
                                             max(output["feed_con1"],0) / problem_description.scaling_factors["feed_con1"] + \
                                             max(output["feed_con2"],0) / problem_description.scaling_factors["feed_con2"] + \
                                             max(output["purity_con"],0) / problem_description.scaling_factors["purity_con"] + \
                                             max(output["drain_con"],0) / problem_description.scaling_factors["drain_con"]
                                     )
        plt.figure()
        plt.contourf(x_points, y_points, z_points, 20)
        plt.colorbar()
        plt.scatter(rto_algorithm.input_history_data[iter_no+1]['TA_Flow'], \
                    rto_algorithm.input_history_data[iter_no+1]['MA_Flow'], s=10, c='black')
        plt.savefig("pic/HPC3_subgrad_MINLP/HPC3_iter" + str(iter_no) + ".png")
        plt.close()

        if print_iter_data:
            print(rto_algorithm.model_history_data)
            print(rto_algorithm.plant_history_data)
            print(rto_algorithm.input_history_data)

    # save data
    save_iteration_data_in_dict(rto_algorithm.model_history_data, result_filename_header + "model_data.txt")
    save_iteration_data_in_dict(rto_algorithm.plant_history_data, result_filename_header + "plant_data.txt")
    save_iteration_data_in_dict(rto_algorithm.input_history_data, result_filename_header + "input_data.txt")


def do_test_QCPWL_Subgrad_PenaltyTR(perturbation_stepsize, starting_point, filtering_factor, \
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
        model_dc_cpwl_functions={'obj':HPC3_3d_cost_QCPWL_subgrad(),
                                "feed_con1":HPC3_3d_feed_con1_QCPWL_subgrad(),
                                "feed_con2":HPC3_3d_feed_con2_QCPWL_subgrad(),
                                "purity_con":HPC3_3d_purity_con_QCPWL_subgrad(),
                                "drain_con":HPC3_3d_drain_con_QCPWL_subgrad(),
                                 "validity_con":HPC3_3d_validity_CPWL_subgrad(),
                                 },
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        nlp_solver_executable=nlp_solver_executable,
        qcqp_solver_executable=qcqp_solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        parameter_set=[],
        model_mvs=['TA_Flow', 'MA_Flow', 'GAN_Flow'],
        )

    homotopy_var=['TA_Flow', 'MA_Flow', 'GAN_Flow', 'LN_Flow']
    rto_algorithm.plant_simulator.homotopy_var = homotopy_var

    rto_algorithm.spec_function = specification_func

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    # draw 2d contour
    num_point = 40
    x_points = np.linspace(200, 400, num_point)
    y_points = np.linspace(1600, 2000, num_point)
    x_points, y_points = np.meshgrid(x_points, y_points)
    z_points = np.zeros((num_point, num_point))
    gan_flow = specification_func(0)['GAN_Flow']
    for i in range(num_point):
        for j in range(num_point):
            output, _ = rto_algorithm.DC_CPWL_RTO_model.simulate(
                {"TA_Flow": x_points[i, j], "MA_Flow": y_points[i, j],
                 "GAN_Flow": gan_flow}, with_modifier=True)
            z_points[i, j] = output["obj"] / problem_description.scaling_factors["obj"] + \
                             global_parameter['sigma'] * ( \
                                         max(output["feed_con1"],0) / problem_description.scaling_factors["feed_con1"] + \
                                         max(output["feed_con2"],0) / problem_description.scaling_factors["feed_con2"] + \
                                         max(output["purity_con"],0) / problem_description.scaling_factors["purity_con"] + \
                                         max(output["drain_con"],0) / problem_description.scaling_factors["drain_con"]
                                 )
    plt.figure()
    plt.contourf(x_points, y_points, z_points, 20)
    plt.colorbar()
    plt.scatter(rto_algorithm.input_history_data[0]['TA_Flow'], \
                rto_algorithm.input_history_data[0]['MA_Flow'], s=10, c='black')
    plt.savefig("pic/HPC3/HPC3_init.png")
    plt.close()

    while rto_algorithm.iter_count <= max_iter:
        rto_algorithm.one_step_simulation()

        # draw 2d contour
        num_point = 40
        x_points = np.linspace(200, 400, num_point)
        y_points = np.linspace(1600, 2000, num_point)
        x_points, y_points = np.meshgrid(x_points, y_points)
        z_points = np.zeros((num_point, num_point))
        gan_flow = specification_func(iter_no)['GAN_Flow']
        for i in range(num_point):
            for j in range(num_point):
                output, _ = rto_algorithm.DC_CPWL_RTO_model.simulate(
                    {"TA_Flow": x_points[i, j], "MA_Flow": y_points[i, j],
                     "GAN_Flow":gan_flow}, with_modifier=True)
                z_points[i, j] = output["obj"]/problem_description.scaling_factors["obj"]+\
                                 global_parameter['sigma']*(\
                    max(output["feed_con1"],0) / problem_description.scaling_factors["feed_con1"]+\
                    max(output["feed_con2"],0) / problem_description.scaling_factors["feed_con2"]+\
                    max(output["purity_con"],0) / problem_description.scaling_factors["purity_con"]+\
                    max(output["drain_con"],0) / problem_description.scaling_factors["drain_con"]
                )
        plt.figure()
        plt.contourf(x_points, y_points, z_points, 20)
        plt.colorbar()
        plt.scatter(rto_algorithm.input_history_data[iter_no+1]['TA_Flow'], \
                    rto_algorithm.input_history_data[iter_no+1]['MA_Flow'], s=10, c='black')
        plt.savefig("pic/HPC3/HPC3_iter" + str(iter_no) + ".png")
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
    # noise_filename = "noise/hpc-noise0.txt"
    noise_filename = "noise/hpc-noise-for-paper.txt"
    nlp_solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    qcqp_solver_executable = r"D:/Softwares/IBM/ILOG/CPLEX_Studio221/cplex/bin/x64_win64/cplex.exe"
    minlp_solver_executable = r"F:\Research\GasNetwork\kazda2020\build_model\scipampl.exe"
    result_filename_folder="data/hpc3d/"
    if not os.path.exists(result_filename_folder):
        os.makedirs(result_filename_folder)

    # ------------------------------------
    filtering_factor = 0.8
    starting_point = {'TA_Flow': 300,
                      'MA_Flow': 1750,
                      }

    # ------------------------------------
    perturbation_stepsize = {'TA_Flow': 10,
                'MA_Flow': 10,
                }
    # ------------------------------------
    print_iter_data = False
    max_iter = 30

    timpstep = int(time.time())
    sys.stdout = open("report/HPC3d_report_%d.txt"%timpstep, mode='w')
    # ------------------------------------
    # print("\nTesting MA")
    # result_filename_header = result_filename_folder + "MA_"
    # do_test_MA(perturbation_stepsize, starting_point, filtering_factor, \
    #                  noise_filename, nlp_solver_executable, \
    #                  print_iter_data, max_iter, \
    #                  result_filename_header)
    # ------------------------------------
    # print("\nTesting QCPWL-MA")
    # result_filename_header = result_filename_folder + "QCPWL_MA_"
    # do_test_QCPWL_MA(perturbation_stepsize, starting_point, filtering_factor, \
    #                  noise_filename, nlp_solver_executable, qcqp_solver_executable, \
    #                  print_iter_data, max_iter, \
    #                  result_filename_header)
    # ------------------------------------
    # print("\nTesting QCPWL-Subgrad")
    # result_filename_header = result_filename_folder + "QCPWL_Subgrad_"
    # do_test_QCPWL_Subgrad(perturbation_stepsize, starting_point, filtering_factor, \
    #                                 noise_filename, nlp_solver_executable, qcqp_solver_executable, \
    #                                 print_iter_data, max_iter, \
    #                                 result_filename_header)
    # ------------------------------------
    # print("\nTesting QCPWL-Subgrad-TR")
    # result_filename_header = result_filename_folder + "QCPWL_Subgrad_PenaltyTR_"
    # do_test_QCPWL_Subgrad_PenaltyTR(perturbation_stepsize, starting_point, filtering_factor, \
    #                  noise_filename, nlp_solver_executable, qcqp_solver_executable, \
    #                  print_iter_data, max_iter, \
    #                  result_filename_header)
    # ------------------------------------
    print("\nTesting QCPWL-Subgrad-MINLP")
    result_filename_header = result_filename_folder + "QCPWL_Subgrad_MINLP_"
    do_test_QCPWL_Subgrad_MINLP(perturbation_stepsize, starting_point, filtering_factor, \
                                    noise_filename, nlp_solver_executable, minlp_solver_executable, \
                                    print_iter_data, max_iter, \
                                    result_filename_header)
    # ------------------------------------
    

def compare_MA_and_QCPWL(pic_name):
    max_iter_no=30
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

    # algorighms = ['QCPWL_MA', 'MA','QCPWL_Subgrad','QCPWL_Subgrad_MINLP']
    algorighms = ['MA','QCPWL_Subgrad','QCPWL_Subgrad_MINLP']
    algo_name_dict = {
        # 'QCPWL_MA':'QCPWL-MA',
        'MA':"MA",
        'QCPWL_Subgrad':'QCPWL_Subgrad',
        'QCPWL_Subgrad_MINLP':'QCPWL_Subgrad_MINLP',
    }
    algo_color = {
        # 'QCPWL_MA':'red',
        'MA':"blue",
        'QCPWL_Subgrad':"black",
        'QCPWL_Subgrad_MINLP':"magenta"
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
            model_data = pandas.read_csv("data/hpc3d/" + algo + '_model_data' + ".txt", sep='\t', index_col=0, header=0)
            plant_data = pandas.read_csv("data/hpc3d/" + algo + '_plant_data' + ".txt", sep='\t', index_col=0, header=0)
            time = model_data.index
            plt.plot(time[0:max_iter_no], plant_data.iloc[0:(max_iter_no),:].loc[:,op], linewidth=global_linewidth, label=algo,
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
            model_data = pandas.read_csv("data/hpc3d/" + algo + '_model_data' + ".txt", sep='\t', index_col=0, header=0)
            input_data = pandas.read_csv("data/hpc3d/" + algo + '_input_data' + ".txt", sep='\t', index_col=0, header=0)
            time = model_data.index
            plt.plot(time[0:max_iter_no], input_data.iloc[1:(max_iter_no+1),:].loc[:,op], linewidth=global_linewidth, label=algo,
                                       color=algo_color[algo],
                                       linestyle='-')
        if op_index == 0:
            plt.legend(prop=font_legend)
        plt.ylabel(mv_label_dict[op], font2)
        if op_index == 4:
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
    # get_plant_optimum()
    # plant_simulator_test()
    # generate_noise_file()
    # draw_plant_model_mismatch()
    algo_test_main()
    # get_plant_contour_data()
    # draw_plant_model_mismatch2()
    compare_MA_and_QCPWL("pic/HPC3/compare_MA_and_QCPWL_noised.png")
