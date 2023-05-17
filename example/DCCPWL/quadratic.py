import numpy as np
import pickle
from rtolib.core.dc_cpwl_model import QuadraticBoostedDCCPWLFunction, \
    QuadraticBoostedDCCPWLFunctionSubgrad
from rtolib.core.cpwl_minlp import QCPWLFunction_MINLP
from rtolib.core.algo.dccpwl_ma import DCCPWL_ModifierAdaptation
from rtolib.core.algo.subgrad import DCCPWL_ModifierAdaptationSubgrad,\
        DCCPWL_ModifierAdaptationTRPenalty, QCPWL_Subgrad_MINLP
from rtolib.model.qcpwl import constrained_quadratic_problem_description,\
    unconstrained_quadratic_problem_description,unconstrained_quadratic_problem_plant,\
    constrained_quadratic_problem_plant
from rtolib.core import NoiseGenerator, ModifierType,\
                    SimpleFiniteDiffPerturbation
from rtolib.core.solve import PyomoSimulator
import os
import copy
from rtolib.util.misc import save_iteration_data_in_dict
import matplotlib.pyplot as plt
from pyomo.environ import value

global_parameter={
        "eta1": 0.01,
        "eta2": 0.9,
        "gamma1": 0.5,
        "gamma2": 1,
        "gamma3": 2,
        "sigma": 10,
        "max_trust_radius": 2,
        "initial_trust_radius": 1,
}

class cost_QCPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        self.input_variables = ['u1']
        self.parameters = []
        self.dimension = 1
        self.no_pos_segment = 2
        self.no_neg_segment = 3
        self.pos_cpwl_coeff = np.array([[-1,0],[1,0]])
        self.neg_cpwl_coeff = np.array([[0,0],[-0.5,-0.5],[0.5,-0.5]])
        # self.no_pos_segment = 1
        # self.no_neg_segment = 1
        # self.pos_cpwl_coeff = np.array([[0, 0]])
        # self.neg_cpwl_coeff = np.array([[0, 0]])
        self.pos_quadratic_A = np.array([[1]])
        # self.pos_quadratic_A = np.array([[0]])
        self.pos_quadratic_b = np.array([0])
        self.pos_quadratic_c = 0
        self.neg_quadratic_A = np.array([[0]])
        self.neg_quadratic_b = np.zeros((self.dimension,))
        self.neg_quadratic_c = 0

class con_QCPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        self.input_variables = ['u1']
        self.parameters = []
        self.dimension = 1
        self.no_pos_segment = 2
        self.no_neg_segment = 3
        self.pos_cpwl_coeff = np.array([[-1,1],[1,-1]])
        self.neg_cpwl_coeff = np.array([[0,0],[-0.5,0],[0.5,-1]])
        self.pos_quadratic_A = np.array([[1]])
        self.pos_quadratic_b = np.array([-2])
        # self.pos_quadratic_A = np.array([[0]])
        # self.pos_quadratic_b = np.array([0])
        self.pos_quadratic_c = -9
        self.neg_quadratic_A = np.array([[0]])
        self.neg_quadratic_b = np.zeros((self.dimension,))
        self.neg_quadratic_c = 0

class validity_con_QCPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        self.input_variables = ['u1']
        self.parameters = []
        self.dimension = 1
        self.no_pos_segment = 1
        self.no_neg_segment = 1
        self.pos_cpwl_coeff = np.array([[0,0]])
        self.neg_cpwl_coeff = np.array([[0,0]])
        self.pos_quadratic_A = np.array([[0]])
        self.pos_quadratic_b = np.array([0])
        self.pos_quadratic_c = -1
        self.neg_quadratic_A = np.array([[0]])
        self.neg_quadratic_b = np.zeros((self.dimension,))
        self.neg_quadratic_c = 0

class cost_QCPWL_subgrad(QuadraticBoostedDCCPWLFunctionSubgrad):
    def __init__(self):
        self.input_variables = ['u1']
        self.parameters = []
        self.dimension = 1
        self.no_pos_segment = 2
        self.no_neg_segment = 3
        self.pos_cpwl_coeff = np.array([[-1,0],[1,0]])
        self.neg_cpwl_coeff = np.array([[0,0],[-0.5,-0.5],[0.5,-0.5]])
        # self.no_pos_segment = 1
        # self.no_neg_segment = 1
        # self.pos_cpwl_coeff = np.array([[0, 0]])
        # self.neg_cpwl_coeff = np.array([[0, 0]])
        self.pos_quadratic_A = np.array([[1]])
        self.pos_quadratic_b = np.array([0])
        self.pos_quadratic_c = 0
        self.neg_quadratic_A = np.array([[0]])
        self.neg_quadratic_b = np.zeros((self.dimension,))
        self.neg_quadratic_c = 0

class con_QCPWL_subgrad(QuadraticBoostedDCCPWLFunctionSubgrad):
    def __init__(self):
        self.input_variables = ['u1']
        self.parameters = []
        self.dimension = 1
        self.no_pos_segment = 2
        self.no_neg_segment = 3
        self.pos_cpwl_coeff = np.array([[-1,1],[1,-1]])
        self.neg_cpwl_coeff = np.array([[0,0],[-0.5,0],[0.5,-1]])
        self.pos_quadratic_A = np.array([[1]])
        self.pos_quadratic_b = np.array([-2])
        self.pos_quadratic_c = -9
        self.neg_quadratic_A = np.array([[0]])
        self.neg_quadratic_b = np.zeros((self.dimension,))
        self.neg_quadratic_c = 0

class validity_con_QCPWL_subgrad(QuadraticBoostedDCCPWLFunctionSubgrad):
    def __init__(self):
        self.input_variables = ['u1']
        self.parameters = []
        self.dimension = 1
        self.no_pos_segment = 1
        self.no_neg_segment = 1
        self.pos_cpwl_coeff = np.array([[0,0]])
        self.neg_cpwl_coeff = np.array([[0,0]])
        self.pos_quadratic_A = np.array([[0]])
        self.pos_quadratic_b = np.array([0])
        self.pos_quadratic_c = -1
        self.neg_quadratic_A = np.array([[0]])
        self.neg_quadratic_b = np.zeros((self.dimension,))
        self.neg_quadratic_c = 0

class cost_QCPWL_subgrad_MINLP(QCPWLFunction_MINLP):
    def __init__(self):
        self.input_variables = ['u1']
        self.parameters = []
        self.dimension = 1
        self.no_pos_segment = 2
        self.no_neg_segment = 3
        self.pos_cpwl_coeff = np.array([[-1,0],[1,0]])
        self.neg_cpwl_coeff = np.array([[0,0],[-0.5,-0.5],[0.5,-0.5]])
        # self.no_pos_segment = 1
        # self.no_neg_segment = 1
        # self.pos_cpwl_coeff = np.array([[0, 0]])
        # self.neg_cpwl_coeff = np.array([[0, 0]])
        self.pos_quadratic_A = np.array([[1]])
        self.pos_quadratic_b = np.array([0])
        self.pos_quadratic_c = 0
        self.neg_quadratic_A = np.array([[0]])
        self.neg_quadratic_b = np.zeros((self.dimension,))
        self.neg_quadratic_c = 0
        self.bigM = 1000

class con_QCPWL_subgrad_MINLP(QCPWLFunction_MINLP):
    def __init__(self):
        self.input_variables = ['u1']
        self.parameters = []
        self.dimension = 1
        self.no_pos_segment = 2
        self.no_neg_segment = 3
        self.pos_cpwl_coeff = np.array([[-1,1],[1,-1]])
        self.neg_cpwl_coeff = np.array([[0,0],[-0.5,0],[0.5,-1]])
        self.pos_quadratic_A = np.array([[1]])
        self.pos_quadratic_b = np.array([-2])
        self.pos_quadratic_c = -9
        self.neg_quadratic_A = np.array([[0]])
        self.neg_quadratic_b = np.zeros((self.dimension,))
        self.neg_quadratic_c = 0
        self.bigM = 1000

class validity_con_QCPWL_subgrad_MINLP(QCPWLFunction_MINLP):
    def __init__(self):
        self.input_variables = ['u1']
        self.parameters = []
        self.dimension = 1
        self.no_pos_segment = 1
        self.no_neg_segment = 1
        self.pos_cpwl_coeff = np.array([[0,0]])
        self.neg_cpwl_coeff = np.array([[0,0]])
        self.pos_quadratic_A = np.array([[0]])
        self.pos_quadratic_b = np.array([0])
        self.pos_quadratic_c = -1
        self.neg_quadratic_A = np.array([[0]])
        self.neg_quadratic_b = np.zeros((self.dimension,))
        self.neg_quadratic_c = 0
        self.bigM = 1000

def generate_noise_file():
    noise_level = {
        "cost": 0,
        "con": 0,
    }
    max_iter = 200
    noise_generator = NoiseGenerator(noise_level)
    noise_generator.generate_noise(max_iter, 6)
    noise_generator.save_noise("noise/quadratic0.txt")

def do_test_QCPWL_MA_uncon(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, nlp_solver_executable, qcqp_solver_executable,\
                     print_iter_data, max_iter,\
               result_filename_header):
    problem_description = copy.deepcopy(unconstrained_quadratic_problem_description)

    rto_algorithm = DCCPWL_ModifierAdaptation()
    options = {
        "homotopy_simulation": False,
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
        plant=unconstrained_quadratic_problem_plant(),
        model_dc_cpwl_functions={'cost':cost_QCPWL(),
                                 },
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        nlp_solver_executable=nlp_solver_executable,
        qcqp_solver_executable=qcqp_solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        parameter_set=[],
        )

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    num_point = 100
    x_points=np.linspace(start=-5, stop=5, num=num_point)
    plant_cost = np.zeros(x_points.shape)
    model_cost = np.zeros(x_points.shape)
    for i in range(num_point):
        y,_ = rto_algorithm.plant_simulator.simulate({"u1":x_points[i]}, param_values=None,
                                            use_homo=False)
        plant_cost[i] = y['cost']
        y,_ = rto_algorithm.DC_CPWL_RTO_model.simulate({"u1":x_points[i]}, with_modifier=True)
        model_cost[i] = y['cost']
    plt.figure()
    plt.plot(x_points, plant_cost, label="plant")
    plt.plot(x_points, model_cost, label="model")
    plt.scatter(rto_algorithm.input_history_data[0]['u1'], 0, s=5,c='black')
    plt.legend()
    plt.show()
    plt.close()


    for iter_no in range(max_iter):
        rto_algorithm.one_step_simulation()
        if print_iter_data:
            print(rto_algorithm.model_history_data)
            print(rto_algorithm.plant_history_data)
            print(rto_algorithm.input_history_data)
        for i in range(num_point):
            y, _ = rto_algorithm.DC_CPWL_RTO_model.simulate({"u1": x_points[i]}, with_modifier=True)
            model_cost[i] = y['cost']
        plt.figure()
        plt.plot(x_points, plant_cost, label="plant")
        plt.plot(x_points, model_cost, label="model")
        plt.scatter(rto_algorithm.input_history_data[iter_no]['u1'], 0, s=5,c='black')
        plt.scatter(rto_algorithm.input_history_data[iter_no+1]['u1'], 0,s=5,c='black')
        plt.legend()
        plt.show()
        plt.close()

    # save data
    save_iteration_data_in_dict(rto_algorithm.model_history_data, result_filename_header + "model_data.txt")
    save_iteration_data_in_dict(rto_algorithm.plant_history_data, result_filename_header + "plant_data.txt")
    save_iteration_data_in_dict(rto_algorithm.input_history_data, result_filename_header + "input_data.txt")

def algo_test_main_uncon():
    # ------------------------------------
    noise_filename = "noise/quadratic0.txt"
    nlp_solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    qcqp_solver_executable = None
    result_filename_folder="data/quadratic0/"
    if not os.path.exists(result_filename_folder):
        os.makedirs(result_filename_folder)

    # ------------------------------------
    filtering_factor = 1
    starting_point = {'u1': 4,
                      }

    # ------------------------------------
    perturbation_stepsize = {'u1': 0.01,
                }
    # ------------------------------------
    print_iter_data = False
    max_iter = 10

    # ------------------------------------
    print("\nTesting QCPWL-MA")
    result_filename_header = result_filename_folder + "QCPWL_MA_"
    do_test_QCPWL_MA_uncon(perturbation_stepsize, starting_point, filtering_factor, \
                     noise_filename, nlp_solver_executable, qcqp_solver_executable, \
                     print_iter_data, max_iter, \
                     result_filename_header)

def do_test_QCPWL_MA_con(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, nlp_solver_executable, qcqp_solver_executable,\
                     print_iter_data, max_iter,\
               result_filename_header):
    problem_description = copy.deepcopy(constrained_quadratic_problem_description)

    rto_algorithm = DCCPWL_ModifierAdaptation()
    options = {
        "homotopy_simulation": False,
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
        plant=constrained_quadratic_problem_plant(),
        model_dc_cpwl_functions={'cost':cost_QCPWL(),
                                 'con':con_QCPWL(),
                                 'validity_con': validity_con_QCPWL(),
                                 },
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        nlp_solver_executable=nlp_solver_executable,
        qcqp_solver_executable=qcqp_solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        parameter_set=[],
        )

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    num_point = 100
    x_points=np.linspace(start=-5, stop=5, num=num_point)
    plant_cost = np.zeros(x_points.shape)
    plant_con = np.zeros(x_points.shape)
    model_cost = np.zeros(x_points.shape)
    model_con = np.zeros(x_points.shape)
    for i in range(num_point):
        y,_ = rto_algorithm.plant_simulator.simulate({"u1":x_points[i]}, param_values=None,
                                            use_homo=False)
        plant_cost[i] = y['cost']
        plant_con[i] = y['con']
        y,_ = rto_algorithm.DC_CPWL_RTO_model.simulate({"u1":x_points[i]}, with_modifier=True)
        model_cost[i] = y['cost']
        model_con[i] = y['con']
    plt.figure()
    plt.subplot(211)
    plt.title("cost")
    plt.plot(x_points, plant_cost, label="plant")
    plt.plot(x_points, model_cost, label="model")
    plt.scatter(rto_algorithm.input_history_data[0]['u1'], 0, s=5,c='black')
    plt.legend()
    plt.subplot(212)
    plt.title("constraint")
    plt.plot(x_points, plant_con, label="plant")
    plt.plot(x_points, model_con, label="model")
    plt.scatter(rto_algorithm.input_history_data[0]['u1'], 0, s=5, c='black')
    plt.legend()
    plt.show()
    plt.close()


    for iter_no in range(max_iter):
        rto_algorithm.one_step_simulation()
        if print_iter_data:
            print(rto_algorithm.model_history_data)
            print(rto_algorithm.plant_history_data)
            print(rto_algorithm.input_history_data)
        for i in range(num_point):
            y, _ = rto_algorithm.DC_CPWL_RTO_model.simulate({"u1": x_points[i]}, with_modifier=True)
            model_cost[i] = y['cost']
            model_con[i] = y['con']
        plt.figure()
        plt.subplot(211)
        plt.title("cost")
        plt.plot(x_points, plant_cost, label="plant")
        plt.plot(x_points, model_cost, label="model")
        plt.scatter(rto_algorithm.input_history_data[iter_no]['u1'], 0, s=5, c='black')
        plt.scatter(rto_algorithm.input_history_data[iter_no+1]['u1'], 0, s=5, c='black')
        plt.legend()
        plt.subplot(212)
        plt.title("constraint")
        plt.plot(x_points, plant_con, label="plant")
        plt.plot(x_points, model_con, label="model")
        plt.scatter(rto_algorithm.input_history_data[iter_no]['u1'], 0, s=5, c='black')
        plt.scatter(rto_algorithm.input_history_data[iter_no+1]['u1'], 0, s=5, c='black')
        plt.legend()
        plt.show()
        plt.close()

    # save data
    save_iteration_data_in_dict(rto_algorithm.model_history_data, result_filename_header + "model_data.txt")
    save_iteration_data_in_dict(rto_algorithm.plant_history_data, result_filename_header + "plant_data.txt")
    save_iteration_data_in_dict(rto_algorithm.input_history_data, result_filename_header + "input_data.txt")

def algo_test_main_con():
    # ------------------------------------
    noise_filename = "noise/quadratic0.txt"
    nlp_solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    qcqp_solver_executable = None
    result_filename_folder="data/quadratic0/"
    if not os.path.exists(result_filename_folder):
        os.makedirs(result_filename_folder)

    # ------------------------------------
    filtering_factor = 1
    # starting from 4, the input converges to a nonoptimal point
    starting_point = {'u1': -4,#4,
                      }

    # ------------------------------------
    perturbation_stepsize = {'u1': 0.0001,
                }
    # ------------------------------------
    print_iter_data = False
    max_iter = 20

    # ------------------------------------
    print("\nTesting QCPWL-MA")
    result_filename_header = result_filename_folder + "QCPWL_MA_"
    do_test_QCPWL_MA_con(perturbation_stepsize, starting_point, filtering_factor, \
                     noise_filename, nlp_solver_executable, qcqp_solver_executable, \
                     print_iter_data, max_iter, \
                     result_filename_header)

def do_test_QCPWL_MA_subgrad_con(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, nlp_solver_executable, qcqp_solver_executable,\
                     print_iter_data, max_iter,\
               result_filename_header):
    problem_description = copy.deepcopy(constrained_quadratic_problem_description)

    rto_algorithm = DCCPWL_ModifierAdaptationSubgrad()
    options = {
        "homotopy_simulation": False,
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
        plant=constrained_quadratic_problem_plant(),
        model_dc_cpwl_functions={'cost':cost_QCPWL_subgrad(),
                                 'con':con_QCPWL_subgrad(),
                                 'validity_con':validity_con_QCPWL_subgrad(),
                                 },
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        nlp_solver_executable=nlp_solver_executable,
        qcqp_solver_executable=qcqp_solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        parameter_set=[],
        )

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    num_point = 100
    x_points=np.linspace(start=-5, stop=5, num=num_point)
    plant_cost = np.zeros(x_points.shape)
    plant_con = np.zeros(x_points.shape)
    model_cost = np.zeros(x_points.shape)
    model_con = np.zeros(x_points.shape)
    for i in range(num_point):
        y,_ = rto_algorithm.plant_simulator.simulate({"u1":x_points[i]}, param_values=None,
                                            use_homo=False)
        plant_cost[i] = y['cost']
        plant_con[i] = y['con']
        y,_ = rto_algorithm.DC_CPWL_RTO_model.simulate({"u1":x_points[i]}, with_modifier=True)
        model_cost[i] = y['cost']
        model_con[i] = y['con']
    plt.figure()
    plt.subplot(211)
    plt.plot(x_points, plant_cost, label="plant")
    plt.plot(x_points, model_cost, label="model")
    plt.scatter(rto_algorithm.input_history_data[0]['u1'], 0, s=5,c='black')
    plt.legend()
    plt.subplot(212)
    plt.plot(x_points, plant_con, label="plant")
    plt.plot(x_points, model_con, label="model")
    plt.scatter(rto_algorithm.input_history_data[0]['u1'], 0, s=5, c='black')
    plt.legend()
    plt.show()
    plt.close()


    for iter_no in range(max_iter):
        rto_algorithm.one_step_simulation()
        if print_iter_data:
            print(rto_algorithm.model_history_data)
            print(rto_algorithm.plant_history_data)
            print(rto_algorithm.input_history_data)
        modifier_k = np.array([value(rto_algorithm.DC_CPWL_RTO_model.subproblem_model.output['cost'].linear_correction_k[0]),\
               value(rto_algorithm.DC_CPWL_RTO_model.subproblem_model.output['con'].linear_correction_k[0])])
        plant_k = np.array([value(rto_algorithm.DC_CPWL_RTO_model.dc_cpwl_functions['cost'].plant_k[0]), \
               value(rto_algorithm.DC_CPWL_RTO_model.dc_cpwl_functions['con'].plant_k[0])])
        concave_k = np.array([value(rto_algorithm.DC_CPWL_RTO_model.dc_cpwl_functions['cost'].concave_k[0]), \
                            value(rto_algorithm.DC_CPWL_RTO_model.dc_cpwl_functions['con'].concave_k[0])])
        convex_k = np.array([value(rto_algorithm.DC_CPWL_RTO_model.dc_cpwl_functions['cost'].convex_k[0]), \
                            value(rto_algorithm.DC_CPWL_RTO_model.dc_cpwl_functions['con'].convex_k[0])])
        base_point = np.array(
            [value(rto_algorithm.DC_CPWL_RTO_model.subproblem_model.output['cost'].base_point[0]), \
             value(rto_algorithm.DC_CPWL_RTO_model.subproblem_model.output['con'].base_point[0])])
        # print(convex_k)
        # print(concave_k)
        # print(modifier_k)
        # print(plant_k)
        print("base point")
        print(base_point)
        print("gradient matching error")
        print(convex_k-concave_k+modifier_k-plant_k)
        # rto_algorithm.DC_CPWL_RTO_model.subproblem_model.display("temp.txt")
        for i in range(num_point):
            y, _ = rto_algorithm.DC_CPWL_RTO_model.simulate({"u1": x_points[i]}, with_modifier=True)
            model_cost[i] = y['cost']
            model_con[i] = y['con']
        plt.figure()
        plt.subplot(211)
        plt.plot(x_points, plant_cost, label="plant")
        plt.plot(x_points, model_cost, label="model")
        plt.ylim([-2,15])
        plt.scatter(rto_algorithm.input_history_data[iter_no]['u1'], 0, s=5, c='black')
        plt.scatter(rto_algorithm.input_history_data[iter_no+1]['u1'], 0, s=5, c='black')
        plt.legend()
        plt.subplot(212)
        plt.plot(x_points, plant_con, label="plant")
        plt.plot(x_points, model_con, label="model")
        plt.scatter(rto_algorithm.input_history_data[iter_no]['u1'], 0, s=5, c='black')
        plt.scatter(rto_algorithm.input_history_data[iter_no+1]['u1'], 0, s=5, c='black')
        plt.legend()
        plt.show()
        plt.close()

    # save data
    save_iteration_data_in_dict(rto_algorithm.model_history_data, result_filename_header + "model_data.txt")
    save_iteration_data_in_dict(rto_algorithm.plant_history_data, result_filename_header + "plant_data.txt")
    save_iteration_data_in_dict(rto_algorithm.input_history_data, result_filename_header + "input_data.txt")

def do_test_QCPWL_MA_subgrad_con_penalty_tr(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, nlp_solver_executable, qcqp_solver_executable,\
                     print_iter_data, max_iter,\
               result_filename_header):
    problem_description = copy.deepcopy(constrained_quadratic_problem_description)

    rto_algorithm = DCCPWL_ModifierAdaptationTRPenalty()
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
        plant=constrained_quadratic_problem_plant(),
        model_dc_cpwl_functions={'cost':cost_QCPWL_subgrad(),
                                 'con':con_QCPWL_subgrad(),
                                 'validity_con':validity_con_QCPWL_subgrad(),
                                 },
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        nlp_solver_executable=nlp_solver_executable,
        qcqp_solver_executable=qcqp_solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        parameter_set=[],
        )

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    num_point = 100
    x_points=np.linspace(start=-5, stop=5, num=num_point)
    plant_cost = np.zeros(x_points.shape)
    plant_con = np.zeros(x_points.shape)
    model_cost = np.zeros(x_points.shape)
    model_con = np.zeros(x_points.shape)
    for i in range(num_point):
        y,_ = rto_algorithm.plant_simulator.simulate({"u1":x_points[i]}, param_values=None,
                                            use_homo=False)
        plant_cost[i] = y['cost']
        plant_con[i] = y['con']
        y,_ = rto_algorithm.DC_CPWL_RTO_model.simulate({"u1":x_points[i]}, with_modifier=True)
        model_cost[i] = y['cost']
        model_con[i] = y['con']
    plt.figure()
    plt.subplot(211)
    plt.plot(x_points, plant_cost, label="plant")
    plt.plot(x_points, model_cost, label="model")
    plt.scatter(rto_algorithm.input_history_data[0]['u1'], 0, s=5,c='black')
    plt.legend()
    plt.subplot(212)
    plt.plot(x_points, plant_con, label="plant")
    plt.plot(x_points, model_con, label="model")
    plt.scatter(rto_algorithm.input_history_data[0]['u1'], 0, s=5, c='black')
    plt.legend()
    plt.show()
    plt.close()

    for iter_no in range(max_iter):
        rto_algorithm.one_step_simulation()
        if print_iter_data:
            print(rto_algorithm.model_history_data)
            print(rto_algorithm.plant_history_data)
            print(rto_algorithm.input_history_data)

        # rto_algorithm.DC_CPWL_RTO_model.subproblem_model.display("temp.txt")
        for i in range(num_point):
            y, _ = rto_algorithm.DC_CPWL_RTO_model.simulate({"u1": x_points[i]}, with_modifier=True)
            model_cost[i] = y['cost']
            model_con[i] = y['con']
        plt.figure()
        plt.subplot(211)
        plt.plot(x_points, plant_cost, label="plant")
        plt.plot(x_points, model_cost, label="model")
        plt.ylim([-2,15])
        plt.scatter(rto_algorithm.input_history_data[iter_no]['u1'], 0, s=5, c='black')
        plt.scatter(rto_algorithm.input_history_data[iter_no+1]['u1'], 0, s=5, c='black')
        plt.legend()
        plt.subplot(212)
        plt.plot(x_points, plant_con, label="plant")
        plt.plot(x_points, model_con, label="model")
        plt.scatter(rto_algorithm.input_history_data[iter_no]['u1'], 0, s=5, c='black')
        plt.scatter(rto_algorithm.input_history_data[iter_no+1]['u1'], 0, s=5, c='black')
        plt.legend()
        plt.show()
        plt.close()

    # save data
    save_iteration_data_in_dict(rto_algorithm.model_history_data, result_filename_header + "model_data.txt")
    save_iteration_data_in_dict(rto_algorithm.plant_history_data, result_filename_header + "plant_data.txt")
    save_iteration_data_in_dict(rto_algorithm.input_history_data, result_filename_header + "input_data.txt")

def do_test_QCPWL_MA_subgrad_MINLP_con(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, nlp_solver_executable, minlp_solver_executable,\
                     print_iter_data, max_iter,\
               result_filename_header):
    problem_description = copy.deepcopy(constrained_quadratic_problem_description)

    rto_algorithm = QCPWL_Subgrad_MINLP()
    options = {
        "homotopy_simulation": False,
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
        plant=constrained_quadratic_problem_plant(),
        model_dc_cpwl_functions={'cost':cost_QCPWL_subgrad_MINLP(),
                                 'con':con_QCPWL_subgrad_MINLP(),
                                 'validity_con':validity_con_QCPWL_subgrad_MINLP(),
                                 },
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        nlp_solver_executable=nlp_solver_executable,
        minlp_solver_executable=minlp_solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        parameter_set=[],
        )

    rto_algorithm.initialize_simulation(starting_point, initial_parameter_value={})

    num_point = 100
    x_points=np.linspace(start=-5, stop=5, num=num_point)
    plant_cost = np.zeros(x_points.shape)
    plant_con = np.zeros(x_points.shape)
    model_cost = np.zeros(x_points.shape)
    model_con = np.zeros(x_points.shape)
    for i in range(num_point):
        y,_ = rto_algorithm.plant_simulator.simulate({"u1":x_points[i]}, param_values=None,
                                            use_homo=False)
        plant_cost[i] = y['cost']
        plant_con[i] = y['con']
        y,_ = rto_algorithm.DC_CPWL_RTO_model.simulate({"u1":x_points[i]}, with_modifier=True)
        model_cost[i] = y['cost']
        model_con[i] = y['con']
    plt.figure()
    plt.subplot(211)
    plt.plot(x_points, plant_cost, label="plant")
    plt.plot(x_points, model_cost, label="model")
    plt.scatter(rto_algorithm.input_history_data[0]['u1'], 0, s=5,c='black')
    plt.legend()
    plt.subplot(212)
    plt.plot(x_points, plant_con, label="plant")
    plt.plot(x_points, model_con, label="model")
    plt.scatter(rto_algorithm.input_history_data[0]['u1'], 0, s=5, c='black')
    plt.legend()
    plt.show()
    plt.close()


    for iter_no in range(max_iter):
        rto_algorithm.one_step_simulation()
        if print_iter_data:
            print(rto_algorithm.model_history_data)
            print(rto_algorithm.plant_history_data)
            print(rto_algorithm.input_history_data)
        modifier_k = np.array([value(rto_algorithm.DC_CPWL_RTO_model.subproblem_model.output['cost'].linear_correction_k[0]),\
               value(rto_algorithm.DC_CPWL_RTO_model.subproblem_model.output['con'].linear_correction_k[0])])
        plant_k = np.array([value(rto_algorithm.DC_CPWL_RTO_model.dc_cpwl_functions['cost'].plant_k[0]), \
               value(rto_algorithm.DC_CPWL_RTO_model.dc_cpwl_functions['con'].plant_k[0])])
        concave_k = np.array([value(rto_algorithm.DC_CPWL_RTO_model.dc_cpwl_functions['cost'].concave_k[0]), \
                            value(rto_algorithm.DC_CPWL_RTO_model.dc_cpwl_functions['con'].concave_k[0])])
        convex_k = np.array([value(rto_algorithm.DC_CPWL_RTO_model.dc_cpwl_functions['cost'].convex_k[0]), \
                            value(rto_algorithm.DC_CPWL_RTO_model.dc_cpwl_functions['con'].convex_k[0])])
        base_point = np.array(
            [value(rto_algorithm.DC_CPWL_RTO_model.subproblem_model.output['cost'].base_point[0]), \
             value(rto_algorithm.DC_CPWL_RTO_model.subproblem_model.output['con'].base_point[0])])
        # print(convex_k)
        # print(concave_k)
        # print(modifier_k)
        # print(plant_k)
        print("base point")
        print(base_point)
        print("gradient matching error")
        print(convex_k-concave_k+modifier_k-plant_k)
        # rto_algorithm.DC_CPWL_RTO_model.subproblem_model.display("temp.txt")
        for i in range(num_point):
            y, _ = rto_algorithm.DC_CPWL_RTO_model.simulate({"u1": x_points[i]}, with_modifier=True)
            model_cost[i] = y['cost']
            model_con[i] = y['con']
        plt.figure()
        plt.subplot(211)
        plt.plot(x_points, plant_cost, label="plant")
        plt.plot(x_points, model_cost, label="model")
        plt.ylim([-2,15])
        plt.scatter(rto_algorithm.input_history_data[iter_no]['u1'], 0, s=5, c='black')
        plt.scatter(rto_algorithm.input_history_data[iter_no+1]['u1'], 0, s=5, c='black')
        plt.legend()
        plt.subplot(212)
        plt.plot(x_points, plant_con, label="plant")
        plt.plot(x_points, model_con, label="model")
        plt.scatter(rto_algorithm.input_history_data[iter_no]['u1'], 0, s=5, c='black')
        plt.scatter(rto_algorithm.input_history_data[iter_no+1]['u1'], 0, s=5, c='black')
        plt.legend()
        plt.show()
        plt.close()

    # save data
    save_iteration_data_in_dict(rto_algorithm.model_history_data, result_filename_header + "model_data.txt")
    save_iteration_data_in_dict(rto_algorithm.plant_history_data, result_filename_header + "plant_data.txt")
    save_iteration_data_in_dict(rto_algorithm.input_history_data, result_filename_header + "input_data.txt")


def algo_test_subgrad_main_con():
    # ------------------------------------
    noise_filename = "noise/quadratic0.txt"
    nlp_solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    qcqp_solver_executable = r"D:/Softwares/IBM/ILOG/CPLEX_Studio221/cplex/bin/x64_win64/cplex.exe"
    # qcqp_solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    minlp_solver_executable = r"F:\Research\GasNetwork\kazda2020\build_model\scipampl.exe"

    result_filename_folder="data/quadratic0/"
    if not os.path.exists(result_filename_folder):
        os.makedirs(result_filename_folder)

    # ------------------------------------
    filtering_factor = 1
    # starting from 4, the input converges to a nonoptimal point
    starting_point = {'u1': 4,#-4,
                      }

    # ------------------------------------
    perturbation_stepsize = {'u1': 0.00001,
                }
    # ------------------------------------
    print_iter_data = False
    max_iter = 20

    # ------------------------------------
    # print("\nTesting QCPWL-MA")
    # result_filename_header = result_filename_folder + "QCPWL_MA_"
    # do_test_QCPWL_MA_con(perturbation_stepsize, starting_point, filtering_factor, \
    #                  noise_filename, nlp_solver_executable, qcqp_solver_executable, \
    #                  print_iter_data, max_iter, \
    #                  result_filename_header)
    # ------------------------------------
    print("\nTesting QCPWL-MA-Subgrad")
    result_filename_header = result_filename_folder + "QCPWL_MA_Subgrad_"
    do_test_QCPWL_MA_subgrad_con(perturbation_stepsize, starting_point, filtering_factor, \
                     noise_filename, nlp_solver_executable, qcqp_solver_executable, \
                     print_iter_data, max_iter, \
                     result_filename_header)

    # ------------------------------------
    # print("\nTesting QCPWL-MA-Penalty-TR")
    # result_filename_header = result_filename_folder + "QCPWL_MA_PenaltyTR_"
    # do_test_QCPWL_MA_subgrad_con_penalty_tr(perturbation_stepsize, starting_point, filtering_factor, \
    #                              noise_filename, nlp_solver_executable, qcqp_solver_executable, \
    #                              print_iter_data, max_iter, \
    #                              result_filename_header)
    # ------------------------------------
    print("\nTesting QCPWL-subgrad-MINLP")
    result_filename_header = result_filename_folder + "QCPWL_MINLP_"
    do_test_QCPWL_MA_subgrad_MINLP_con(perturbation_stepsize, starting_point, filtering_factor, \
                                            noise_filename, nlp_solver_executable, minlp_solver_executable, \
                                            print_iter_data, max_iter, \
                                            result_filename_header)


if __name__ == "__main__":
    # generate_noise_file()
    # algo_test_main_uncon()
    # algo_test_main_con()
    algo_test_subgrad_main_con()