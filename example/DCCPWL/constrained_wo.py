import numpy as np
import pickle
from rtolib.core.dc_cpwl_model import QuadraticBoostedDCCPWLFunction
from rtolib.core.algo.dccpwl_ma import DCCPWL_ModifierAdaptation
from rtolib.model.wo_reactor_con import RTO_Plant_WO_reactor, default_WOR_description
from rtolib.core import NoiseGenerator, ModifierType,\
                    SimpleFiniteDiffPerturbation
import os
import copy
from rtolib.util.misc import save_iteration_data_in_dict

class ConstrainedWO_obj_CPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "model/Constrained_WO_3d_processed_output_0.bin"
        self.load_from_cpwl_model_file(cpwl_file)
        self.input_variables = ['Tr', 'Fa', 'Fb']
        self.parameters = []

class ConstrainedWO_obj_QCPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "qmodel/Constrained_WO_3d_cpwl_processed_output_0.bin"
        quadr_file = "qmodel/Constrained_WO_3d_quadr_processed_output_0.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['Tr', 'Fa', 'Fb']
        self.parameters = []

class ConstrainedWO_con_CPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "model/Constrained_WO_3d_processed_output_1.bin"
        self.load_from_cpwl_model_file(cpwl_file)
        self.input_variables = ['Tr', 'Fa', 'Fb']
        self.parameters = []

class ConstrainedWO_con_QCPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "qmodel/Constrained_WO_3d_cpwl_processed_output_1.bin"
        quadr_file = "qmodel/Constrained_WO_3d_quadr_processed_output_1.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['Tr', 'Fa', 'Fb']
        self.parameters = []

def generate_noise_file():
    noise_level = {
            'cost': 0,
            'con': 0,
    }
    max_iter = 200
    noise_generator = NoiseGenerator(noise_level)
    noise_generator.generate_noise(max_iter, 5)
    noise_generator.save_noise("noise/con-wo-noise0.txt")

def do_test_CPWL_MA(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, nlp_solver_executable, qcqp_solver_executable,\
                     print_iter_data, max_iter,\
               result_filename_header):
    problem_description = copy.deepcopy(default_WOR_description)
    problem_description.symbol_list['CV']=('cost','con')

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
        plant=RTO_Plant_WO_reactor(),
        model_dc_cpwl_functions={"cost":ConstrainedWO_obj_CPWL(),
                                 "con": ConstrainedWO_con_CPWL()},
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        nlp_solver_executable=nlp_solver_executable,
        qcqp_solver_executable=qcqp_solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        parameter_set=[],
        )

    homotopy_var=['Fa','Fb','Tr']
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
    problem_description.symbol_list['CV']=('cost','con')

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
        plant=RTO_Plant_WO_reactor(),
        model_dc_cpwl_functions={"cost":ConstrainedWO_obj_QCPWL(),
                                 "con": ConstrainedWO_con_QCPWL()},
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        nlp_solver_executable=nlp_solver_executable,
        qcqp_solver_executable=qcqp_solver_executable,
        spec_function=None,
        modifier_type=ModifierType.RTO,
        parameter_set=[],
        )

    homotopy_var=['Fa','Fb','Tr']
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

def algo_test_main():
    # ------------------------------------
    noise_filename = "noise/con-wo-noise0.txt"
    nlp_solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    qcqp_solver_executable = None
    result_filename_folder="data/con-wo/"
    if not os.path.exists(result_filename_folder):
        os.makedirs(result_filename_folder)

    # ------------------------------------
    filtering_factor = 0.1
    starting_point = {
        "Fa": 3.6,
        "Fb": 10,
        "Tr": 85,
    }

    # ------------------------------------
    perturbation_stepsize = {
        "Fa": 0.02,
        "Fb": 0.02,
        "Tr": 0.2,
    }
    # ------------------------------------
    print_iter_data = False
    max_iter = 10

    # ------------------------------------
    print("\nTesting CPWL-MA")
    result_filename_header = result_filename_folder + "CPWL_MA_"
    do_test_CPWL_MA(perturbation_stepsize, starting_point, filtering_factor, \
               noise_filename, nlp_solver_executable, qcqp_solver_executable,\
               print_iter_data, max_iter, \
               result_filename_header)

    # print("\nTesting QCPWL-MA")
    # result_filename_header = result_filename_folder + "QCPWL_MA_"
    # do_test_QCPWL_MA(perturbation_stepsize, starting_point, filtering_factor, \
    #                  noise_filename, nlp_solver_executable, qcqp_solver_executable, \
    #                  print_iter_data, max_iter, \
    #                  result_filename_header)

if __name__ == "__main__":
    # generate_noise_file()
    algo_test_main()