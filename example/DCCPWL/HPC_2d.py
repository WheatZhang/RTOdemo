import numpy as np
import pickle
from rtolib.core.dc_cpwl_model import QuadraticBoostedDCCPWLFunction
from rtolib.core.algo.dccpwl_ma import DCCPWL_ModifierAdaptation
from rtolib.model.HPC import RTO_Plant_HPC, default_HPC_2d_description, RTO_Mismatched_HPC
from rtolib.core import NoiseGenerator, ModifierType,\
                    SimpleFiniteDiffPerturbation
from rtolib.core.solve import PyomoSimulator
import os
import copy
from rtolib.util.misc import save_iteration_data_in_dict
from rtolib.core.pyomo_model import SolverFactory

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

class ConstrainedWO_cost_QCPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "qmodel/HPC_2d_cpwl_processed_output_0.bin"
        quadr_file = "qmodel/HPC_2d_quadr_processed_output_0.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow']
        self.parameters = []

class ConstrainedWO_feed_con1_QCPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "qmodel/HPC_2d_cpwl_processed_output_1.bin"
        quadr_file = "qmodel/HPC_2d_quadr_processed_output_1.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow']
        self.parameters = []

class ConstrainedWO_feed_con2_QCPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "qmodel/HPC_2d_cpwl_processed_output_2.bin"
        quadr_file = "qmodel/HPC_2d_quadr_processed_output_2.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow']
        self.parameters = []

class ConstrainedWO_purity_con_QCPWL(QuadraticBoostedDCCPWLFunction):
    def __init__(self):
        cpwl_file = "qmodel/HPC_2d_cpwl_processed_output_3.bin"
        quadr_file = "qmodel/HPC_2d_quadr_processed_output_3.bin"
        self.load_from_qcpwl_model_file(cpwl_file, quadr_file)
        self.input_variables = ['TA_Flow', 'MA_Flow']
        self.parameters = []

class ConstrainedWO_drain_con_QCPWL(QuadraticBoostedDCCPWLFunction):
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
        model_dc_cpwl_functions={'obj':ConstrainedWO_cost_QCPWL(),
                                "feed_con1":ConstrainedWO_feed_con1_QCPWL(),
                                "feed_con2":ConstrainedWO_feed_con2_QCPWL(),
                                "purity_con":ConstrainedWO_purity_con_QCPWL(),
                                "drain_con":ConstrainedWO_drain_con_QCPWL(),
                                 },
        perturbation_method=ffd_perturb,
        noise_generator=noise_generator,
        nlp_solver_executable=nlp_solver_executable,
        qcqp_solver_executable=qcqp_solver_executable,
        spec_function=None,
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

def algo_test_main():
    # ------------------------------------
    noise_filename = "noise/hpc-noise0.txt"
    nlp_solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    qcqp_solver_executable = None
    result_filename_folder="data/hpc2d/"
    if not os.path.exists(result_filename_folder):
        os.makedirs(result_filename_folder)

    # ------------------------------------
    filtering_factor = 0.1
    starting_point = {'TA_Flow': 300,
                      'MA_Flow': 1750,
                      }

    # ------------------------------------
    perturbation_stepsize = {'TA_Flow': 20,
                'MA_Flow': 20,
                }
    # ------------------------------------
    print_iter_data = False
    max_iter = 10

    # ------------------------------------
    print("\nTesting QCPWL-MA")
    result_filename_header = result_filename_folder + "QCPWL_MA_"
    do_test_QCPWL_MA(perturbation_stepsize, starting_point, filtering_factor, \
                     noise_filename, nlp_solver_executable, qcqp_solver_executable, \
                     print_iter_data, max_iter, \
                     result_filename_header)

if __name__ == "__main__":
    # plant_simulator_test()
    # generate_noise_file()
    algo_test_main()