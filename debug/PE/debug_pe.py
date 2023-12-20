from rtolib.core import SimpleFiniteDiffPerturbation, \
    PyomoSimulator, \
    PyomoGradientParamEstimator,\
    NoiseGenerator
from rtolib.model.wo_reactor import RTO_Plant_WO_reactor, \
    RTO_Mismatched_WO_reactor, \
    default_WOR_description
import copy
import numpy
from pyomo.environ import SolverFactory

def model_pe_model():
    problem_description = copy.deepcopy(default_WOR_description)
    pertubation_stepsize = {
        "Fb": 0.2,
        "Tr": 2,
    }
    ffd_perturb = SimpleFiniteDiffPerturbation(pertubation_stepsize, problem_description)
    model_simulator = PyomoSimulator(RTO_Mismatched_WO_reactor())
    pe_estimator = PyomoGradientParamEstimator(RTO_Mismatched_WO_reactor(), ffd_perturb.number_of_data_points)

    model_simulator.build(problem_description)
    pe_estimator.build(problem_description)
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    solver1 = SolverFactory('ipopt', executable=solver_executable)
    sim_options = {'max_iter': 100}
    model_simulator.set_solver(solver1, tee=False, default_options=sim_options)
    solver2 = SolverFactory('ipopt', executable=solver_executable)
    pe_options = {'max_iter': 100,
                  'obj_scaling_factor':1}
    pe_estimator.set_solver(solver2, tee=False, default_options=pe_options)

    output_weight = {
        "XFr_A": 1,
        "XFr_B": 1,
        "XFr_E": 1,
        "XFr_P": 1,
        "XFr_G": 1,
    }
    parameter_weight = {
        'Ka1': 0,
        'Ka2': 0,
        'Kb1': 0,
        'Kb2': 0,
    }
    pe_estimator.set_weight(output_weight, parameter_weight)
    default_para = {
        'Ka1': 2.189e8,
        'Ka2': 4.310e13,
        'Kb1': 8077.6,
        'Kb2': 12438,
    }
    input_values = {
        "Fb": 5,
        "Tr": 80,
    }
    noise_generator = NoiseGenerator()
    noise_generator.load_noise("noise/noise1.txt")
    no = 0
    with open('result.txt', 'w') as fp:
        fp.write("Ka1_True\tKa2_True\tKb1_True\tKb2_True\tKa1_est\tKa2_est\tKb1_est\tKb2_est\n")
        for p1 in numpy.linspace(start=0.9, stop=1.1, num=3):
            for p2 in numpy.linspace(start=0.9, stop=1.1, num=3):
                for p3 in numpy.linspace(start=0.9, stop=1.1, num=3):
                    for p4 in numpy.linspace(start=0.9, stop=1.1, num=3):
                        print(no)
                        no+=1

                        para = {
                            'Ka1': default_para['Ka1'] * p1,
                            'Ka2': default_para['Ka2'] * p2,
                            'Kb1': default_para['Kb1'] * p3,
                            'Kb2': default_para['Kb2'] * p4,
                        }
                        fp.write("%.6e\t%.6e\t%.6e\t%.6e\t" % (para['Ka1'], para['Ka2'], para['Kb1'], para['Kb2']))
                        trial_points = ffd_perturb.get_trial_points(input_values)
                        plant_data = [{} for i in range(len(trial_points))]
                        try:
                            for i, p in enumerate(trial_points):
                                for k, v in p.items():
                                    plant_data[i][k] = v
                                outputs, solve_status = model_simulator.simulate(p, param_values=para, use_homo=True)

                                for k in problem_description.symbol_list['CV']:
                                    plant_data[i][k] = outputs[k]
                        except Exception as e:
                            print(e)
                            fp.write("sim_error\n")
                            continue

                        # add noise to the plant data
                        for trial_point_no in range(len(trial_points)):
                            for k in problem_description.symbol_list['CV']:
                                noise = noise_generator.get_noise(0, trial_point_no, k)
                                plant_data[trial_point_no][k] += noise

                        try:
                            para_ret = pe_estimator.estimate_parameter(plant_data, fixed_param_values={}, use_homo=True,
                                                                       pre_simulation=True)
                        except Exception as e:
                            print(e)
                            fp.write("opt_error\n")
                            continue
                        fp.write("%.6e\t%.6e\t%.6e\t%.6e\n" % (
                            para_ret['Ka1'], para_ret['Ka2'], para_ret['Kb1'], para_ret['Kb2']))


def model_pe_plant():
    problem_description = copy.deepcopy(default_WOR_description)
    pertubation_stepsize = {
        "Fb": 0.2,
        "Tr": 2,
    }
    ffd_perturb = SimpleFiniteDiffPerturbation(pertubation_stepsize, problem_description)
    plant_simulator = PyomoSimulator(RTO_Plant_WO_reactor())
    pe_estimator = PyomoGradientParamEstimator(RTO_Mismatched_WO_reactor(), ffd_perturb.number_of_data_points)

    plant_simulator.build(problem_description)
    pe_estimator.build(problem_description)
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    solver1 = SolverFactory('ipopt', executable=solver_executable)
    sim_options = {'max_iter': 100}
    plant_simulator.set_solver(solver1, tee=False, default_options=sim_options)
    solver2 = SolverFactory('ipopt', executable=solver_executable)
    pe_options = {'max_iter': 100,
                  'obj_scaling_factor': 1}
    pe_estimator.set_solver(solver2, tee=False, default_options=pe_options)

    output_weight = {
        "XFr_A": 1,
        "XFr_B": 1,
        "XFr_E": 1,
        "XFr_P": 1,
        "XFr_G": 1,
    }
    parameter_weight = {
        'Ka1': 0,
        'Ka2': 0,
        'Kb1': 0,
        'Kb2': 0,
    }
    pe_estimator.set_weight(output_weight, parameter_weight)
    default_para = {
        'Ka1': 2.189e8,
        'Ka2': 4.310e13,
        'Kb1': 8077.6,
        'Kb2': 12438,
    }
    input_values = {
        "Fb": 5,
        "Tr": 80,
    }
    noise_generator = NoiseGenerator()
    noise_generator.load_noise("noise/noise1.txt")
    no = 0
    with open('result_plant.txt', 'w') as fp:
        fp.write("Tr\tFb\tKa1_est\tKa2_est\tKb1_est\tKb2_est\n")
        for Tr in numpy.linspace(start=70, stop=100, num=9):
            for Fb in numpy.linspace(start=3, stop=6, num=9):
                    print(no)
                    no += 1
                    fp.write("%.6e\t%.6e\t"%(Tr,Fb))

                    input_values={
                        'Tr':Tr,
                        'Fb':Fb,
                    }
                    trial_points = ffd_perturb.get_trial_points(input_values)
                    plant_data = [{} for i in range(len(trial_points))]
                    try:
                        for i, p in enumerate(trial_points):
                            for k, v in p.items():
                                plant_data[i][k] = v
                            outputs, solve_status = plant_simulator.simulate(p, param_values={}, use_homo=True)

                            for k in problem_description.symbol_list['CV']:
                                plant_data[i][k] = outputs[k]
                    except Exception as e:
                        print(e)
                        fp.write("sim_error\n")
                        continue

                    # add noise to the plant data
                    for trial_point_no in range(len(trial_points)):
                        for k in problem_description.symbol_list['CV']:
                            noise = noise_generator.get_noise(0, trial_point_no, k)
                            plant_data[trial_point_no][k] += noise

                    try:
                        para_ret = pe_estimator.estimate_parameter(plant_data, fixed_param_values={}, use_homo=True,
                                                                   pre_simulation=True)
                    except Exception as e:
                        print(e)
                        fp.write("opt_error\n")
                        continue
                    fp.write("%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n" % (
                        para_ret['Ka1'], para_ret['Ka2'], para_ret['Kb1'], para_ret['Kb2'], para_ret['$residue']))


def model_pe_plant_with_profit():
    problem_description = copy.deepcopy(default_WOR_description)
    pertubation_stepsize = {
        "Fb": 0.2,
        "Tr": 2,
    }
    ffd_perturb = SimpleFiniteDiffPerturbation(pertubation_stepsize, problem_description)
    plant_simulator = PyomoSimulator(RTO_Plant_WO_reactor())
    pe_estimator = PyomoGradientParamEstimator(RTO_Mismatched_WO_reactor(), ffd_perturb.number_of_data_points)

    plant_simulator.build(problem_description)
    pe_estimator.build(problem_description)
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    solver1 = SolverFactory('ipopt', executable=solver_executable)
    sim_options = {'max_iter': 100}
    plant_simulator.set_solver(solver1, tee=False, default_options=sim_options)
    solver2 = SolverFactory('ipopt', executable=solver_executable)
    pe_options = {'max_iter': 100,
                  'obj_scaling_factor': 1}
    pe_estimator.set_solver(solver2, tee=False, default_options=pe_options)

    output_weight = {
        "XFr_A": 1,
        "XFr_B": 1,
        "XFr_E": 1,
        "XFr_P": 1,
        "XFr_G": 1,
        'profit': 1e-4,
    }
    parameter_weight = {
        'Ka1': 0,
        'Ka2': 0,
        'Kb1': 0,
        'Kb2': 0,
    }
    pe_estimator.set_weight(output_weight, parameter_weight)
    default_para = {
        'Ka1': 2.189e8,
        'Ka2': 4.310e13,
        'Kb1': 8077.6,
        'Kb2': 12438,
    }
    input_values = {
        "Fb": 5,
        "Tr": 80,
    }
    noise_generator = NoiseGenerator()
    noise_generator.load_noise("noise/noise1.txt")
    no = 0
    with open('result_plant.txt', 'w') as fp:
        fp.write("Tr\tFb\tKa1_est\tKa2_est\tKb1_est\tKb2_est\n")
        for Tr in numpy.linspace(start=70, stop=100, num=9):
            for Fb in numpy.linspace(start=3, stop=6, num=9):
                    print(no)
                    no += 1
                    fp.write("%.6e\t%.6e\t"%(Tr,Fb))

                    input_values={
                        'Tr':Tr,
                        'Fb':Fb,
                    }
                    trial_points = ffd_perturb.get_trial_points(input_values)
                    plant_data = [{} for i in range(len(trial_points))]
                    try:
                        for i, p in enumerate(trial_points):
                            for k, v in p.items():
                                plant_data[i][k] = v
                            outputs, solve_status = plant_simulator.simulate(p, param_values={}, use_homo=True)

                            for k in problem_description.symbol_list['CV']:
                                plant_data[i][k] = outputs[k]
                    except Exception as e:
                        print(e)
                        fp.write("sim_error\n")
                        continue

                    # add noise to the plant data
                    for trial_point_no in range(len(trial_points)):
                        for k in problem_description.symbol_list['CV']:
                            noise = noise_generator.get_noise(0, trial_point_no, k)
                            plant_data[trial_point_no][k] += noise

                    try:
                        para_ret = pe_estimator.estimate_parameter(plant_data, fixed_param_values={}, use_homo=True,
                                                                   pre_simulation=True)
                    except Exception as e:
                        print(e)
                        fp.write("opt_error\n")
                        continue
                    fp.write("%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n" % (
                        para_ret['Ka1'], para_ret['Ka2'], para_ret['Kb1'], para_ret['Kb2'], para_ret['$residue']))

if __name__ == "__main__":
    model_pe_plant_with_profit()