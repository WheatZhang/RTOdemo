from rtolib.core.solve import PyomoSimulator
from rtolib.model.wo_reactor import RTO_Plant_WO_reactor, RTO_Mismatched_WO_reactor, default_WOR_description
from pyomo.environ import SolverFactory
from draw_pic_lib import generate_contour_data
import pandas
import os
import numpy

def generate_all_contour_data(algo_results_folder,resolution):
    datafile_folder = "data/contour_data/"
    if not os.path.exists(algo_results_folder):
        os.makedirs(algo_results_folder)
    data_files = {"Plant": datafile_folder + "_Plant",
                  "Model": datafile_folder + "_Model",
                  "PE": datafile_folder + "_Plant",
                  "ISOPE": datafile_folder + "_ISOPE",
                  "GPE": datafile_folder + "_GPE", }

    plant_simulator = PyomoSimulator(RTO_Plant_WO_reactor())
    model_simulator = PyomoSimulator(RTO_Mismatched_WO_reactor())
    plant_simulator.build(default_WOR_description)
    model_simulator.build(default_WOR_description)
    default_options = {'max_iter': 100}
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    solver = SolverFactory('ipopt', executable=solver_executable)
    plant_simulator.set_solver(solver, tee=False, default_options=default_options)
    model_simulator.set_solver(solver, tee=False, default_options=default_options)


    # --------read data--------------------
    PE_model_data = pandas.read_csv(algo_results_folder + "PE_model_data.txt", sep='\t', index_col=0, header=0)
    GPE_model_data = pandas.read_csv(algo_results_folder + "GPE_model_data.txt", sep='\t', index_col=0, header=0)
    ISOPE_model_data = pandas.read_csv(algo_results_folder + "ISOPE_model_data.txt", sep='\t', index_col=0, header=0)

    #------------------------------------
    Fb_points = numpy.linspace(start=3, stop=6, num=resolution, endpoint=True)
    Tr_points = numpy.linspace(start=70, stop=100, num=resolution, endpoint=True)

    def plant_simulation_callback(Tr,Fb):
        outputs, solve_status = plant_simulator.simulate({'Tr':Tr,'Fb':Fb}, param_values=None,
                                                             use_homo=True)
        return outputs['profit']

    def model_simulation_callback(Tr,Fb,param_values):
        outputs, solve_status = model_simulator.simulate({'Tr':Tr,'Fb':Fb}, param_values=param_values,
                                                             use_homo=True)
        return outputs['profit']

    #---------------------------
    func = lambda Tr, Fb: plant_simulation_callback(Tr, Fb)
    generate_contour_data(Tr_points, Fb_points, func, data_files['Plant'])

    param_values = {}
    for k in model_simulator.pyomo_model.parameters.keys():
        param_values[k] = model_simulator.pyomo_model.default_value[k]
    func = lambda Tr, Fb: model_simulation_callback(Tr, Fb, param_values)
    generate_contour_data(Tr_points, Fb_points, func, data_files['Model'])

    param_values = {}
    for k in model_simulator.pyomo_model.parameters.keys():
        param_values[k] = PE_model_data.loc[20, k]
    func = lambda Tr,Fb:model_simulation_callback(Tr,Fb,param_values)
    generate_contour_data(Tr_points, Fb_points, func, data_files['PE'])

    param_values = {}
    for k in model_simulator.pyomo_model.parameters.keys():
        param_values[k] = ISOPE_model_data.loc[20, k]
    func = lambda Tr, Fb: model_simulation_callback(Tr, Fb, param_values)
    generate_contour_data(Tr_points, Fb_points, func, data_files['ISOPE'])

    param_values = {}
    for k in model_simulator.pyomo_model.parameters.keys():
        param_values[k] = GPE_model_data.loc[20, k]
    func = lambda Tr, Fb: model_simulation_callback(Tr, Fb, param_values)
    generate_contour_data(Tr_points, Fb_points, func, data_files['GPE'])


if __name__ == "__main__":
    plot_model_maintenance(algo_results_folder, datafile_folder, pic_filename)