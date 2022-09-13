from rtolib.model.wo_reactor_con import RTO_Plant_WO_reactor,\
    RTO_Mismatched_WO_reactor, RTO_Mismatched_WO_reactor_QC, default_WOR_description
from rtolib.core.solve import PyomoSimulator
from pyomo.environ import SolverFactory
import numpy
import rtolib.util.init_value as ivt


plant_simulator = PyomoSimulator(RTO_Plant_WO_reactor())
plant_simulator.build(default_WOR_description)
default_options = {'max_iter': 100}
solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
solver = SolverFactory('ipopt', executable=solver_executable)
plant_simulator.set_solver(solver, tee=True, default_options=default_options)

plant_simulator.simulate({'Fa': 4, 'Fb': 9, 'Tr':90}, param_values=None,
                                                     use_homo=False)
ivt.to_template(plant_simulator.model, "data/plant_wo_init_value.txt")

#------------------------------------------------------------
model_simulator = PyomoSimulator(RTO_Mismatched_WO_reactor())
model_simulator.build(default_WOR_description)
default_options = {'max_iter': 100}
solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
solver = SolverFactory('ipopt', executable=solver_executable)
model_simulator.set_solver(solver, tee=True, default_options=default_options)

model_simulator.simulate({'Fa': 4, 'Fb': 9, 'Tr':90}, param_values=None,
                                                     use_homo=False)
ivt.to_template(model_simulator.model, "data/model_wo_init_value.txt")

