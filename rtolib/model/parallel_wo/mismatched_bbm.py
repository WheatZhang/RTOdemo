from rtolib.core.black_box_model import BlackBoxModel
from rtolib.core.solve import PyomoSimulator
from pyomo.environ import SolverFactory
from .mismatched_model import RTO_Mismatched_Parallel_WO
from .description import default_WOR_description

class RTO_Mismatched_BlackBox_Parallel_WO(BlackBoxModel):
    def build(self):
        self.model_simulator = PyomoSimulator(RTO_Mismatched_Parallel_WO())
        self.model_simulator.build(default_WOR_description)
        default_options = {'max_iter': 100}
        solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
        self.simulation_solver = SolverFactory('ipopt', executable=solver_executable)
        self.model_simulator.set_solver(self.simulation_solver, tee=False,\
                                        default_options=default_options)

    def simulate(self, input_dict):
        output_dict, solve_status = self.model_simulator.simulate(input_dict, param_values={},
                                                             use_homo=False)
        return output_dict