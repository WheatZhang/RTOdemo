from rtolib.core.black_box_model import BlackBoxModel
from rtolib.core.solve import PyomoSimulator
from pyomo.environ import SolverFactory
from .quadraticmodelParallelWO import QuadraticModel_Parallel_WO_reactor
from .description import default_WOR_description

class QuadraticBlackBoxModel_Parallel_WO(BlackBoxModel):
    def build(self):
        pass
        # self.model_simulator = PyomoSimulator(QuadraticModel_Parallel_WO_reactor())
        # self.model_simulator.build(default_WOR_description)
        # default_options = {'max_iter': 100}
        # solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
        # self.simulation_solver = SolverFactory('ipopt', executable=solver_executable)
        # self.model_simulator.set_solver(self.simulation_solver, tee=False,\
        #                                 default_options=default_options)

    def simulate(self, input_dict):
        output_dict = {}
        output_dict['cost'] = 2.882206e+03\
            +-1.028687e+01*input_dict["Fb1"]\
            +-4.871776e+01*input_dict["Fb2"]\
            +-7.968971e+01*input_dict["Fb3"]\
            +-2.712656e+01*input_dict["Tr1"]\
            +-2.302494e+01*input_dict["Tr2"]\
            +-2.932710e+01*input_dict["Tr3"]\
            +1.207003e+01*input_dict["Fb1"]*input_dict["Fb1"]\
            +6.674691e-03*input_dict["Fb1"]*input_dict["Fb2"]\
            +1.778422e-01*input_dict["Fb1"]*input_dict["Fb3"]\
            +-4.808340e-01*input_dict["Fb1"]*input_dict["Tr1"]\
            +-8.529660e-03*input_dict["Fb1"]*input_dict["Tr2"]\
            +-2.959244e-02*input_dict["Fb1"]*input_dict["Tr3"]\
            +6.674691e-03*input_dict["Fb2"]*input_dict["Fb1"]\
            +1.270260e+01*input_dict["Fb2"]*input_dict["Fb2"]\
            +6.543841e-01*input_dict["Fb2"]*input_dict["Fb3"]\
            +-2.597251e-03*input_dict["Fb2"]*input_dict["Tr1"]\
            +-5.218354e-01*input_dict["Fb2"]*input_dict["Tr2"]\
            +7.992147e-03*input_dict["Fb2"]*input_dict["Tr3"]\
            +1.778422e-01*input_dict["Fb3"]*input_dict["Fb1"]\
            +6.543841e-01*input_dict["Fb3"]*input_dict["Fb2"]\
            +1.014296e+01*input_dict["Fb3"]*input_dict["Fb3"]\
            +-9.907168e-04*input_dict["Fb3"]*input_dict["Tr1"]\
            +1.234329e-02*input_dict["Fb3"]*input_dict["Tr2"]\
            +-2.970644e-01*input_dict["Fb3"]*input_dict["Tr3"]\
            +-4.808340e-01*input_dict["Tr1"]*input_dict["Fb1"]\
            +-2.597251e-03*input_dict["Tr1"]*input_dict["Fb2"]\
            +-9.907168e-04*input_dict["Tr1"]*input_dict["Fb3"]\
            +2.017913e-01*input_dict["Tr1"]*input_dict["Tr1"]\
            +-6.195867e-03*input_dict["Tr1"]*input_dict["Tr2"]\
            +-1.092047e-04*input_dict["Tr1"]*input_dict["Tr3"]\
            +-8.529660e-03*input_dict["Tr2"]*input_dict["Fb1"]\
            +-5.218354e-01*input_dict["Tr2"]*input_dict["Fb2"]\
            +1.234329e-02*input_dict["Tr2"]*input_dict["Fb3"]\
            +-6.195867e-03*input_dict["Tr2"]*input_dict["Tr1"]\
            +1.818400e-01*input_dict["Tr2"]*input_dict["Tr2"]\
            +6.653902e-03*input_dict["Tr2"]*input_dict["Tr3"]\
            +-2.959244e-02*input_dict["Tr3"]*input_dict["Fb1"]\
            +7.992147e-03*input_dict["Tr3"]*input_dict["Fb2"]\
            +-2.970644e-01*input_dict["Tr3"]*input_dict["Fb3"]\
            +-1.092047e-04*input_dict["Tr3"]*input_dict["Tr1"]\
            +6.653902e-03*input_dict["Tr3"]*input_dict["Tr2"]\
            +2.072300e-01*input_dict["Tr3"]*input_dict["Tr3"]

        output_dict['con'] = -1.200000e+01\
            +1.000000e+00*input_dict["Fb1"]\
            +1.000000e+00*input_dict["Fb2"]\
            +1.000000e+00*input_dict["Fb3"]\
            +-2.316181e-11*input_dict["Tr1"]\
            +3.927039e-11*input_dict["Tr2"]\
            +5.943451e-11*input_dict["Tr3"]\
            +8.310227e-09*input_dict["Fb1"]*input_dict["Fb1"]\
            +-3.822154e-10*input_dict["Fb1"]*input_dict["Fb2"]\
            +1.308946e-10*input_dict["Fb1"]*input_dict["Fb3"]\
            +2.427230e-12*input_dict["Fb1"]*input_dict["Tr1"]\
            +-8.805013e-12*input_dict["Fb1"]*input_dict["Tr2"]\
            +-1.041662e-11*input_dict["Fb1"]*input_dict["Tr3"]\
            +-3.822154e-10*input_dict["Fb2"]*input_dict["Fb1"]\
            +2.160638e-09*input_dict["Fb2"]*input_dict["Fb2"]\
            +1.237324e-09*input_dict["Fb2"]*input_dict["Fb3"]\
            +5.460131e-12*input_dict["Fb2"]*input_dict["Tr1"]\
            +-4.415342e-12*input_dict["Fb2"]*input_dict["Tr2"]\
            +-1.188666e-11*input_dict["Fb2"]*input_dict["Tr3"]\
            +1.308946e-10*input_dict["Fb3"]*input_dict["Fb1"]\
            +1.237324e-09*input_dict["Fb3"]*input_dict["Fb2"]\
            +4.248171e-09*input_dict["Fb3"]*input_dict["Fb3"]\
            +-2.867380e-12*input_dict["Fb3"]*input_dict["Tr1"]\
            +-1.395684e-11*input_dict["Fb3"]*input_dict["Tr2"]\
            +-1.427420e-11*input_dict["Fb3"]*input_dict["Tr3"]\
            +2.427230e-12*input_dict["Tr1"]*input_dict["Fb1"]\
            +5.460131e-12*input_dict["Tr1"]*input_dict["Fb2"]\
            +-2.867380e-12*input_dict["Tr1"]*input_dict["Fb3"]\
            +8.744075e-13*input_dict["Tr1"]*input_dict["Tr1"]\
            +-2.601797e-13*input_dict["Tr1"]*input_dict["Tr2"]\
            +4.941721e-14*input_dict["Tr1"]*input_dict["Tr3"]\
            +-8.805013e-12*input_dict["Tr2"]*input_dict["Fb1"]\
            +-4.415342e-12*input_dict["Tr2"]*input_dict["Fb2"]\
            +-1.395684e-11*input_dict["Tr2"]*input_dict["Fb3"]\
            +-2.601797e-13*input_dict["Tr2"]*input_dict["Tr1"]\
            +4.597747e-13*input_dict["Tr2"]*input_dict["Tr2"]\
            +-5.091232e-13*input_dict["Tr2"]*input_dict["Tr3"]\
            +-1.041662e-11*input_dict["Tr3"]*input_dict["Fb1"]\
            +-1.188666e-11*input_dict["Tr3"]*input_dict["Fb2"]\
            +-1.427420e-11*input_dict["Tr3"]*input_dict["Fb3"]\
            +4.941721e-14*input_dict["Tr3"]*input_dict["Tr1"]\
            +-5.091232e-13*input_dict["Tr3"]*input_dict["Tr2"]\
            +1.645048e-12*input_dict["Tr3"]*input_dict["Tr3"]
        return output_dict