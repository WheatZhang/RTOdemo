from rtolib.core.pyomo_model import PyomoModel
from pyomo.environ import *
import os

class QuadraticModel_Parallel_WO_reactor(PyomoModel):
    def __init__(self):
        self.output_variables = {
            'cost': (lambda m: m.cost),
            'con': (lambda m: m.Fb_con),
        }
        self.noised_outputs = {
        }
        self.input_variables = {
            "Fb1": (lambda m: m.reactor_block[1].Fb),
            "Fb2": (lambda m: m.reactor_block[2].Fb),
            "Fb3": (lambda m: m.reactor_block[3].Fb),
            "Tr1": (lambda m: m.reactor_block[1].Tr),
            "Tr2": (lambda m: m.reactor_block[2].Tr),
            "Tr3": (lambda m: m.reactor_block[3].Tr),
        }
        self.parameters = {
        }
        self.default_value = {
        }
        self.parameter_scaling_factors = {
        }
        self.initial_value_file = os.path.join(os.path.dirname(__file__) + r"\plant_init_value.txt")

    def build_body(self, parallel_wocstr):
        parallel_wocstr.block_index = Set(initialize=[1, 2, 3])
        def block_rule(wocstr):
            wocstr.Fb = Var(initialize=4, within=NonNegativeReals, bounds=(0, 10))  # kg/s
            wocstr.Fb.fixed = True
            wocstr.Tr = Var(initialize=90, bounds=(0, 200))  # Celsius degree
            wocstr.Tr.fixed = True

        parallel_wocstr.reactor_block = Block(parallel_wocstr.block_index, rule=block_rule)

        parallel_wocstr.obj = Var(initialize=0)
        parallel_wocstr.e = Var(initialize=0)

        def eq1(m):
            return m.e == 0.0*(\
                sum([m.reactor_block[i].Fb for i in m.block_index])+\
                sum([m.reactor_block[i].Tr for i in m.block_index])\
                    )
        parallel_wocstr.eq1 = Constraint(rule=eq1)

        def _obj(m):
            return 0  #10*(m.Fb**2+m.Tr**2)
        parallel_wocstr._obj = Expression(rule=_obj)

    def build_rto(self, wocstr, cv_func):
        def cost(m):
            return m._obj+2.882206e+03\
            +-1.028687e+01*m.reactor_block[1].Fb\
            +-4.871776e+01*m.reactor_block[2].Fb\
            +-7.968971e+01*m.reactor_block[3].Fb\
            +-2.712656e+01*m.reactor_block[1].Tr\
            +-2.302494e+01*m.reactor_block[2].Tr\
            +-2.932710e+01*m.reactor_block[3].Tr\
            +1.207003e+01*m.reactor_block[1].Fb*m.reactor_block[1].Fb\
            +6.674691e-03*m.reactor_block[1].Fb*m.reactor_block[2].Fb\
            +1.778422e-01*m.reactor_block[1].Fb*m.reactor_block[3].Fb\
            +-4.808340e-01*m.reactor_block[1].Fb*m.reactor_block[1].Tr\
            +-8.529660e-03*m.reactor_block[1].Fb*m.reactor_block[2].Tr\
            +-2.959244e-02*m.reactor_block[1].Fb*m.reactor_block[3].Tr\
            +6.674691e-03*m.reactor_block[2].Fb*m.reactor_block[1].Fb\
            +1.270260e+01*m.reactor_block[2].Fb*m.reactor_block[2].Fb\
            +6.543841e-01*m.reactor_block[2].Fb*m.reactor_block[3].Fb\
            +-2.597251e-03*m.reactor_block[2].Fb*m.reactor_block[1].Tr\
            +-5.218354e-01*m.reactor_block[2].Fb*m.reactor_block[2].Tr\
            +7.992147e-03*m.reactor_block[2].Fb*m.reactor_block[3].Tr\
            +1.778422e-01*m.reactor_block[3].Fb*m.reactor_block[1].Fb\
            +6.543841e-01*m.reactor_block[3].Fb*m.reactor_block[2].Fb\
            +1.014296e+01*m.reactor_block[3].Fb*m.reactor_block[3].Fb\
            +-9.907168e-04*m.reactor_block[3].Fb*m.reactor_block[1].Tr\
            +1.234329e-02*m.reactor_block[3].Fb*m.reactor_block[2].Tr\
            +-2.970644e-01*m.reactor_block[3].Fb*m.reactor_block[3].Tr\
            +-4.808340e-01*m.reactor_block[1].Tr*m.reactor_block[1].Fb\
            +-2.597251e-03*m.reactor_block[1].Tr*m.reactor_block[2].Fb\
            +-9.907168e-04*m.reactor_block[1].Tr*m.reactor_block[3].Fb\
            +2.017913e-01*m.reactor_block[1].Tr*m.reactor_block[1].Tr\
            +-6.195867e-03*m.reactor_block[1].Tr*m.reactor_block[2].Tr\
            +-1.092047e-04*m.reactor_block[1].Tr*m.reactor_block[3].Tr\
            +-8.529660e-03*m.reactor_block[2].Tr*m.reactor_block[1].Fb\
            +-5.218354e-01*m.reactor_block[2].Tr*m.reactor_block[2].Fb\
            +1.234329e-02*m.reactor_block[2].Tr*m.reactor_block[3].Fb\
            +-6.195867e-03*m.reactor_block[2].Tr*m.reactor_block[1].Tr\
            +1.818400e-01*m.reactor_block[2].Tr*m.reactor_block[2].Tr\
            +6.653902e-03*m.reactor_block[2].Tr*m.reactor_block[3].Tr\
            +-2.959244e-02*m.reactor_block[3].Tr*m.reactor_block[1].Fb\
            +7.992147e-03*m.reactor_block[3].Tr*m.reactor_block[2].Fb\
            +-2.970644e-01*m.reactor_block[3].Tr*m.reactor_block[3].Fb\
            +-1.092047e-04*m.reactor_block[3].Tr*m.reactor_block[1].Tr\
            +6.653902e-03*m.reactor_block[3].Tr*m.reactor_block[2].Tr\
            +2.072300e-01*m.reactor_block[3].Tr*m.reactor_block[3].Tr
        wocstr.cost = Expression(rule=cost)

        def Fb_con(m):
            return -1.200000e+01\
            +1.000000e+00*m.reactor_block[1].Fb\
            +1.000000e+00*m.reactor_block[2].Fb\
            +1.000000e+00*m.reactor_block[3].Fb\
            +-2.316181e-11*m.reactor_block[1].Tr\
            +3.927039e-11*m.reactor_block[2].Tr\
            +5.943451e-11*m.reactor_block[3].Tr\
            +8.310227e-09*m.reactor_block[1].Fb*m.reactor_block[1].Fb\
            +-3.822154e-10*m.reactor_block[1].Fb*m.reactor_block[2].Fb\
            +1.308946e-10*m.reactor_block[1].Fb*m.reactor_block[3].Fb\
            +2.427230e-12*m.reactor_block[1].Fb*m.reactor_block[1].Tr\
            +-8.805013e-12*m.reactor_block[1].Fb*m.reactor_block[2].Tr\
            +-1.041662e-11*m.reactor_block[1].Fb*m.reactor_block[3].Tr\
            +-3.822154e-10*m.reactor_block[2].Fb*m.reactor_block[1].Fb\
            +2.160638e-09*m.reactor_block[2].Fb*m.reactor_block[2].Fb\
            +1.237324e-09*m.reactor_block[2].Fb*m.reactor_block[3].Fb\
            +5.460131e-12*m.reactor_block[2].Fb*m.reactor_block[1].Tr\
            +-4.415342e-12*m.reactor_block[2].Fb*m.reactor_block[2].Tr\
            +-1.188666e-11*m.reactor_block[2].Fb*m.reactor_block[3].Tr\
            +1.308946e-10*m.reactor_block[3].Fb*m.reactor_block[1].Fb\
            +1.237324e-09*m.reactor_block[3].Fb*m.reactor_block[2].Fb\
            +4.248171e-09*m.reactor_block[3].Fb*m.reactor_block[3].Fb\
            +-2.867380e-12*m.reactor_block[3].Fb*m.reactor_block[1].Tr\
            +-1.395684e-11*m.reactor_block[3].Fb*m.reactor_block[2].Tr\
            +-1.427420e-11*m.reactor_block[3].Fb*m.reactor_block[3].Tr\
            +2.427230e-12*m.reactor_block[1].Tr*m.reactor_block[1].Fb\
            +5.460131e-12*m.reactor_block[1].Tr*m.reactor_block[2].Fb\
            +-2.867380e-12*m.reactor_block[1].Tr*m.reactor_block[3].Fb\
            +8.744075e-13*m.reactor_block[1].Tr*m.reactor_block[1].Tr\
            +-2.601797e-13*m.reactor_block[1].Tr*m.reactor_block[2].Tr\
            +4.941721e-14*m.reactor_block[1].Tr*m.reactor_block[3].Tr\
            +-8.805013e-12*m.reactor_block[2].Tr*m.reactor_block[1].Fb\
            +-4.415342e-12*m.reactor_block[2].Tr*m.reactor_block[2].Fb\
            +-1.395684e-11*m.reactor_block[2].Tr*m.reactor_block[3].Fb\
            +-2.601797e-13*m.reactor_block[2].Tr*m.reactor_block[1].Tr\
            +4.597747e-13*m.reactor_block[2].Tr*m.reactor_block[2].Tr\
            +-5.091232e-13*m.reactor_block[2].Tr*m.reactor_block[3].Tr\
            +-1.041662e-11*m.reactor_block[3].Tr*m.reactor_block[1].Fb\
            +-1.188666e-11*m.reactor_block[3].Tr*m.reactor_block[2].Fb\
            +-1.427420e-11*m.reactor_block[3].Tr*m.reactor_block[3].Fb\
            +4.941721e-14*m.reactor_block[3].Tr*m.reactor_block[1].Tr\
            +-5.091232e-13*m.reactor_block[3].Tr*m.reactor_block[2].Tr\
            +1.645048e-12*m.reactor_block[3].Tr*m.reactor_block[3].Tr
        wocstr.Fb_con = Expression(rule=Fb_con)