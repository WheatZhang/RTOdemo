from rtolib.core.pyomo_model import PyomoModel
from pyomo.environ import *
import os

class quadr_con_wrong_curv_model(PyomoModel):
    def __init__(self):
        self.output_variables = {
            'cost': (lambda m: m.cost_var),
            'con': (lambda m: m.con_var),
        }
        self.input_variables = {
            'u1': (lambda m: m.u1),
            'u2': (lambda m: m.u2),
        }
        self.parameters = {
        }
        self.default_value = {}
        self.parameter_scaling_factors = {}
        self.output_noise = {
        }
        self.initial_value_file = os.path.join(os.path.dirname(__file__) + r"\quadratic_init.txt")

    def build_body(self, model):
        model.u1 = Var(initialize=1, bounds=(-2, 4))
        model.u2 = Var(initialize=1, bounds=(-4, 2))
        model.cost_var = Var(initialize=0)
        model.con_var = Var(initialize=0)

        def cost(m):
            # return 2*m.u1 * m.u1 + 2*m.u2 * m.u2 == m.cost_var
            return m.u1 * m.u1 + m.u2 * m.u2 == m.cost_var

        model.cost = Constraint(rule=cost)

        def constr(m):
            return 1 - m.u1 - m.u2 * m.u2 == m.con_var
            # return 1 - m.u1 + 2*m.u2 * m.u2 +m.u1*m.u1+ 2*m.u2 == m.con_var

        model.constr = Constraint(rule=constr)


    def build_rto(self, model, cv_func):
        pass
