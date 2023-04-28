from rtolib.core.pyomo_model import PyomoModel
from pyomo.environ import *
import os

class constrained_quadratic_problem_plant(PyomoModel):
    def __init__(self):
        self.output_variables = {
            'cost': (lambda m: m.cost),
            'con': (lambda m: m.con),
        }
        self.noised_outputs = {
            'cost': (lambda m: m.cost),
            'con': (lambda m: m.con),
        }
        self.input_variables = {
            'u1': (lambda m: m.u1),
        }
        self.parameters = {
        }
        self.default_value={}
        self.parameter_scaling_factors={}
        self.output_noise = {
        }
        self.initial_value_file = os.path.join(os.path.dirname(__file__) + r"\quadratic_init.txt")

    def build_body(self, model):
        model.u1 = Var(initialize=1)
        model.cost_var = Var(initialize=10)
        model.cost_var2 = Var(initialize=10)
        def cost_var_con(m):
            return m.cost_var == m.cost_var2
        model.cost_var_con = Constraint(rule=cost_var_con)

    def build_rto(self, model, cv_func):
        def cost(m):
            # return 1/50*(m.u1+1)**4+m.cost_var**2-m.cost_var2**2
            return 1/2*(m.u1+1)**2+m.cost_var**2-m.cost_var2**2
        model.cost = Expression(rule=cost)
        def con(m):
            return 1/50*(m.u1-2)**4-10
        model.con = Expression(rule=con)
