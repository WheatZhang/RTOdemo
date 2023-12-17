from rtolib.core.pyomo_model import PyomoModel
from pyomo.environ import *
import os

class ZeroModel_WO_reactor(PyomoModel):
    def __init__(self):
        self.output_variables = {
            'cost': (lambda m: m.cost),
            'con': (lambda m: m.purity_constraint),
            'XFr_A': (lambda m: m.XFr['A']),
            'XFr_B': (lambda m: m.XFr['B']),
            'XFr_E': (lambda m: m.XFr['E']),
            'XFr_P': (lambda m: m.XFr['P']),
            'XFr_G': (lambda m: m.XFr['G']),
        }
        self.input_variables = {
            'Fa': (lambda m: m.Fa),
            'Fb': (lambda m: m.Fb),
            'Tr': (lambda m: m.Tr),
        }
        self.parameters = {
        }
        self.default_value = {
        }
        self.parameter_scaling_factors = {
        }
        self.initial_value_file = os.path.join(os.path.dirname(__file__) + r"\wo_reactor_model_init.txt")

    def build_body(self, wocstr):
        wocstr.Fa = Var(initialize=4, within=NonNegativeReals, bounds=(3, 4.5))  # kg/s
        wocstr.Fa.fixed = True
        wocstr.Fb = Var(initialize=9, within=NonNegativeReals, bounds=(6, 11))  # kg/s
        wocstr.Fb.fixed = True
        wocstr.Tr = Var(initialize=90, bounds=(80, 105))  # Celsius degree
        wocstr.Tr.fixed = True
        wocstr.obj = Var(initialize=0)
        wocstr.e = Var(initialize=0)
        wocstr.Components = Set(initialize=['A', 'B', 'C', 'E', 'G', 'P'])
        wocstr.XFr = Var(wocstr.Components, initialize=0)

        def eq1(m):
            return m.e == 0.0*(m.Fb + m.Tr)
        wocstr.eq1 = Constraint(rule=eq1)

        def obj(m):
            return 0  #10*(m.Fb**2+m.Tr**2)
        wocstr.obj = Expression(rule=obj)

    def build_rto(self, wocstr, cv_func):
        def cost(m):
            return m.obj
        wocstr.cost = Expression(rule=cost)

        def purity_constraint(m):
            return 0
        wocstr.purity_constraint = Expression(rule=purity_constraint)