from rtolib.core.pyomo_model import PyomoModel
from pyomo.environ import *
import os

class QuadraticModel_WO_reactor(PyomoModel):
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
            return m.obj+3.032142e+03\
        +-3.161973e+02*m.Fa\
        +8.356869e+01*m.Fb\
        +-7.445638e+01*m.Tr\
        +2.566075e+01*m.Fa*m.Fa\
        +-1.190035e+01*m.Fa*m.Fb\
        +1.946919e+00*m.Fa*m.Tr\
        +-1.190035e+01*m.Fb*m.Fa\
        +5.815342e+00*m.Fb*m.Fb\
        +-6.087520e-01*m.Fb*m.Tr\
        +1.946919e+00*m.Tr*m.Fa\
        +-6.087520e-01*m.Tr*m.Fb\
        +4.395533e-01*m.Tr*m.Tr
        wocstr.cost = Expression(rule=cost)

        def purity_constraint(m):
            return 7.374878e-01\
            +-2.455611e-01*m.Fa\
            +-3.100480e-02*m.Fb\
            +-7.467196e-03*m.Tr\
            +1.623000e-02*m.Fa*m.Fa\
            +2.480823e-03*m.Fa*m.Fb\
            +5.541484e-04*m.Fa*m.Tr\
            +2.480823e-03*m.Fb*m.Fa\
            +1.869133e-03*m.Fb*m.Fb\
            +-1.768741e-04*m.Fb*m.Tr\
            +5.541484e-04*m.Tr*m.Fa\
            +-1.768741e-04*m.Tr*m.Fb\
            +6.484424e-05*m.Tr*m.Tr
        wocstr.purity_constraint = Expression(rule=purity_constraint)