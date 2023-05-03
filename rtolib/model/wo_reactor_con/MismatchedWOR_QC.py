from rtolib.core.pyomo_model import PyomoModel
from pyomo.environ import *
import os

class RTO_Mismatched_WO_reactor_QC(PyomoModel):
    '''
    QC = With Quadratic Correction terms
    '''
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
            'Ka1': (lambda m: m.Ka[1]),
            'Ka2': (lambda m: m.Ka[2]),
            'Kb1': (lambda m: m.Kb[1]),
            'Kb2': (lambda m: m.Kb[2]),
        }
        self.default_value = {
            'Ka1': 2.189e8,
            'Ka2': 4.310e13,
            'Kb1': 8077.6,
            'Kb2': 12438,
        }
        self.parameter_scaling_factors = {
            'Ka1': 2.189e8,
            'Ka2': 4.310e13,
            'Kb1': 8077.6,
            'Kb2': 12438,
        }
        self.initial_value_file = os.path.join(os.path.dirname(__file__) + r"\wo_reactor_model_init.txt")

    def build_body(self, wocstr):
        wocstr.Reactions = Set(initialize=[1, 2])
        wocstr.Components = Set(initialize=['A', 'B', 'C', 'E', 'G', 'P'])
        wocstr.Fa = Var(initialize=4, within=NonNegativeReals, bounds=(3, 4.5))  # kg/s
        wocstr.Fa.fixed = True
        wocstr.XFa = Param(wocstr.Components, initialize={'A': 1, 'B': 0, 'C': 0, 'E': 0, 'G': 0, 'P': 0})
        wocstr.XFb = Param(wocstr.Components, initialize={'A': 0, 'B': 1, 'C': 0, 'E': 0, 'G': 0, 'P': 0})
        wocstr.MassHoldup = Param(initialize=2105)  # kg
        wocstr.Ka = Var(wocstr.Reactions, initialize={1: 2.189e8, 2: 4.310e13})
        wocstr.Ka[1].lb = 2.189e8 * 0.2
        wocstr.Ka[1].ub = 2.189e8 * 5
        wocstr.Ka[2].lb = 4.310e13 * 0.2
        wocstr.Ka[2].ub = 4.310e13 * 5
        wocstr.Kb = Var(wocstr.Reactions, initialize={1: 8077.6, 2: 12438})
        wocstr.Kb[1].lb = 8077.6 * 0.5
        wocstr.Kb[1].ub = 8077.6 * 1.5
        wocstr.Kb[2].lb = 12438 * 0.5
        wocstr.Kb[2].ub = 12438 * 1.5
        for temp in wocstr.Reactions:
            wocstr.Ka[temp].fixed = True
            wocstr.Kb[temp].fixed = True
        wocstr.MolarMass = Param(wocstr.Components,
                                 initialize={'A': 100, 'B': 100, 'C': 200, 'E': 200, 'G': 300, 'P': 100})  # kg/mol
        stoichiometry = {(1, 'A'): -1, (1, 'B'): -2, (1, 'C'): 0, (1, 'E'): 1, (1, 'G'): 0, (1, 'P'): 1,
                         (2, 'A'): -1, (2, 'B'): -1, (2, 'C'): 0, (2, 'E'): 0, (2, 'G'): 1, (2, 'P'): -1
                         }
        wocstr.stoichiometry = Param(wocstr.Reactions, wocstr.Components, initialize=stoichiometry)
        wocstr.Fb = Var(initialize=9, within=NonNegativeReals, bounds=(6, 11))  # kg/s
        wocstr.Fb.fixed = True
        wocstr.Tr = Var(initialize=90, bounds=(80, 105))  # Celsius degree
        wocstr.Tr.fixed = True
        wocstr.Fr = Var(initialize=13, within=NonNegativeReals)  # kg/s
        wocstr.XFr = Var(wocstr.Components, bounds=(0, 1))
        wocstr.K = Var(wocstr.Reactions, initialize=0.1)

        def Kmxx(m, r):
            if r == 1:
                return m.K[1] * m.XFr['A'] * m.XFr['B'] * m.XFr['B'] / m.MolarMass['A']
            elif r == 2:
                return m.K[2] * m.XFr['A'] * m.XFr['B'] * m.XFr['P'] / m.MolarMass['A']
        wocstr.Kmxx = Expression(wocstr.Reactions, rule=Kmxx)

        def react_massblnc(m, comp):
            return m.Fa * m.XFa[comp] + m.Fb * m.XFb[comp] - m.Fr * m.XFr[comp] + sum(
                [m.MolarMass[comp] * m.stoichiometry[r, comp] * m.Kmxx[r] * m.MassHoldup for r in m.Reactions]) == 0
        wocstr.react_massblnc = Constraint(wocstr.Components, rule=react_massblnc)

        def overall_massblnc(m):
            return m.Fa + m.Fb == m.Fr
        wocstr.overall_massblnc = Constraint(rule=overall_massblnc)

        def kinetic_coff(m, r):
            return (m.K[r] - m.Ka[r] * exp(-m.Kb[r] / (m.Tr + 273.15))) * 1e1 == 0
        wocstr.kinetic_coff = Constraint(wocstr.Reactions, rule=kinetic_coff)

    def build_rto(self, wocstr, cv_func):
        def cost(m):
            return -(1143.38 * m.Fr * cv_func['XFr_P'].__call__(m) +\
                     25.92 * m.Fr * cv_func['XFr_E'].__call__(m) -\
                     76.23 * m.Fa - 114.34 * m.Fb)+\
                   74.38/2*(m.Fa**2+m.Fb**2+m.Tr**2)
            # return -(1143.38 * m.Fr * cv_func['XFr_P'].__call__(m) + \
            #          25.92 * m.Fr * cv_func['XFr_E'].__call__(m) - \
            #          76.23 * m.Fa - 114.34 * m.Fb) + \
            #        0 / 2 * (m.Fa ** 2 + m.Fb ** 2 + m.Tr ** 2) #0 #74.38
        wocstr.cost = Expression(rule=cost)

        def purity_constraint(m):
            return cv_func['XFr_G'].__call__(m) - 0.08+\
                   0.0563/2*(m.Fa**2+m.Fb**2+m.Tr**2)

            # return cv_func['XFr_G'].__call__(m) - 0.08 + \
            #        0.01 / 2 * (m.Fa ** 2 + m.Fb ** 2 + m.Tr ** 2) # 0.01 # 0.0563
        wocstr.purity_constraint = Expression(rule=purity_constraint)