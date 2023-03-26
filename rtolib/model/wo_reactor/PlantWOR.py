from rtolib.core.pyomo_model import PyomoModel
from pyomo.environ import *
import os

class RTO_Plant_WO_reactor(PyomoModel):
    def __init__(self):
        self.output_variables = {
            'cost': (lambda m: m.cost_noise_free),
            'XFr_A': (lambda m: m.XFr['A']),
            'XFr_B': (lambda m: m.XFr['B']),
            'XFr_E': (lambda m: m.XFr['E']),
            'XFr_P': (lambda m: m.XFr['P']),
            'XFr_G': (lambda m: m.XFr['G']),
        }
        self.noised_outputs = {
            'cost': (lambda m: m.cost),
            'XFr_A': (lambda m: m.XFr['A'] + m.XFr_noise['A']),
            'XFr_B': (lambda m: m.XFr['B'] + m.XFr_noise['B']),
            'XFr_E': (lambda m: m.XFr['E'] + m.XFr_noise['E']),
            'XFr_P': (lambda m: m.XFr['P'] + m.XFr_noise['P']),
            'XFr_G': (lambda m: m.XFr['G'] + m.XFr_noise['G']),
        }
        self.input_variables = {
            'Fb': (lambda m: m.Fb),
            'Tr': (lambda m: m.Tr),
        }
        self.parameters = {
        }
        self.default_value={}
        self.parameter_scaling_factors={}
        self.output_noise = {
            'XFr_A': (lambda m: m.XFr_noise['A']),
            'XFr_B': (lambda m: m.XFr_noise['B']),
            'XFr_E': (lambda m: m.XFr_noise['E']),
            'XFr_P': (lambda m: m.XFr_noise['P']),
            'XFr_G': (lambda m: m.XFr_noise['G']),
        }
        self.initial_value_file = os.path.join(os.path.dirname(__file__) + r"\wo_reactor_plant_init.txt")

    def build_body(self, wocstr):
        wocstr.Reactions = Set(initialize=[1, 2, 3])
        wocstr.Components = Set(initialize=['A', 'B', 'C', 'E', 'G', 'P'])
        wocstr.Fa = Param(initialize=1.8275)  # kg/s
        wocstr.XFa = Param(wocstr.Components, initialize={'A': 1, 'B': 0, 'C': 0, 'E': 0, 'G': 0, 'P': 0})
        wocstr.XFb = Param(wocstr.Components, initialize={'A': 0, 'B': 1, 'C': 0, 'E': 0, 'G': 0, 'P': 0})
        wocstr.MassHoldup = Param(initialize=2105)  # kg
        wocstr.Ka = Param(wocstr.Reactions, initialize={1: 1.6599e6, 2: 7.2117e8, 3: 2.6745e12})
        wocstr.Kb = Param(wocstr.Reactions, initialize={1: 6666.7, 2: 8333.3, 3: 11111})
        wocstr.MolarMass = Param(wocstr.Components,
                                 initialize={'A': 100, 'B': 100, 'C': 200, 'E': 200, 'G': 300, 'P': 100})  # kg/mol
        stoichiometry = {(1, 'A'): -1, (1, 'B'): -1, (1, 'C'): 1, (1, 'E'): 0, (1, 'G'): 0, (1, 'P'): 0,
                         (2, 'A'): 0, (2, 'B'): -1, (2, 'C'): -1, (2, 'E'): 1, (2, 'G'): 0, (2, 'P'): 1,
                         (3, 'A'): 0, (3, 'B'): 0, (3, 'C'): -1, (3, 'E'): 0, (3, 'G'): 1, (3, 'P'): -1
                         }
        wocstr.stoichiometry = Param(wocstr.Reactions, wocstr.Components, initialize=stoichiometry)
        wocstr.Fb = Var(initialize=4, within=NonNegativeReals, bounds=(3, 6))  # kg/s
        wocstr.Fb.fixed = True
        wocstr.Tr = Var(initialize=90, bounds=(70, 100))  # Celsius degree
        wocstr.Tr.fixed = True
        wocstr.Fr = Var(initialize=5, within=NonNegativeReals)  # kg/s
        wocstr.XFr = Var(wocstr.Components, bounds=(0, 1))
        wocstr.XFr_noise = Var(wocstr.Components, initialize=0)
        for c in wocstr.Components:
            wocstr.XFr_noise[c].fixed=True
        wocstr.K = Var(wocstr.Reactions, initialize=0.1)

        def Kmxx(m, r):
            if r == 1:
                return m.K[1] * m.XFr['A'] * m.XFr['B'] / m.MolarMass['A']
            elif r == 2:
                return m.K[2] * m.XFr['B'] * m.XFr['C'] / m.MolarMass['B']
            elif r == 3:
                return m.K[3] * m.XFr['C'] * m.XFr['P'] / m.MolarMass['C']
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
        def cost_noise_free(m):
            # return -(1143.38 * m.Fr * m.XFr['P'] + 25.92 * m.Fr * m.XFr[
            #     'E'] - 76.23 * m.Fa - 114.34 * m.Fb)
            return -(1143.38 * m.Fr * (cv_func['XFr_P'].__call__(m)) + \
                     25.92 * m.Fr * (cv_func['XFr_E'].__call__(m)) - \
                     76.23 * m.Fa - 114.34 * m.Fb)
        wocstr.cost_noise_free = Expression(rule=cost_noise_free)

        def cost(m):
            # return -(1143.38 * m.Fr * m.XFr['P'] + 25.92 * m.Fr * m.XFr[
            #     'E'] - 76.23 * m.Fa - 114.34 * m.Fb)
            return -(1143.38 * m.Fr * (cv_func['XFr_P'].__call__(m)+m.XFr_noise['P']) + \
                     25.92 * m.Fr * (cv_func['XFr_E'].__call__(m)+m.XFr_noise['E']) - \
                     76.23 * m.Fa - 114.34 * m.Fb)
        wocstr.cost = Expression(rule=cost)