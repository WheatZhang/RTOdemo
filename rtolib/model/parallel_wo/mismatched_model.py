from rtolib.core.pyomo_model import PyomoModel
from pyomo.environ import *
import os

class RTO_Mismatched_Parallel_WO(PyomoModel):
    def __init__(self):
        self.output_variables = {
            'cost': (lambda m: m.cost),
            'con': (lambda m: m.Fb_con),
        }
        self.noised_outputs = {
        }
        self.input_variables = {
            "Fb1": (lambda m: m.Fb1),
            "Fb2": (lambda m: m.Fb2),
            "Fb3": (lambda m: m.Fb3),
            "Tr1": (lambda m: m.Tr1),
            "Tr2": (lambda m: m.Tr2),
            "Tr3": (lambda m: m.Tr3),
        }
        self.parameters = {
        }
        self.default_value={}
        self.parameter_scaling_factors={}
        self.output_noise = {
        }
        self.initial_value_file = os.path.join(os.path.dirname(__file__) + r"\plant_init_value.txt")

    def build_body(self, wocstr):
        parallel_wocstr.block_index = Set(initialize=[1, 2, 3])

        def block_rule(wocstr):
            wocstr.Reactions = Set(initialize=[1, 2])
            wocstr.Components = Set(initialize=['A', 'B', 'C', 'E', 'G', 'P'])
            wocstr.Fa = Var(initialize=1.8275)  # kg/s
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
            wocstr.Fb = Var(initialize=4, within=NonNegativeReals, bounds=(0, 10))  # kg/s
            wocstr.Tr = Var(initialize=90, bounds=(0, 200))  # Celsius degree
            wocstr.Fr = Var(initialize=5, within=NonNegativeReals)  # kg/s
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

            def profit(m):
                return 1143.38 * m.Fr * m.XFr['P'] + 25.92 * m.Fr * m.XFr['E'] - 76.23 * m.Fa - 114.34 * m.Fb

            wocstr.profit = Expression(rule=profit)

        parallel_wocstr.reactor_block = Block(parallel_wocstr.block_index, rule=block_rule)
        Fa = {1: 1.5, 2: 1.8, 3: 2.1}
        for i in parallel_wocstr.block_index:
            parallel_wocstr.reactor_block[i].Fa = Fa[i]

    def build_rto(self, wocstr, cv_func):
        def Fb_con(m):
            return sum([parallel_wocstr.reactor_block[i].Fb for i in parallel_wocstr.block_index]) - 12

        parallel_wocstr.Fb_con = Expression(rule=Fb_con)

        def cost(m):
            return -sum([parallel_wocstr.reactor_block[i].profit for i in parallel_wocstr.block_index])

        parallel_wocstr.cost = Expression(rule=cost)