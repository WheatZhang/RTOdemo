#!/usr/bin/env python
#-*- coding:utf-8 -*-
import numpy as np
from pyomo.environ import *
import rtolib.util.init_value as init_value

def plant_simulator():
    wocstr = ConcreteModel()
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
    wocstr.Fb = Var(initialize=4, within=NonNegativeReals,bounds=(3,6))  # kg/s
    wocstr.Fb.fixed = True
    wocstr.Tr = Var(initialize=90,bounds=(70,100))  # Celsius degree
    wocstr.Tr.fixed = True
    wocstr.Fr = Var(initialize=5, within=NonNegativeReals)  # kg/s
    wocstr.XFr = Var(wocstr.Components, bounds=(0, 1))
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
    def profit(m):
        return 1143.38 * m.Fr * m.XFr['P'] + 25.92 * m.Fr * m.XFr['E'] - 76.23 * m.Fa - 114.34 * m.Fb
    wocstr.profit = Expression(rule=profit)
    def obj(m):
        return m.profit
    wocstr.objective = Objective(rule=obj, sense=maximize)
    return wocstr

def plant_simulate(simulator, solver, Fb, Tr):
    simulator.Fb.fixed = True
    simulator.Tr.fixed = True
    simulator.Fb = Fb
    simulator.Tr = Tr
    solver.solve(simulator, tee=True)
    return value(simulator.profit)

def get_mismatched_model():
    wocstr = ConcreteModel()
    wocstr.Reactions = Set(initialize=[1, 2])
    wocstr.Components = Set(initialize=['A', 'B', 'C', 'E', 'G', 'P'])
    wocstr.Fa = Param(initialize=1.8275)  # kg/s
    wocstr.XFa = Param(wocstr.Components, initialize={'A': 1, 'B': 0, 'C': 0, 'E': 0, 'G': 0, 'P': 0})
    wocstr.XFb = Param(wocstr.Components, initialize={'A': 0, 'B': 1, 'C': 0, 'E': 0, 'G': 0, 'P': 0})
    wocstr.MassHoldup = Param(initialize=2105)  # kg
    wocstr.Ka = Var(wocstr.Reactions, initialize={1: 2.189e8, 2: 4.310e13})
    wocstr.Kb = Var(wocstr.Reactions, initialize={1: 8077.6, 2: 12438})
    for temp in wocstr.Reactions:
        wocstr.Ka[temp].fixed = True
        wocstr.Kb[temp].fixed = True
    wocstr.MolarMass = Param(wocstr.Components,
                             initialize={'A': 100, 'B': 100, 'C': 200, 'E': 200, 'G': 300, 'P': 100})  # kg/mol
    stoichiometry = {(1, 'A'): -1, (1, 'B'): -2, (1, 'C'): 0, (1, 'E'): 1, (1, 'G'): 0, (1, 'P'): 1,
                     (2, 'A'): -1, (2, 'B'): -1, (2, 'C'): 0, (2, 'E'): 0, (2, 'G'): 1, (2, 'P'): -1
                     }
    wocstr.stoichiometry = Param(wocstr.Reactions, wocstr.Components, initialize=stoichiometry)
    wocstr.Fb0 = Var(initialize=4.77)
    wocstr.Fb0.fixed = True
    wocstr.Tr0 = Var(initialize=78.2)
    wocstr.Tr0.fixed = True
    wocstr.Fb = Var(initialize=4, within=NonNegativeReals,bounds=(3,6))  # kg/s
    wocstr.Fb.fixed = True
    wocstr.Tr = Var(initialize=90,bounds=(70,100))  # Celsius degree
    wocstr.Tr.fixed = True
    wocstr.Fr = Var(initialize=5, within=NonNegativeReals)  # kg/s
    wocstr.XFr = Var(wocstr.Components, bounds=(0, 1))
    wocstr.K = Var(wocstr.Reactions, initialize=0.1)
    wocstr.phi_0th_term = Var(initialize=0)
    wocstr.phi_0th_term.fixed = True
    wocstr.phi_1st_Fb_term = Var(initialize=0)
    wocstr.phi_1st_Fb_term.fixed = True
    wocstr.phi_1st_Tr_term = Var(initialize=0)
    wocstr.phi_1st_Tr_term.fixed = True
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
    def profit_no_modifiers(m):
        return 1143.38 * m.Fr * m.XFr['P'] + 25.92 * m.Fr * m.XFr['E'] - 76.23 * m.Fa - 114.34 * m.Fb
    wocstr.profit_no_modifiers = Expression(rule=profit_no_modifiers)
    def profit(m):
        return 1143.38 * m.Fr * m.XFr['P'] + 25.92 * m.Fr * m.XFr['E'] - 76.23 * m.Fa - 114.34 * m.Fb+ m.phi_0th_term + m.phi_1st_Fb_term*(m.Fb-m.Fb0)+ m.phi_1st_Tr_term*(m.Tr-m.Tr0)
    wocstr.profit = Expression(rule=profit)
    def obj(m):
        return m.profit
    wocstr.objective = Objective(rule=obj, sense=maximize)
    return wocstr

def model_simulate(simulator, solver, Fb, Tr):
    simulator.Fb.fixed = True
    simulator.Tr.fixed = True
    simulator.Fb = Fb
    simulator.Tr = Tr
    solver.solve(simulator, tee=False)

solver = SolverFactory('ipopt', executable=r"../external/bin/ipopt.exe")

wocstr = plant_simulator()
plant_simulate(wocstr, solver, Fb=5, Tr=80)
init_value.to_template(wocstr, "wo_reactor_plant_init.txt")

# wocstr_model = get_mismatched_model()
# model_simulate(wocstr_model, solver, Fb=5, Tr=80)
# init_value.to_template(wocstr_model, "wo_reactor_model_init.txt")