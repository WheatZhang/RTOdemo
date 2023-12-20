from pyomo.environ import *
from pyomo.dae import *
import init_value

process = ConcreteModel()
process.tf = Param(initialize=192)
process.t = ContinuousSet(bounds=(0, process.tf))
setattr(process, 'mu_x', Param(initialize=0.092))
setattr(process, 'mu_p', Param(initialize=0.005))
setattr(process, 'K_x', Param(initialize=0.15))
setattr(process, 'K_p', Param(initialize=0.0002))
setattr(process, 'K_I', Param(initialize=0.1))
setattr(process, 'K_H', Param(initialize=0.04))
setattr(process, 'Y_xs', Param(initialize=0.45))
setattr(process, 'Y_ps', Param(initialize=0.9))
setattr(process, 'm_x', Param(initialize=0.014))
setattr(process, 'S_F', Param(initialize=600))
setattr(process, 'Beta', Param(initialize=0.00062))

setattr(process, 'X', Var(process.t, within=NonNegativeReals))
setattr(process, 'P', Var(process.t, within=NonNegativeReals))
setattr(process, 'S', Var(process.t, within=NonNegativeReals))
setattr(process, 'V', Var(process.t, within=NonNegativeReals))
setattr(process, 'F', Var(process.t))
setattr(process, 'F_mv', Var(initialize=0.11,within=NonNegativeReals))
setattr(process, 'dXdt', DerivativeVar(process.X))
setattr(process, 'dPdt', DerivativeVar(process.P))
setattr(process, 'dSdt', DerivativeVar(process.S))
setattr(process, 'dVdt', DerivativeVar(process.V))

def _eq22(m, i):
    if i == 0:
        return Constraint.Skip
    return m.dXdt[i] == m.mu_x * m.S[i] * m.X[i] / (m.K_x * m.X[i] + m.S[i]) - m.X[i] * \
           (m.F[i] - m.Beta * m.V[i]) / m.V[i]
setattr(process, 'eq22', Constraint(process.t, rule=_eq22))

def _eq23(m, i):
    if i == 0:
        return Constraint.Skip
    return m.dPdt[i] == m.mu_p * m.S[i] * m.X[i] / (m.K_p + m.S[i] + m.S[i] ** 2 / m.K_I) - m.P[i] * \
           (m.F[i] - m.Beta * m.V[i]) / m.V[i] - m.K_H * m.P[i]
setattr(process, 'eq23', Constraint(process.t, rule=_eq23))

def _eq24(m, i):
    if i == 0:
        return Constraint.Skip
    return m.dSdt[i] == -m.mu_x * m.S[i] * m.X[i] / (m.Y_xs * (m.K_x * m.X[i] + m.S[i])) - m.mu_p * \
           m.S[i] * m.X[i] / (m.Y_ps * (m.K_p + m.S[i] + m.S[i] ** 2 / m.K_I)) - m.m_x * m.X[i] + m.S_F * \
           m.F[i] / m.V[i] - m.S[i] * (m.F[i] - m.Beta * m.V[i]) / m.V[i]
setattr(process, 'eq24', Constraint(process.t, rule=_eq24))

def _eq25(m, i):
    if i == 0:
        return Constraint.Skip
    return m.dVdt[i] == m.F[i] - m.Beta * m.V[i]
setattr(process, 'eq25', Constraint(process.t, rule=_eq25))

def _volume_con(m):
    return m.V[192] <= 120
setattr(process, 'volume_con', Constraint(rule=_volume_con))

def _feed_con(m, i):
    return m.F[i] == m.F_mv
setattr(process, 'feed_con', Constraint(process.t, rule=_feed_con))

def _obj(m):
    return -m.P[192]
setattr(process, 'obj', Objective(rule=_obj))

process.X[0].fix(0.1)
process.P[0].fix(0)
process.V[0].fix(100)

discretizer = TransformationFactory('dae.collocation')
discretizer.apply_to(process, nfe=192/6, ncp=3, scheme='LAGRANGE-RADAU')

# init_value.to_template(process, 'init_template.txt')

init_value.load_init_from_template(process, 'init_solved_plant.txt', ignore_init_mismatch=True)


solver = SolverFactory('ipopt', executable=r"../../external/bin/ipopt.exe")



solver.solve(process, tee=True)

init_value.to_template(process, 'init_solved.txt')