# coding=utf-8

from pyomo.environ import *
from pyomo.dae import *


class FBR:
    def __init__(self):
        pass

    def build(self, block, option):
        setattr(block, 'mu_x', Param(initialize=0.092))
        setattr(block, 'mu_p', Param(initialize=0.005))
        setattr(block, 'K_x', Param(initialize=0.15))
        setattr(block, 'K_p', Param(initialize=0.0002))
        setattr(block, 'K_I', Param(initialize=0.1))
        setattr(block, 'K_H', Param(initialize=0.04))
        setattr(block, 'Y_xs', Param(initialize=0.45))
        setattr(block, 'Y_ps', Param(initialize=0.9))
        setattr(block, 'm_x', Param(initialize=0.014))
        setattr(block, 'S_F', Param(initialize=600))
        setattr(block, 'Beta', Param(initialize=0.00062))

        setattr(block, 'X', Var(block.t))
        setattr(block, 'P', Var(block.t))
        setattr(block, 'S', Var(block.t))
        setattr(block, 'V', Var(block.t))
        setattr(block, 'F', Var(block.t))
        setattr(block, 'dXdt', DerivativeVar(block.X))
        setattr(block, 'dPdt', DerivativeVar(block.P))
        setattr(block, 'dSdt', DerivativeVar(block.S))
        setattr(block, 'dVdt', DerivativeVar(block.V))

        def _eq22(m, i):
            if i == 0:
                return Constraint.Skip
            return m.dXdt[i] == m.mu_x * m.S[i] * m.X[i] / (m.K_x * m.X[i] + m.S[i]) - m.X[i] * \
                   (m.F[i] - m.Beta * m.V[i]) / m.V[i]
        setattr(block, 'eq22', Constraint(block.t, rule=_eq22))

        def _eq23(m, i):
            if i == 0:
                return Constraint.Skip
            elif option == 'RTO':
                return m.dPdt[i] == m.mu_p * m.S[i] * m.X[i] / (m.K_p + m.S[i] + m.S[i] ** 2 / m.K_I) - m.P[i] * \
                   (m.F[i] - m.Beta * m.V[i]) / m.V[i]
            return m.dPdt[i] == m.mu_p * m.S[i] * m.X[i] / (m.K_p + m.S[i] + m.S[i] ** 2 / m.K_I) - m.P[i] * \
                   (m.F[i] - m.Beta * m.V[i]) / m.V[i] - m.K_H * m.P[i]
        setattr(block, 'eq23', Constraint(block.t, rule=_eq23))

        def _eq24(m, i):
            if i == 0:
                return Constraint.Skip
            return m.dSdt[i] == -m.mu_x * m.S[i] * m.X[i] / (m.Y_xs * (m.K_x * m.X[i] + m.S[i])) - m.mu_p * \
                   m.S[i] * m.X[i] / (m.Y_ps * (m.K_p + m.S[i] + m.S[i] ** 2 / m.K_I)) - m.m_x * m.X[i] + m.S_F * \
                   m.F[i] / m.V[i] - m.S[i] * (m.F[i] - m.Beta * m.V[i]) / m.V[i]
        setattr(block, 'eq24', Constraint(block.t, rule=_eq24))

        def _eq25(m, i):
            if i == 0:
                return Constraint.Skip
            return m.dVdt[i] == m.F[i] - m.Beta * m.V[i]
        setattr(block, 'eq25', Constraint(block.t, rule=_eq25))

    def initialize(self, block, option=None):
        for i in block.t:
            getattr(block, 'X')[i] = 0.1
            getattr(block, 'P')[i] = 0.0
            getattr(block, 'S')[i] = 6.0
            getattr(block, 'V')[i] = 100.0
            getattr(block, 'F')[i] = 0.11
            # getattr(block, 'dXdt')[i] = -0.01206845160154677
            # getattr(block, 'dPdt')[i] = -1.946472755523073e-06
            # getattr(block, 'dSdt')[i] = 0.6739444207652815
            # getattr(block, 'dVdt')[i] = 0.04743532983366011
