#!/usr/bin/env python
#-*- coding:utf-8 -*-
from pyomo.environ import *
from pyomo.dae import *
import rtolib.util.init_value as ivt
import pyomo
from bisect import bisect_right, bisect_left
import copy

def get_init_val_through_optimization():
    model = ConcreteModel()

    model.t = ContinuousSet(bounds=[0, 200])
    model.mu_x = Param(initialize=0.092)
    model.K_x = Param(initialize=0.15)
    model.mu_P = Param(initialize=0.005)
    model.K_P = Param(initialize=0.0002)
    model.K_I = Param(initialize=0.1)
    model.K_H = Param(initialize=0.04)
    model.Y_XS = Param(initialize=0.45)
    model.Y_PS = Param(initialize=0.9)
    model.m_X = Param(initialize=0.014)
    model.s_f = Param(initialize=600)

    model.X = Var(model.t, initialize=0, within=NonNegativeReals)
    model.P = Var(model.t, initialize=0, within=NonNegativeReals)
    model.S = Var(model.t, initialize=20, within=NonNegativeReals)
    model.V = Var(model.t, initialize=100, within=NonNegativeReals)

    model.dXdt = DerivativeVar(model.X)
    model.dPdt = DerivativeVar(model.P)
    model.dSdt = DerivativeVar(model.S)
    model.dVdt = DerivativeVar(model.V)

    model.F = Var(initialize=0.04)
    model.F.fix()
    model.S0 = Var(initialize=0.1)
    model.S0.fix()

    def _expr_1(m, t):
        return m.mu_x * m.S[t] * m.X[t] / (m.K_x * m.X[t] + m.S[t])

    model.expr_1 = Expression(model.t, rule=_expr_1)

    def _expr_2(m, t):
        return m.mu_P * m.S[t] * m.X[t] / (m.K_P + m.S[t] + m.S[t] * m.S[t] / m.K_I)

    model.expr_2 = Expression(model.t, rule=_expr_2)

    def _equ_1(m, t):
        if t == 0:
            return m.X[0] == 0.1
        else:
            return m.dXdt[t] == m.expr_1[t] - m.X[t] / m.V[t] * m.dVdt[t]

    model.equ_1 = Constraint(model.t, rule=_equ_1)

    def _equ_2(m, t):
        if t == 0:
            return m.P[0] == 0
        else:
            return m.dPdt[t] == m.expr_2[t] - m.K_H * m.P[t] - m.P[t] / m.V[t] * m.dVdt[t]

    model.equ_2 = Constraint(model.t, rule=_equ_2)

    def _equ_3(m, t):
        if t == 0:
            return m.S[0] == m.S0
        else:
            return m.dSdt[t] == -m.expr_1[t] / m.Y_XS - m.expr_2[t] / m.Y_PS - m.m_X * m.X[t] + m.F * m.s_f / m.V[t] - \
                m.S[t] / m.V[t] * m.dVdt[t]

    model.equ_3 = Constraint(model.t, rule=_equ_3)

    def _equ_4(m, t):
        if t == 0:
            return m.V[0] == 100
        else:
            return m.dVdt[t] == m.F - 6.226e-4 * m.V[t]

    model.equ_4 = Constraint(model.t, rule=_equ_4)

    # def _terminal_con(m):
    #     return m.V[200] <= 120
    # model.terminal_con = Constraint(rule=_terminal_con)

    def obj(m):
        return -m.P[200]

    model.obj = Objective(rule=obj)

    discretizer = TransformationFactory('dae.collocation')
    discretizer.apply_to(model, nfe=20, ncp=3, scheme='LAGRANGE-RADAU')

    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    solver = SolverFactory('ipopt', executable=solver_executable)
    results = solver.solve(model, tee=True)
    ivt.to_template(model, "init_value/opt_init.txt")

def get_simulation_result():
    model = ConcreteModel()

    model.t = ContinuousSet(bounds=[0,200])
    model.mu_x = Param(initialize=0.092)
    model.K_x = Param(initialize=0.15)
    model.mu_P=Param(initialize=0.005)
    model.K_P=Param(initialize=0.0002)
    model.K_I=Param(initialize=0.1)
    model.K_H=Param(initialize=0.04)
    model.Y_XS=Param(initialize=0.45)
    model.Y_PS=Param(initialize=0.9)
    model.m_X=Param(initialize=0.014)
    model.s_f=Param(initialize=600)

    model.X = Var(model.t, initialize=0, within=NonNegativeReals)
    model.P = Var(model.t, initialize=0, within=NonNegativeReals)
    model.S = Var(model.t, initialize=20, within=NonNegativeReals)
    model.V = Var(model.t, initialize=100, within=NonNegativeReals)

    model.dXdt = DerivativeVar(model.X)
    model.dPdt = DerivativeVar(model.P)
    model.dSdt = DerivativeVar(model.S)
    model.dVdt = DerivativeVar(model.V)

    model.F = Var(initialize=0.1728)
    model.F.fix()
    model.S0 = Var(initialize=54.72)
    model.S0.fix()

    def _expr_1(m, t):
        return m.mu_x*m.S[t]*m.X[t]/(m.K_x*m.X[t]+m.S[t])
    model.expr_1 = Expression(model.t, rule=_expr_1)
    def _expr_2(m, t):
        return m.mu_P*m.S[t]*m.X[t]/(m.K_P+m.S[t]+m.S[t]*m.S[t]/m.K_I)
    model.expr_2 = Expression(model.t, rule=_expr_2)

    def _equ_1(m, t):
        if t == 0:
            return m.X[0] == 0.1
        else:
            return m.dXdt[t] == m.expr_1[t] -m.X[t]/m.V[t]*m.dVdt[t]
    model.equ_1 = Constraint(model.t, rule=_equ_1)
    def _equ_2(m, t):
        if t == 0:
            return m.P[0] == 0
        else:
            return m.dPdt[t] == m.expr_2[t] -m.K_H*m.P[t]-m.P[t]/m.V[t]*m.dVdt[t]
    model.equ_2 = Constraint(model.t, rule=_equ_2)
    def _equ_3(m, t):
        if t == 0:
            return m.S[0] == m.S0
        else:
            return m.dSdt[t] == -m.expr_1[t]/m.Y_XS-m.expr_2[t]/m.Y_PS -m.m_X*m.X[t]+m.F*m.s_f/m.V[t]-m.S[t]/m.V[t]*m.dVdt[t]
    model.equ_3 = Constraint(model.t, rule=_equ_3)
    def _equ_4(m, t):
        if t == 0:
            return m.V[0] == 100
        else:
            return m.dVdt[t]==m.F-6.226e-4*m.V[t]
    model.equ_4 = Constraint(model.t, rule=_equ_4)

    discretizer = TransformationFactory('dae.collocation')
    discretizer.apply_to(model,nfe=200,ncp=3,scheme='LAGRANGE-RADAU')


    init_value_old = {}
    for var_name, dict, index in ivt.load_general_txt_init("init_value/init_from_ode.txt"):
        init_value_old[var_name] = dict
    time_old=[]
    for t in init_value_old["X"].keys():
        time_old.append(t)
    time_old.sort()

    all_attr = dir(model)
    for var_name in all_attr:
        if var_name in init_value_old.keys():
            if (isinstance(getattr(model, var_name),pyomo.core.base.var.IndexedVar) and \
                getattr(model, var_name)._index.name == "t") or \
                    (isinstance(getattr(model, var_name),pyomo.dae.DerivativeVar)):
                print(var_name)
                index = getattr(model, var_name)._index
                for idx in index:
                    if isinstance(idx, int) or isinstance(idx, float):
                        old_time_index = bisect_left(time_old, idx)
                        if old_time_index >= len(time_old):
                            old_time_index = len(time_old)-1
                        idx_old = time_old[old_time_index]
                    else:
                        new_t=idx[0]
                        old_time_index=bisect_left(time_old, new_t)
                        if old_time_index >= len(time_old):
                            old_time_index = len(time_old)-1
                        old_t=time_old[old_time_index]
                        idx_old=copy(idx)
                        idx_old[0]=old_t
                    getattr(model, var_name)[idx]=init_value_old[var_name][idx_old]
            else:
                ivt.set_initials(model, var_name, init_value_old[var_name])

    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    solver = SolverFactory('ipopt', executable=solver_executable)
    results = solver.solve(model, tee=True)
    ivt.to_template(model, "init_value/simulation50.txt")

def get_optimization_result():
    model = ConcreteModel()

    model.t = ContinuousSet(bounds=[0, 200])
    model.mu_x = Param(initialize=0.092)
    model.K_x = Param(initialize=0.15)
    model.mu_P = Param(initialize=0.005)
    model.K_P = Param(initialize=0.0002)
    model.K_I = Param(initialize=0.1)
    model.K_H = Param(initialize=0.04)
    model.Y_XS = Param(initialize=0.45)
    model.Y_PS = Param(initialize=0.9)
    model.m_X = Param(initialize=0.014)
    model.s_f = Param(initialize=600)

    model.X = Var(model.t, initialize=0, bounds=(0,100))
    model.P = Var(model.t, initialize=0, bounds=(0,10))
    model.S = Var(model.t, initialize=20, bounds=(0,200))
    model.V = Var(model.t, initialize=100, bounds=(0,150))

    model.dXdt = DerivativeVar(model.X)
    model.dPdt = DerivativeVar(model.P)
    model.dSdt = DerivativeVar(model.S)
    model.dVdt = DerivativeVar(model.V)

    model.F = Var(initialize=0.1728, bounds=(0,0.4))
    model.S0 = Var(initialize=54.72, bounds=(0,100))

    def _expr_1(m, t):
        return m.mu_x * m.S[t] * m.X[t] / (m.K_x * m.X[t] + m.S[t])

    model.expr_1 = Expression(model.t, rule=_expr_1)

    def _expr_2(m, t):
        return m.mu_P * m.S[t] * m.X[t] / (m.K_P + m.S[t] + m.S[t] * m.S[t] / m.K_I)

    model.expr_2 = Expression(model.t, rule=_expr_2)

    def _equ_1(m, t):
        if t == 0:
            return m.X[0] == 0.1
        else:
            return m.dXdt[t] == m.expr_1[t] - m.X[t] / m.V[t] * m.dVdt[t]

    model.equ_1 = Constraint(model.t, rule=_equ_1)

    def _equ_2(m, t):
        if t == 0:
            return m.P[0] == 0
        else:
            return m.dPdt[t] == m.expr_2[t] - m.K_H * m.P[t] - m.P[t] / m.V[t] * m.dVdt[t]

    model.equ_2 = Constraint(model.t, rule=_equ_2)

    def _equ_3(m, t):
        if t == 0:
            return m.S[0] == m.S0
        else:
            return m.dSdt[t] == -m.expr_1[t] / m.Y_XS - m.expr_2[t] / m.Y_PS - m.m_X * m.X[t] + m.F * m.s_f / m.V[t] - \
                m.S[t] / m.V[t] * m.dVdt[t]

    model.equ_3 = Constraint(model.t, rule=_equ_3)

    def _equ_4(m, t):
        if t == 0:
            return m.V[0] == 100
        else:
            return m.dVdt[t] == m.F - 6.226e-4 * m.V[t]

    model.equ_4 = Constraint(model.t, rule=_equ_4)

    def _terminal_con(m):
        return m.V[200] <= 120
    model.terminal_con = Constraint(rule=_terminal_con)

    def obj(m):
        return -m.P[200]*1e6
    model.obj = Objective(rule=obj)

    discretizer = TransformationFactory('dae.collocation')
    discretizer.apply_to(model, nfe=200, ncp=3, scheme='LAGRANGE-RADAU')

    ivt.load_init_from_template(model, "init_value/simulation50.txt")

    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    solver = SolverFactory('ipopt', executable=solver_executable)
    results = solver.solve(model, tee=True)
    solver.options['mu_strategy']="adaptive"
    ivt.to_template(model, "init_value/optimization50.txt")

if __name__ == "__main__":
    # get_init_val_through_optimization()
    # get_simulation_result()
    get_optimization_result()