#!/usr/bin/env python
#-*- coding:utf-8 -*-
from pyomo.environ import *
from pyomo.dae import *
import rtolib.util.init_value as ivt

def get_model():
    '''
    Reference: Hille, R. and Budman, H.M., 2018. Simultaneous identification and optimization of biochemical processes
     under model-plant mismatch using output uncertainty bounds. Computers & Chemical Engineering, 113, pp.125-138.
    :return:
    '''
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

    model.X = Var(model.t,initialize=0)
    model.P = Var(model.t,initialize=0)
    model.S = Var(model.t,initialize=20)
    model.V = Var(model.t,initialize=100)

    model.dXdt = DerivativeVar(model.X)
    model.dPdt = DerivativeVar(model.P)
    model.dSdt = DerivativeVar(model.S)
    model.dVdt = DerivativeVar(model.V)

    model.F = Var(initialize=0.04)
    model.S0 = Var(initialize=0.1)

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

    def _terminal_con(m):
        return m.V[200] <= 120
    model.terminal_con = Constraint(rule=_terminal_con)

    def obj(m):
        return -m.P[200]
    model.obj = Objective(rule=obj)

    discretizer = TransformationFactory('dae.collocation')
    discretizer.apply_to(model, nfe=10, ncp=3, scheme='LAGRANGE-RADAU')

    return model

def fix_input_variable(model):
    model.F.fix()
    model.S0.fix()

def disable_inequality_con(model):
    model.terminal_con.deactivate()

def disable_obj(model):
    model.obj.deactivate()

def enable_all(model):
    model.F.fixed = False
    model.S0.fixed = False
    model.terminal_con.activate()
    model.obj.activate()

def get_solver():
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    solver = SolverFactory('ipopt', executable=solver_executable)
    solver.options['max_iter']=500
    return solver

def main():
    '''
    simulation problem solved as a simulation problem: fails
    simulation problem solved as an optimization problem: successes
    Conclution: the objective function or the ineq constraint function is important, even for the simulation problem
                Possible explainations: it gives good bounds during the iteraiton
    solve the optimization problem directly: fails
    Hypothesis: the input variable is better to be fixed at the early stage of the solving process, to prevent going
                to a bad value
    increase the number of finite elements: all problems fail
    Conclution: such skills are only of limited help
    :return:
    '''
    solver = get_solver()

    # simulation problem solved as a simulation problem
    print("\n=============================\nsimulation problem solved as a simulation problem\n=============================")
    model = get_model()
    fix_input_variable(model)
    disable_inequality_con(model)
    disable_obj(model)
    solver.solve(model, tee=True)

    # simulation problem solved as an optimization problem
    print(
        "\n=============================\nsimulation problem solved as an optimization problem 1\n=============================")
    model = get_model()
    fix_input_variable(model)
    disable_obj(model)
    solver.solve(model, tee=True)

    print(
        "\n=============================\nsimulation problem solved as an optimization problem 2\n=============================")
    model = get_model()
    fix_input_variable(model)
    disable_inequality_con(model)
    solver.solve(model, tee=True)

    # optimization problem solved directly
    print(
        "\n=============================\nsolve the optimization problem directly\n=============================")
    model = get_model()
    solver.solve(model, tee=True)

    # optimization problem solved with good initial value
    print(
        "\n=============================\nsolve the optimization problem with good initial value\n=============================")
    model = get_model()
    fix_input_variable(model)
    disable_inequality_con(model)
    solver.solve(model, tee=True)
    enable_all(model)
    solver.solve(model, tee=True)
    ivt.to_template(model, "optimization_result.txt")

if __name__ == "__main__":
    main()