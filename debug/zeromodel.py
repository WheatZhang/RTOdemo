from pyomo.environ import *


def ipopt_bug():
    # problem 1
    model = ConcreteModel()
    model.e = Var()
    model.con = Constraint(rule=lambda m:m.e==2)

    solver = SolverFactory('ipopt', executable=r"../external/bin/ipopt.exe")
    solver.solve(model, tee=True)

    # problem 2
    model = ConcreteModel()
    model.e = Var(initialize=1)
    model.con = Constraint(rule=lambda m:m.e**2==2)

    solver = SolverFactory('ipopt', executable=r"../external/bin/ipopt.exe")
    solver.solve(model, tee=True)

    # IPOPT bug：
    # 上述两个问题求解结果EXIT: Restoration Failed!
