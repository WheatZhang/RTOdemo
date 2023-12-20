from pyomo.environ import *
from pyomo.dae import *
import init_value
import numpy
import pandas
import matplotlib.pyplot as plt

def get_plant_simulator():
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
    setattr(process, 'F_mv', Var(initialize=0.11, within=NonNegativeReals))
    setattr(process, 'dXdt', DerivativeVar(process.X))
    setattr(process, 'dPdt', DerivativeVar(process.P))
    setattr(process, 'dSdt', DerivativeVar(process.S))
    setattr(process, 'dVdt', DerivativeVar(process.V))

    process.F_homo=Var()
    process.F_homo.fixed=True
    process.S0_homo=Var()
    process.S0_homo.fixed=True

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

    # def _volume_con(m):
    #     return m.V[192] <= 120
    # setattr(process, 'volume_con', Constraint(rule=_volume_con))

    def _feed_con(m, i):
        return m.F[i] == m.F_mv

    setattr(process, 'feed_con', Constraint(process.t, rule=_feed_con))

    def _eco_obj(m):
        return m.P[192] * m.V[192]

    setattr(process, 'eco_obj', Expression(rule=_eco_obj))

    def _obj(m):
        return (((m.S[0] - m.S0_homo) / 50) ** 2 + ((m.F_mv - m.F_homo) / 0.1) ** 2) * 1000

    setattr(process, 'obj', Objective(rule=_obj))

    process.X[0].fix(0.1)
    process.P[0].fix(0)
    process.V[0].fix(100)

    discretizer = TransformationFactory('dae.collocation')
    discretizer.apply_to(process, nfe=192, ncp=3, scheme='LAGRANGE-RADAU')

    init_value.load_init_from_tpl_free_collocation(process, 'init_solved_mesh192.txt', ignore_init_mismatch=True)
    return process


def plant_simulate(process, F, S0):
    process.F_homo.fix(F)
    process.S0_homo.fix(S0)

    solver = SolverFactory('ipopt', executable=r"../../external/bin/ipopt.exe")
    solver.solve(process, tee=False)

    return value(process.F_mv),value(process.S[0]),value(process.eco_obj),value(process.V[192])


def model_simulate(F, S0):
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
               (m.F[i] - m.Beta * m.V[i]) / m.V[i]
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

    # def _volume_con(m):
    #     return m.V[192] <= 120
    # setattr(process, 'volume_con', Constraint(rule=_volume_con))

    def _feed_con(m, i):
        return m.F[i] == m.F_mv
    setattr(process, 'feed_con', Constraint(process.t, rule=_feed_con))

    def _eco_obj(m):
        return m.P[192] * m.V[192]

    setattr(process, 'eco_obj', Expression(rule=_eco_obj))

    def _obj(m):
        return (((m.S[0] - S0) / 50) ** 2 + ((m.F_mv - F) / 0.1) ** 2) * 1000

    setattr(process, 'obj', Objective(rule=_obj))

    process.X[0].fix(0.1)
    process.P[0].fix(0)
    process.V[0].fix(100)

    discretizer = TransformationFactory('dae.collocation')
    discretizer.apply_to(process, nfe=192, ncp=3, scheme='LAGRANGE-RADAU')

    init_value.load_init_from_template(process, 'init_solved_mesh192_model.txt', ignore_init_mismatch=True)

    solver = SolverFactory('ipopt', executable=r"../../external/bin/ipopt.exe")
    solver.solve(process, tee=False)

    return value(process.F_mv),value(process.S[0]),value(process.eco_obj),value(process.V[192])


def get_simulation_data_plant():
    process = get_plant_simulator()
    n=0
    data = pandas.DataFrame(columns=['F','S0','obj','con'])
    for F in numpy.linspace(start=0.2, stop=0.3, num=4):
        for S0 in numpy.linspace(start=30, stop=60, num=4):
            print(n)
            F, S0, obj, con = plant_simulate(process, F, S0)
            data = data.append({'F':F,'S0':S0,'obj':obj,'con':con}, ignore_index=True)
            n+=1
    data.to_csv("data/plant_data_4.txt", sep='\t')


def plot_contour(data_file, x_list, y_list):
    data_stored = pandas.read_csv(data_file, index_col=0, header=0, sep='\t')
    x_len = len(x_list)
    y_len = len(y_list)
    Z = numpy.zeros(shape=(x_len, y_len))

    for i in range(x_len):
        for j in range(y_len):
            Z[j, i] = data_stored.iloc[i*y_len+j].loc['obj']

    CS = plt.contour(x_list,y_list,Z)
    plt.clabel(CS, inline=True, fmt='%d')
    plt.show()


if __name__ == "__main__":
    get_simulation_data_plant()

    F = numpy.linspace(start=0.1, stop=0.3, num=4)
    S0 = numpy.linspace(start=0, stop=60, num=4)
    plot_contour("data/plant_data_4.txt", F, S0)