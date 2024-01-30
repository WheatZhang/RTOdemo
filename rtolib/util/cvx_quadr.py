from pyomo.environ import *
import numpy as np
import scipy.stats.qmc

def fit_convex_quadratic_model(X, y, eps_hessian=1e-4):
    '''

    :param X:
    :param y:
    :param eps_hessian:
    :return: fitted_func = c0 + np.dot(c,x) + np.dot(x,np.dot(A,x))
    '''

    dimension = X.shape[1]
    n_data = X.shape[0]

    model_lp = ConcreteModel()
    model_lp.Dimension = Set(initialize=range(dimension))
    model_lp.DataIndex = Set(initialize=range(n_data))
    model_lp.c0 = Var()
    model_lp.c = Var(model_lp.Dimension)
    model_lp.h = Var(model_lp.Dimension,model_lp.Dimension)
    model_lp.beta = Var(model_lp.Dimension,model_lp.Dimension)
    def q_expr(m,i):
        return m.c0+sum([m.c[d]*X[i,d] for d in m.Dimension])+\
            sum([m.h[d1,d2]*X[i,d1]*X[i,d2] for d1 in m.Dimension for d2 in m.Dimension if d2>d1])+\
            0.5*sum([m.h[d,d]*X[i,d]*X[i,d] for d in m.Dimension])
    model_lp.q = Expression(model_lp.DataIndex,rule=q_expr)
    def beta_con(m,d):
        return m.h[d,d] >= 2*sum([m.beta[d,d1] for d1 in m.Dimension if d1 > d])+ \
               sum([m.beta[d2, d] for d2 in m.Dimension if d2 < d]) + eps_hessian
    model_lp.beta_con = Constraint(model_lp.Dimension,rule=beta_con)
    def beta_def_con_1(m,d1,d2):
        if d2 <= d1:
            return Constraint.Skip
        return m.beta[d1,d2] >= m.h[d1,d2]
    model_lp.beta_def_con_1 = Constraint(model_lp.Dimension,model_lp.Dimension,rule=beta_def_con_1)
    def beta_def_con_2(m,d1,d2):
        if d2 <= d1:
            return Constraint.Skip
        return -m.beta[d1,d2] <= m.h[d1,d2]
    model_lp.beta_def_con_2 = Constraint(model_lp.Dimension,model_lp.Dimension,rule=beta_def_con_2)
    def y_con(m, i):
        return m.q[i] <= y[i]
    model_lp.y_con = Constraint(model_lp.DataIndex,rule=y_con)
    def obj(m):
        return sum([y[i]-m.q[i] for i in m.DataIndex])
    model_lp.obj = Objective(rule=obj,sense=minimize)

    solver = SolverFactory('gurobi')
    # solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    # solver = SolverFactory('ipopt', executable=solver_executable)
    solver.solve(model_lp, tee=True)

    Hessian = np.zeros((dimension,dimension))
    for i in range(dimension):
        for j in range(dimension):
            if j >= i:
                Hessian[i,j] = model_lp.h[i,j].value
            else:
                Hessian[i,j] = model_lp.h[j,i].value
    L0=np.linalg.cholesky(Hessian)

    model_nlp = ConcreteModel()
    model_nlp.Dimension = Set(initialize=range(dimension))
    model_nlp.DataIndex = Set(initialize=range(n_data))
    model_nlp.c0 = Var(initialize=model_lp.c0.value)
    model_nlp.c = Var(model_nlp.Dimension, initialize=lambda m, d: model_lp.c[d].value)
    model_nlp.l = Var(model_nlp.Dimension, model_nlp.Dimension)
    def h_rule(m, i, j):
        return sum([m.l[i, k] * m.l[j, k] for k in m.Dimension if k >= i and k >= j])+ \
               sum([m.l[i, k] * m.l[k, j] for k in m.Dimension if k >= i and k < j])+\
               sum([m.l[k, i] * m.l[k, j] for k in m.Dimension if k < i and k >= j])+\
                  sum([m.l[k, i] * m.l[k, j] for k in m.Dimension if k < i and k < j])
    model_nlp.h = Expression(model_nlp.Dimension, model_nlp.Dimension, rule=h_rule)

    def q_expr(m, i):
        return m.c0 + sum([m.c[d] * X[i, d] for d in m.Dimension]) + \
               sum([m.h[d1, d2] * X[i, d1] * X[i, d2] for d1 in m.Dimension for d2 in m.Dimension if d2 > d1]) + \
               0.5 * sum([m.h[d, d] * X[i, d] * X[i, d] for d in m.Dimension])
    model_nlp.q = Expression(model_nlp.DataIndex, rule=q_expr)

    def y_con(m, i):
        return m.q[i] <= y[i]
    model_nlp.y_con = Constraint(model_nlp.DataIndex,rule=y_con)
    def obj(m):
        return sum([y[i]-m.q[i] for i in m.DataIndex])
    model_nlp.obj = Objective(rule=obj,sense=minimize)

    for i in range(dimension):
        for j in range(dimension):
            model_nlp.l[i, j] = L0[i, j]

    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    solver = SolverFactory('ipopt', executable=solver_executable)
    solver.options['max_iter'] = 5000
    solver.solve(model_nlp, tee=True)

    c0 = model_nlp.c0.value
    c = np.zeros((dimension,))
    for i in range(dimension):
        c[i] = model_nlp.c[i].value
    Hessian = np.zeros((dimension, dimension))
    for i in range(dimension):
        for j in range(dimension):
            if j >= i:
                Hessian[i, j] = value(model_nlp.h[i, j])
            else:
                Hessian[i, j] = value(model_nlp.h[j, i])
    y_pred = np.zeros((n_data,))
    for i in range(n_data):
        y_pred[i] = value(model_nlp.q[i])
    return c0, c, 0.5*Hessian, y_pred

def get_hypercube_sampling(dimension, n_data, lb, ub, seed=1):
    from rtolib.util.misc import sort_samples
    engine = scipy.stats.qmc.LatinHypercube(d=dimension, seed=seed)
    samples = engine.random(n_data)
    samples = sort_samples(samples)
    sample_X = np.zeros((n_data,dimension))
    permutation = np.random.permutation(n_data)
    for i, s in enumerate(samples):
        for j in range(dimension):
            sample_X[permutation[i],j] = s[j] * (ub[j] - lb[j]) + lb[j]
    return sample_X