from rtolib.model.wo_reactor_con import RTO_Mismatched_WO_reactor as Wo_reactor
from rtolib.model.wo_reactor_con import default_WOR_description as wo_description
from rtolib.model.parallel_wo import RTO_Mismatched_Parallel_WO as Parallel_wo_reactor
from rtolib.model.parallel_wo import default_WOR_description as parallel_wo_description
from rtolib.util.cvx_quadr import fit_convex_quadratic_model, get_hypercube_sampling
from rtolib.core.solve import PyomoSimulator
from pyomo.environ import SolverFactory
import numpy as np

def fit_wo_convex_quadratic_surrogate():
    model_simulator = PyomoSimulator(Wo_reactor())
    model_simulator.build(wo_description)
    default_options = {'max_iter': 100}
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    solver = SolverFactory('ipopt', executable=solver_executable)
    model_simulator.set_solver(solver, tee=False, default_options=default_options)

    mvs = wo_description.symbol_list['MV']
    dimension = len(mvs)
    n_data = 100
    lb = []
    ub = []
    for mv in mvs:
        lb.append(wo_description.bounds[mv][0])
        ub.append(wo_description.bounds[mv][1])
    X = get_hypercube_sampling(dimension, n_data, lb, ub, seed=1)

    cvs=['cost', 'con']
    Y = np.zeros((n_data,len(cvs)), dtype=float)
    param_values = {}
    for k in model_simulator.pyomo_model.parameters.keys():
        param_values[k] = model_simulator.pyomo_model.default_value[k]
    for i in range(n_data):
        input_dict = {}
        for j,mv in enumerate(mvs):
            input_dict[mv] = X[i,j]
        outputs, solve_status = model_simulator.simulate(input_dict, param_values=param_values,
                                                         use_homo=False)
        for j, cv in enumerate(cvs):
            Y[i,j] = outputs[cv]

    mv_strings = ['m.Fa', 'm.Fb', 'm.Tr']
    for j, cv in enumerate(cvs):
        c, b, A, _ = fit_convex_quadratic_model(X, Y[:, j])
        print("fitting %s" % cv)
        print("c=")
        print(c)
        print("b=")
        print(b)
        print("A=")
        print(A)

        model_string = "%.6e"%c+"\\\n"
        for k1, mv1 in enumerate(mvs):
            model_string += "+%.6e" % b[k1] +"*"+ mv_strings[k1]+"\\\n"
        for k1, mv1 in enumerate(mvs):
            for k2, mv2 in enumerate(mvs):
                model_string+="+%.6e"%A[k1,k2]+"*"+mv_strings[k1]+"*"+mv_strings[k2]+"\\\n"
        print(model_string)

def fit_parallel_wo_convex_quadratic_surrogate():
    model_simulator = PyomoSimulator(Parallel_wo_reactor())
    model_simulator.build(parallel_wo_description)
    default_options = {'max_iter': 200}
    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    solver = SolverFactory('ipopt', executable=solver_executable)
    model_simulator.set_solver(solver, tee=False, default_options=default_options)

    mvs = parallel_wo_description.symbol_list['MV']
    dimension = len(mvs)
    n_data = 1000
    lb = []
    ub = []
    for mv in mvs:
        lb.append(parallel_wo_description.bounds[mv][0])
        ub.append(parallel_wo_description.bounds[mv][1])
    X = get_hypercube_sampling(dimension, n_data, lb, ub, seed=1)

    cvs=['cost', 'con']
    Y = np.zeros((n_data,len(cvs)), dtype=float)
    param_values = {}
    for k in model_simulator.pyomo_model.parameters.keys():
        param_values[k] = model_simulator.pyomo_model.default_value[k]
    for i in range(n_data):
        if i%100 == 0:
            print("simulating data %d"%i)
        input_dict = {}
        for j,mv in enumerate(mvs):
            input_dict[mv] = X[i,j]
        outputs, solve_status = model_simulator.simulate(input_dict, param_values=param_values,
                                                         use_homo=False)
        for j, cv in enumerate(cvs):
            Y[i,j] = outputs[cv]

    mv_strings = ['m.reactor_block[1].Fb', 'm.reactor_block[2].Fb', 'm.reactor_block[3].Fb',\
                  'm.reactor_block[1].Tr', 'm.reactor_block[2].Tr', 'm.reactor_block[3].Tr']
    for j, cv in enumerate(cvs):
        c, b, A, _ = fit_convex_quadratic_model(X, Y[:, j])
        print("fitting %s" % cv)
        print("c=")
        print(c)
        print("b=")
        print(b)
        print("A=")
        print(A)

        model_string = "%.6e"%c+"\\\n"
        for k1, mv1 in enumerate(mvs):
            model_string += "+%.6e" % b[k1] +"*"+ mv_strings[k1]+"\\\n"
        for k1, mv1 in enumerate(mvs):
            for k2, mv2 in enumerate(mvs):
                model_string+="+%.6e"%A[k1,k2]+"*"+mv_strings[k1]+"*"+mv_strings[k2]+"\\\n"
        print(model_string)

if __name__ == "__main__":
    # fit_wo_convex_quadratic_surrogate()
    fit_parallel_wo_convex_quadratic_surrogate()
