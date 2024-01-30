from rtolib.model.wo_reactor_con import RTO_Mismatched_WO_reactor as Wo_reactor
from rtolib.model.wo_reactor_con import default_WOR_description as wo_description
from rtolib.model.parallel_wo import RTO_Mismatched_Parallel_WO as Parallel_wo_reactor
from rtolib.model.parallel_wo import default_WOR_description as parallel_wo_description
from rtolib.model.asu17000 import default_ASU_description, ASU_Model, ASU_Plant
from rtolib.util.cvx_quadr import fit_convex_quadratic_model, get_hypercube_sampling
from rtolib.core.solve import PyomoSimulator
from pyomo.environ import SolverFactory
from rtolib.core.pyomo_model import ModelSolvingStatus
import numpy as np
import pickle
from rtolib.util.init_value import load_init_from_template

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


def get_asu_model_samples():
    model_simulator = PyomoSimulator(ASU_Model())
    model_simulator.build(default_ASU_description)
    solver_executable = r"F:\Research\ModularTask\wys_fit\WholeOpt\ipopt38.exe"
    solver = SolverFactory('ipopt', executable=solver_executable)
    default_options = {}
    default_options['tol'] = 1e-8  # this could be smaller if necessary
    default_options['bound_push'] = 1e-10  # this is also crucial, for both convergence and speed
    default_options['max_iter'] = 1000
    solver.options['linear_solver'] = 'ma57'
    # it is sensitive to "obj_scaling_factor"
    # default_options['obj_scaling_factor'] = 10
    model_simulator.set_solver(solver, tee=False, default_options=default_options)

    mvs = default_ASU_description.symbol_list['MV']
    dimension = len(mvs)
    n_data = 200
    shrinking_factor = 0.2
    lb = []
    ub = []
    for mv in mvs:
        lb_origin = default_ASU_description.bounds[mv][0]
        ub_origin = default_ASU_description.bounds[mv][1]
        center = (lb_origin+ub_origin)/2
        lb_new = center-(center-lb_origin)*shrinking_factor
        ub_new = center+(ub_origin-center)*shrinking_factor
        lb.append(lb_new)
        ub.append(ub_new)
    X = get_hypercube_sampling(dimension, n_data, lb, ub, seed=1)

    cvs = list(default_ASU_description.symbol_list['CON'])
    cvs.append(default_ASU_description.symbol_list['OBJ'])

    Y = np.zeros((n_data,len(cvs)), dtype=float)
    success_flag = np.zeros((n_data,), dtype=int)
    param_values = {}
    for k in model_simulator.pyomo_model.parameters.keys():
        param_values[k] = model_simulator.pyomo_model.default_value[k]
    for i in range(n_data):
        if i%1 == 0:
            print("simulating data %d"%i)
        input_dict = {}
        for j,mv in enumerate(mvs):
            input_dict[mv] = X[i,j]
        try:
            outputs, solve_status = model_simulator.simulate(input_dict, param_values=param_values,
                                                         use_homo=True)
        except Exception as e:
            try:
                raise AssertionError("No need to try this")
                outputs, solve_status = model_simulator.simulate(input_dict, param_values=param_values,
                                                       use_homo=False)
            except Exception as e:
                solve_status=False
        else:
            if solve_status == False or solve_status == ModelSolvingStatus.OPTIMIZATION_FAILED:
                try:
                    raise AssertionError("No need to try this")
                    outputs, solve_status = model_simulator.simulate(input_dict, param_values=param_values,
                                                                     use_homo=False)
                except Exception as e:
                    solve_status = False
        if solve_status == False or solve_status == ModelSolvingStatus.OPTIMIZATION_FAILED:
            print("solution fails")
        else:
            for j, cv in enumerate(cvs):
                success_flag[i] = 1
                Y[i,j] = outputs[cv]


    with open("data/asu_quadratic/X_model",'wb') as f:
        pickle.dump(X, f)
    with open("data/asu_quadratic/Y_model", 'wb') as f:
        pickle.dump(Y, f)
    with open("data/asu_quadratic/success_flag_model", 'wb') as f:
        pickle.dump(success_flag, f)


def get_asu_plant_samples():
    model_simulator = PyomoSimulator(ASU_Plant())
    model_simulator.build(default_ASU_description)
    solver_executable = r"F:\Research\ModularTask\wys_fit\WholeOpt\ipopt38.exe"
    solver = SolverFactory('ipopt', executable=solver_executable)
    default_options = {}
    default_options['tol'] = 1e-8  # this could be smaller if necessary
    default_options['bound_push'] = 1e-10  # this is also crucial, for both convergence and speed
    default_options['max_iter'] = 500
    solver.options['linear_solver'] = 'ma57'
    # it is sensitive to "obj_scaling_factor"
    # default_options['obj_scaling_factor'] = 10
    model_simulator.set_solver(solver, tee=False, default_options=default_options)

    mvs = default_ASU_description.symbol_list['MV']
    dimension = len(mvs)
    n_data = 2000
    shrinking_factor = 1
    lb = []
    ub = []
    for mv in mvs:
        lb_origin = default_ASU_description.bounds[mv][0]
        ub_origin = default_ASU_description.bounds[mv][1]
        center = (lb_origin+ub_origin)/2
        lb_new = center-(center-lb_origin)*shrinking_factor
        ub_new = center+(ub_origin-center)*shrinking_factor
        lb.append(lb_new)
        ub.append(ub_new)
    X = get_hypercube_sampling(dimension, n_data, lb, ub, seed=1)

    cvs = list(default_ASU_description.symbol_list['CON'])
    cvs.append(default_ASU_description.symbol_list['OBJ'])

    Y = np.zeros((n_data,len(cvs)), dtype=float)
    success_flag = np.zeros((n_data,), dtype=int)
    param_values = {}
    for k in model_simulator.pyomo_model.parameters.keys():
        param_values[k] = model_simulator.pyomo_model.default_value[k]
    for i in range(n_data):
        if i%10 == 0:
            print("simulating data %d"%i)
        input_dict = {}
        for j,mv in enumerate(mvs):
            input_dict[mv] = X[i,j]
        try:
            outputs, solve_status = model_simulator.simulate(input_dict, param_values=param_values,
                                                         use_homo=True)
        except Exception as e:
            load_init_from_template(model_simulator.model, model_simulator.pyomo_model.initial_value_file, \
                                               ignore_init_mismatch=True)
            solve_status=False
        else:
            if solve_status == False or solve_status == ModelSolvingStatus.OPTIMIZATION_FAILED:
                load_init_from_template(model_simulator.model, model_simulator.pyomo_model.initial_value_file, \
                                        ignore_init_mismatch=True)
                solve_status = False
        if solve_status == False or solve_status == ModelSolvingStatus.OPTIMIZATION_FAILED:
            print("solution fails")
        else:
            for j, cv in enumerate(cvs):
                success_flag[i] = 1
                Y[i,j] = outputs[cv]


        with open("data/asu_quadratic/X_plant",'wb') as f:
            pickle.dump(X[0:i,:], f)
        with open("data/asu_quadratic/Y_plant", 'wb') as f:
            pickle.dump(Y[0:i,:], f)
        with open("data/asu_quadratic/success_flag_plant", 'wb') as f:
            pickle.dump(success_flag[0:i], f)

def fit_asu_convex_quadratic_surrogate():
    with open("data/asu_quadratic/X",'rb') as f:
        X=pickle.load(f)
    with open("data/asu_quadratic/Y", 'rb') as f:
        Y=pickle.load(f)
    with open("data/asu_quadratic/success_flag", 'rb') as f:
        success_flag = pickle.load(f)
    mvs = parallel_wo_description.symbol_list['MV']
    cvs = list(default_ASU_description.symbol_list['CON'])
    cvs.append(default_ASU_description.symbol_list['OBJ'])
    X = X[success_flag==1,:]
    Y = Y[success_flag==1,:]
    mv_strings = ['m.FeedMFlow',
                  'm.FeedSplitterOut1Ratio',
                  'm.FeedSplitterOut2Ratio',
                  'm.HPCCondPrdtMFlow',
                  'm.LPCExt46MFlow',
                  'm.LPCExt15MFlow',
                  'm.ASCCondRefRatio',
                  'm.OxSplitterOutARatio']
    cv_strings = {
        "cost": "OBJ",
        "OxygenPrdtPurityCon": "PurityCon1",
        "OxygenPrdtImpurityCon": "PurityCon2",
        "NitrogenPrdtPurityCon": "PurityCon3",
        "NitrogenPrdtImpurityCon": "PurityCon4",
        "ArgonPrdtImpurityCon": "PurityCon5",
        "CrudeArgonImpurityCon": "PurityCon6",
        "LOX_Con": "ProductCon1",
        "GAN_Con": "ProductCon2",
        "GAR_Con": "ProductCon3",
        "MainCoolerTempDiffCon": "TempDiffCon1",
        "HeatLPCTempDiffCon": "TempDiffCon2",
        "HPA_GOX_LB_Con": "TurbineCon1",
        "HPA_GOX_UB_Con": "TurbineCon2",
        "Feed_LB_Con": "TurbineCon3",
        "Feed_UB_Con": "TurbineCon4",
        "HPA_LB_Con": "TurbineCon7",
        "HPA_UB_Con": "TurbineCon8",
        "TA_LB_Con": "TurbineCon9",
        "TA_UB_Con": "TurbineCon10",
        "OverallHeatCon": "OverallHeatBlnc",
    }
    with open("data/asu_quadratic_model_text.txt", 'w') as fp:
        for j, cv in enumerate(cvs):
            fp.write("def "+cv_strings[cv]+"(m):\n\treturn ")
            c, b, A, _ = fit_convex_quadratic_model(X, Y[:, j])
            print("fitting %s" % cv)
            print("c=")
            print(c)
            print("b=")
            print(b)
            print("A=")
            print(A)

            model_string = "%.6e"%c
            for k1, mv1 in enumerate(mvs):
                model_string += "\\\n+%.6e" % b[k1] +"*"+ mv_strings[k1]
            for k1, mv1 in enumerate(mvs):
                for k2, mv2 in enumerate(mvs):
                    model_string+="\\\n+%.6e"%A[k1,k2]+"*"+mv_strings[k1]+"*"+mv_strings[k2]
            fp.write(model_string)
            fp.write("\ntest_model."+cv_strings[cv]+" = Expression(rule="+cv_strings[cv]+")\n\n")


if __name__ == "__main__":
    # fit_wo_convex_quadratic_surrogate()
    # fit_parallel_wo_convex_quadratic_surrogate()
    get_asu_plant_samples()
    # fit_asu_convex_quadratic_surrogate()
