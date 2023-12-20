from pyomoModel import get_full_asu_model
from pyomo.environ import *
import rtolib.util.init_value as ivt

def norminal_simulation():
    model = get_full_asu_model()
    ivt.load_init_from_template(model, "GOX17000InitValue.txt")

    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    solver = SolverFactory('ipopt', executable=solver_executable)
    results = solver.solve(model, tee=True)
    ivt.to_template(model, "GOX17000SimInitValue.txt")

def get_optimization_result(model):
    ret = {}
    ret['FeedMFlow'] = value(model.FeedMFlow)
    ret['FeedSplitterOut1Ratio'] = value(model.FeedSplitterOut1Ratio) #TA Frac
    ret['FeedSplitterOut2Ratio'] = value(model.FeedSplitterOut2Ratio) #MA_Frac
    ret['HPCCondPrdtMFlow'] = value(model.HPCCondPrdtMFlow) #40-LIN
    ret['LPCExt46MFlow']	 = value(model.LPCExt46MFlow) #52-WN
    ret['LPCExt15MFlow']	 = value(model.LPCExt15MFlow) #47-ARC
    ret['ASCCondRefRatio']	 = value(model.ASCCondRefRatio) #ASC_Reflux
    ret['OxSplitterOutARatio'] = value(model.OxSplitterOutARatio) #LOX_Frac
    ret['OBJ'] = value(model.OBJ) # objective function
    return ret

def unfix_mvs(model):
    model.FeedMFlow.fixed=False
    model.FeedSplitterOut1Ratio.fixed=False
    model.FeedSplitterOut2Ratio.fixed=False
    model.HPCCondPrdtMFlow.fixed=False
    model.LPCExt46MFlow.fixed=False
    model.LPCExt15MFlow.fixed=False
    model.ASCCondRefRatio.fixed=False
    model.OxSplitterOutARatio.fixed=False


def nominal_optimization():
    model = get_full_asu_model(build_constraint=True)
    unfix_mvs(model)
    ivt.load_init_from_template(model, "OriSSOFilled.txt")

    # solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    solver_executable = r"F:\Research\ModularTask\wys_fit\WholeOpt\ipopt38.exe"
    solver = SolverFactory('ipopt', executable=solver_executable)
    solver.options['tol']=1e-8  # this could be smaller if necessary
    solver.options['bound_push']=1e-8   # this is also crucial, for both convergence and speed
    # good value = 1e-10
    solver.options['linear_solver']='ma57'  # mumps fails
    solver.solve(model, tee=True)
    result = get_optimization_result(model)
    print(result)



if __name__ == "__main__":
    # norminal_simulation()
    nominal_optimization()