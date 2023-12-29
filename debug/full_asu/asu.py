from pyomoModel import get_full_asu_model
from pyomo.environ import *
import rtolib.util.init_value as ivt

def norminal_simulation():
    model = get_full_asu_model(build_constraint=False, binary_coeff=-0.015)
    ivt.load_init_from_template(model, "GOX17000InitValue.txt")

    solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
    solver = SolverFactory('ipopt', executable=solver_executable)
    results = solver.solve(model, tee=True)
    ivt.to_template(model, "GOX17000SimInitValue_model.txt")

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
    '''
    binary_coeff=-0.01238, iterations = 81
    {'FeedMFlow': 4030.761542015044, 'FeedSplitterOut1Ratio': 0.18402973851281232,
    'FeedSplitterOut2Ratio': 0.5170164329011614, 'HPCCondPrdtMFlow': 1381.045926899585,
    'LPCExt46MFlow': 1750.8448603091485, 'LPCExt15MFlow': 647.4296890235756,
    'ASCCondRefRatio': 27.966333753879805, 'OxSplitterOutARatio': 0.08354167049712301,
    'OBJ': 7378.754395715636}

    binary_coeff=-0.015, iterations = 302
    {'FeedMFlow': 4086.0037360855003, 'FeedSplitterOut1Ratio': 0.18154168241091428,
    'FeedSplitterOut2Ratio': 0.5235463025900136, 'HPCCondPrdtMFlow': 1402.9400273571484,
    'LPCExt46MFlow': 1793.0412318958624, 'LPCExt15MFlow': 633.524517727619,
    'ASCCondRefRatio': 28.14249703353321, 'OxSplitterOutARatio': 0.09843645583070242,
    'OBJ': 7231.011764329876}

    :return:
    '''
    model = get_full_asu_model(build_constraint=True, binary_coeff=-0.01238)
    # model = get_full_asu_model(build_constraint=True, binary_coeff=-0.015)
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
    # ivt.to_template(model, "GOX17000SimInitValue_model.txt")
    ivt.to_template(model, "GOX17000SimInitValue_plant.txt")



if __name__ == "__main__":
    # norminal_simulation()
    nominal_optimization()