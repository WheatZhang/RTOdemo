import pandas
import matplotlib.pyplot as plt
from rtolib.core.solve import PyomoSimulator
from rtolib.model.quadratic_con import quadratic_con_problem_plant, quadratic_con_problem_model,\
    quadratic_con_problem_description, quadr_obj_con_wrong_curv_model,\
    quadr_con_wrong_curv_model,quadr_obj_wrong_curv_model,quadr_obj_wrong_curv_model2
from pyomo.environ import SolverFactory
import numpy
import os

def generate_contour_data(x,y,func,filename):
    x_len = x.shape[0]
    y_len = y.shape[0]
    with open(filename, 'w') as fp:
        for i in range(x_len):
            if i != 0:
                fp.write('\n')
            for j in range(y_len):
                if j != 0:
                    fp.write('\t')
                fp.write("%.6e" % func(x[i], y[j]))
        fp.write('\n')


def plot_contour(data_file,title):
    global_font_size = 20
    contour_linewidth = 1.5
    global_contour_label_size = 18
    font_factor = numpy.sqrt(1 / 0.39370)
    global_font_size = global_font_size / font_factor
    contour_linewidth = contour_linewidth / font_factor
    global_contour_label_size = global_contour_label_size / font_factor

    font_title = {'family': 'Times New Roman',
             'size': global_font_size,
             }

    data_stored = pandas.read_csv(data_file, index_col=None, header=None, sep='\t')
    x_len = data_stored.shape[0]
    y_len = data_stored.shape[1]
    Z = numpy.zeros(shape=(x_len, y_len))

    for i in range(x_len):
        for j in range(y_len):
            Z[j, i] = data_stored.iloc[i,j]

    a=numpy.linspace(start=-2,stop=4,num=x_len,endpoint=True)
    b=numpy.linspace(start=-4, stop=2, num=y_len,endpoint=True)
    # N = numpy.array([20,60,100,140,160,180,190,200,210])  # 用来指明等高线对应的值为多少时才出图线
    CS = plt.contour(a,b,Z, linewidths=contour_linewidth, cmap=plt.get_cmap('jet'))  # 画出等高线图，cmap表示颜色的图层。
    plt.tick_params(top= ['on', 'true'], right= ['on', 'true'], which='both')
    plt.xlabel('u1')
    plt.ylabel('u2')
    plt.clabel(CS, inline=True, fmt='%d', fontsize=global_contour_label_size)  # 在等高线图里面加入每条线对应的值
    plt.title(title,fontdict=font_title)

plant_simulator = PyomoSimulator(quadratic_con_problem_plant())
norminal_model_simulator = PyomoSimulator(quadratic_con_problem_model())
obj_wrong_model_simulator = PyomoSimulator(quadr_obj_wrong_curv_model())
obj_wrong_model2_simulator = PyomoSimulator(quadr_obj_wrong_curv_model2())
con_wrong_model_simulator = PyomoSimulator(quadr_con_wrong_curv_model())
obj_con_wrong_model_simulator = PyomoSimulator(quadr_obj_con_wrong_curv_model())

plant_simulator.build(quadratic_con_problem_description)
norminal_model_simulator.build(quadratic_con_problem_description)
obj_wrong_model_simulator.build(quadratic_con_problem_description)
obj_wrong_model2_simulator.build(quadratic_con_problem_description)
con_wrong_model_simulator.build(quadratic_con_problem_description)
obj_con_wrong_model_simulator.build(quadratic_con_problem_description)

default_options = {'max_iter': 100}
solver_executable = r"F:\Research\RTOdemo\external\bin\ipopt.exe"
solver = SolverFactory('ipopt', executable=solver_executable)
plant_simulator.set_solver(solver, tee=False, default_options=default_options)
norminal_model_simulator.set_solver(solver, tee=False, default_options=default_options)
obj_wrong_model_simulator.set_solver(solver, tee=False, default_options=default_options)
obj_wrong_model2_simulator.set_solver(solver, tee=False, default_options=default_options)
con_wrong_model_simulator.set_solver(solver, tee=False, default_options=default_options)
obj_con_wrong_model_simulator.set_solver(solver, tee=False, default_options=default_options)

def plant_simulation_callback(u1, u2, var):
    outputs, solve_status = plant_simulator.simulate({'u1': u1, 'u2': u2}, param_values=None,
                                                     use_homo=False)
    return outputs[var]

def model_simulation_callback(model_name, u1, u2, param_values, var):
    if model_name == "norminal":
        outputs, solve_status = norminal_model_simulator.simulate({'u1': u1, 'u2': u2}, param_values=param_values,
                                                         use_homo=False)
    elif model_name == "obj_wrong_curv":
        outputs, solve_status = obj_wrong_model_simulator.simulate({'u1': u1, 'u2': u2}, param_values=param_values,
                                                                  use_homo=False)
    elif model_name == "obj_partial_wrong_curv":
        outputs, solve_status = obj_wrong_model2_simulator.simulate({'u1': u1, 'u2': u2}, param_values=param_values,
                                                                  use_homo=False)
    elif model_name == "con_wrong_curv":
        outputs, solve_status = con_wrong_model_simulator.simulate({'u1': u1, 'u2': u2}, param_values=param_values,
                                                                  use_homo=False)
    elif model_name == "obj_con_wrong_curv":
        outputs, solve_status = obj_con_wrong_model_simulator.simulate({'u1': u1, 'u2': u2}, param_values=param_values,
                                                                  use_homo=False)
    else:
        raise ValueError("Wrong model name.")
    return outputs[var]


def generate_all_contour_data(resolution):
    datafile_folder = "data/contour_data/res%d" % resolution
    if not os.path.exists("data/contour_data"):
        os.makedirs("data/contour_data")
    def get_data_path(name):
        return datafile_folder + "_" + name+ ".txt"

    u1_points = numpy.linspace(start=-2, stop=4, num=resolution, endpoint=True)
    u2_points = numpy.linspace(start=-4, stop=2, num=resolution, endpoint=True)

    # ---------------------------
    func = lambda u1, u2: plant_simulation_callback(u1, u2, 'cost')
    generate_contour_data(u1_points, u2_points, func, get_data_path('cost_plant'))

    func = lambda u1, u2: plant_simulation_callback(u1, u2, 'con')
    generate_contour_data(u1_points, u2_points, func, get_data_path('con_plant'))

    model_names = ["norminal", "obj_wrong_curv", "obj_partial_wrong_curv", "con_wrong_curv",\
            "obj_con_wrong_curv"]
    for model_name in model_names:
        func = lambda u1, u2: model_simulation_callback(model_name, u1, u2, {}, 'cost')
        generate_contour_data(u1_points, u2_points, func, get_data_path('cost_model_'+model_name))

        func = lambda u1, u2: model_simulation_callback(model_name, u1, u2, {}, 'con')
        generate_contour_data(u1_points, u2_points, func, get_data_path('con_model_'+model_name))