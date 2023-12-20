import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
import pandas
from sko.PSO import PSO

def soft_threshold(x, width, x0=0, constant=0):
    if x < x0:
        ret = 0
    elif x > width+x0:
        ret = 1
    else:
        ret = (1-np.cos(-np.pi*(x-x0)/width))/2
    return ret+constant

def test_soft_threshold():
    width=0.2
    x = np.linspace(start=-2*width,stop=2*width, num=100)
    y = np.zeros((100,))
    for i in range(100):
        y[i] = soft_threshold(x[i], width)
    fig = plt.figure()
    plt.plot(x, y)
    plt.show()

def get_solution(F, S0, t_step=0.1, label='plant'):
    mu_x = 0.092
    K_x = 0.15
    mu_P = 0.005
    K_P = 0.0002
    K_I = 0.1
    K_H = 0.04
    Y_XS = 0.45
    Y_PS = 0.9
    m_X = 0.014
    s_f = 600

    def ode_penicillin(t, s, F):
        X = s[0]
        P = s[1]
        S = s[2]
        V = s[3]

        expr_1 =  mu_x * S * X / (K_x * X + S)
        expr_2 =  mu_P * S * X / (K_P + S + S * S / K_I)
        dVdt = F - 6.226e-4 * V
        dXdt = expr_1 - X / V * dVdt
        if label == 'plant':
            dPdt = expr_2 - K_H * P - P / V * dVdt
        elif label == 'model':
            dPdt = expr_2 - P / V * dVdt
        else:
            raise ValueError("unknown label")
        dSdt = -expr_1 / Y_XS - expr_2 / Y_PS - m_X * X + F * s_f / V - \
                    S / V * dVdt

        epsilon=1e-4
        # if V<epsilon and dVdt<0:
        #     dVdt=0
        # if X<epsilon and dXdt<0:
        #     dXdt=0
        # if P<epsilon and dPdt<0:
        #     dPdt=0
        # if S<epsilon and dSdt<0:
        #     dSdt=0
        # dVdt=dVdt*soft_threshold(V,width=epsilon)
        # dXdt=dXdt*soft_threshold(X,width=epsilon)
        # dPdt=dPdt*soft_threshold(P,width=epsilon)
        if dSdt<0:
            if label == 'plant':
                dSdt=dSdt*soft_threshold(S,x0=epsilon, width=1e-2)
            elif label == 'model':
                dSdt = dSdt * soft_threshold(S, x0=epsilon, width=1e-1)

        ret = np.zeros(shape=(4,))
        ret[0] = dXdt
        ret[1] = dPdt
        ret[2] = dSdt
        ret[3] = dVdt

        return ret

    ode_func = lambda t,s: ode_penicillin(t,s,F)

    t_eval = np.arange(0, 200, t_step)
    sol = solve_ivp(ode_func, [0, 200], [0.1, 0, S0, 100], t_eval=t_eval)
    return sol

def plot_result(sol):
    plt.figure(figsize = (9,9))
    plt.subplot(221)
    plt.plot(sol.t, sol.y[0])
    plt.xlabel('t')
    plt.ylabel('X(t)')
    plt.subplot(222)
    plt.plot(sol.t, sol.y[1])
    plt.xlabel('t')
    plt.ylabel('P(t)')
    plt.subplot(223)
    plt.plot(sol.t, sol.y[2])
    plt.xlabel('t')
    plt.ylabel('S(t)')
    plt.subplot(224)
    plt.plot(sol.t, sol.y[3])
    plt.xlabel('t')
    plt.ylabel('V(t)')
    plt.tight_layout()
    plt.show()

def generate_init_value_file(sol):
    var_names=['X', 'P', 'S', 'V']
    with open("init_value/init_from_ode.txt", "w") as init_value_file:
        for i in range(len(var_names)):
            init_value_file.write(var_names[i]+"\nt\n")
            t = sol.t
            val = sol.y[i]
            for j in range(len(t)):
                init_value_file.write(str(t[j]))
                init_value_file.write("\t")
                init_value_file.write(str(val[j]))
                init_value_file.write("\n")
            init_value_file.write("\n")

        t = sol.t
        val = sol.y
        dsdt = np.zeros(shape=(len(t), 4))
        for j in range(len(t)):
            s = np.zeros(shape=(4,))
            s[0] = val[0][j]
            s[1] = val[1][j]
            s[2] = val[2][j]
            s[3] = val[3][j]
            ret=ode_penicillin(t[j], s, F)
            dsdt[j,0] = ret[0]
            dsdt[j,1] = ret[1]
            dsdt[j,2] = ret[2]
            dsdt[j,3] = ret[3]

        for i in range(len(var_names)):
            init_value_file.write("d"+var_names[i] + "dt\nt\n")
            for j in range(len(t)):
                init_value_file.write(str(t[j]))
                init_value_file.write("\t")
                init_value_file.write(str(dsdt[j,i]))
                init_value_file.write("\n")
            init_value_file.write("\n")


def get_plant_contour_data(num_point):
    filename = "data/plant_contour_data%d.txt"%num_point
    x_points = np.linspace(0, 0.5, num_point)
    y_points = np.linspace(0, 100, num_point)
    with open(filename, 'w') as fp:
        fp.write("obj\tcon\n")
        for i in range(num_point):
            for j in range(num_point):
                print(i*num_point+j)
                sol = get_solution(x_points[i], y_points[j],t_step=0.001,label='plant')
                fp.write("%.6e" % sol.y[1][-1])
                fp.write('\t')
                fp.write("%.6e" % sol.y[3][-1])
                fp.write('\n')

def draw_plant_characteristic(num_point):
    data_file = "data/plant_contour_data%d.txt"%num_point

    data_stored = pandas.read_csv(data_file, index_col=None, header=0, sep='\t')
    Z_obj = np.zeros(shape=(num_point, num_point))
    Z_con = np.zeros(shape=(num_point, num_point))

    for i in range(num_point):
        for j in range(num_point):
            Z_obj[j, i] = data_stored.iloc[i*num_point+j, 0]
            Z_con[j, i] = data_stored.iloc[i*num_point+j, 1]

    x_points = np.linspace(0, 0.4, num_point)
    y_points = np.linspace(0, 100, num_point)
    fig=plt.figure(figsize=(10,6))
    plt.subplot(121)
    CS = plt.contourf(x_points,y_points,Z_obj)
    plt.title('obj')
    plt.colorbar()
    plt.subplot(122)
    CS = plt.contourf(x_points, y_points, Z_con)
    plt.title('con')
    plt.colorbar()
    plt.savefig("pic/plant_contour.png", dpi=600)
    plt.show()

def get_model_contour_data(num_point):
    filename = "data/model_contour_data%d.txt"%num_point
    x_points = np.linspace(0, 0.5, num_point)
    y_points = np.linspace(0, 100, num_point)
    with open(filename, 'w') as fp:
        fp.write("obj\tcon\n")
        for i in range(num_point):
            for j in range(num_point):
                print(i*num_point+j)
                sol = get_solution(x_points[i], y_points[j],t_step=0.001, label='model')
                fp.write("%.6e" % sol.y[1][-1])
                fp.write('\t')
                fp.write("%.6e" % sol.y[3][-1])
                fp.write('\n')

def draw_model_characteristic(num_point):
    data_file = "data/model_contour_data%d.txt"%num_point

    data_stored = pandas.read_csv(data_file, index_col=None, header=0, sep='\t')
    Z_obj = np.zeros(shape=(num_point, num_point))
    Z_con = np.zeros(shape=(num_point, num_point))

    for i in range(num_point):
        for j in range(num_point):
            Z_obj[j, i] = data_stored.iloc[i*num_point+j, 0]
            Z_con[j, i] = data_stored.iloc[i*num_point+j, 1]

    x_points = np.linspace(0, 0.4, num_point)
    y_points = np.linspace(0, 100, num_point)
    fig=plt.figure(figsize=(10,6))
    plt.subplot(121)
    CS = plt.contourf(x_points,y_points,Z_obj)
    plt.title('obj')
    plt.colorbar()
    plt.subplot(122)
    CS = plt.contourf(x_points, y_points, Z_con)
    plt.title('con')
    plt.colorbar()
    plt.savefig("pic/model_contour.png", dpi=600)
    plt.show()

def try_pso():
    def obj_func(x):
        F, S0 = x
        sol = get_solution(F, S0, t_step=0.1, label='plant')
        merit_f=-sol.y[1][-1]+100*max(sol.y[3][-1]-120, 0)
        return merit_f

    # pop=population w=inertia
    pso = PSO(func=obj_func, n_dim=2, pop=10, max_iter=50, lb=[0, 0], ub=[0.4, 100], \
              w=0.8, c1=0.5, c2=0.5)
    pso.run()
    print('best_x is ', pso.gbest_x, 'best_y is', pso.gbest_y)

    # %% Plot the result
    import matplotlib.pyplot as plt

    plt.plot(pso.gbest_y_hist)
    plt.show()

if __name__ == "__main__":
    # sol = get_solution(0.1728, 54.72, t_step=1)
    # plot_result(sol)
    # test_soft_threshold()
    # num_point = 30
    # get_plant_contour_data(num_point)
    # draw_plant_characteristic(num_point)
    # num_point = 30
    # get_model_contour_data(num_point)
    # draw_model_characteristic(num_point)
    try_pso()