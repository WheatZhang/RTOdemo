from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy
import pandas
from scipy import ndimage

def plant_simulation(F, S0):
    def ode_fun(t, y):
        # y=[X,P,S,V]
        mu_x=0.092
        mu_p=0.005
        K_x=0.15
        K_p=0.0002
        K_I=0.1
        K_H=0.04
        Y_xs=0.45
        Y_ps=0.9
        m_x=0.014
        S_F=600
        Beta=0.00062

        dydt=[0,0,0,0]

        dydt[0] = mu_x * y[2] * y[0] / (K_x * y[0] + y[2]) - y[0] * \
                   (F - Beta * y[3]) / y[3]
        dydt[1] = mu_p * y[2] * y[0] / (K_p + y[2] + y[2] ** 2 / K_I) - y[1] * \
                   (F - Beta * y[3]) / y[3] - K_H * y[1]
        dydt[2] = -mu_x * y[2] * y[0] / (Y_xs * (K_x * y[0] + y[2])) - mu_p * \
                   y[2] * y[0] / (Y_ps * (K_p + y[2] + y[2] ** 2 / K_I)) - m_x * y[0] + S_F * \
                   F / y[3] - y[2] * (F - Beta * y[3]) / y[3]
        dydt[3] = F - Beta * y[3]
        return dydt

    sol = solve_ivp(ode_fun, [0, 192], [0.1,0,S0,100],dense_output=True,rtol=1e-6)

    # t = numpy.linspace(0, 192, 100)
    # z = sol.sol(t)
    # plt.plot(t, z.T)
    # plt.xlabel('t')
    # plt.legend(['X', 'P', 'S', 'V'], shadow=True)
    # plt.title('penicillin')
    # plt.show()

    y_tf=[0,0,0,0]
    for i in range(4):
        y_tf[i]=sol.y[i][-1]
    return y_tf

def model_simulation(F, S0):
    def ode_fun(t, y):
        # y=[X,P,S,V]
        mu_x=0.092
        mu_p=0.005
        K_x=0.15
        K_p=0.0002
        K_I=0.1
        K_H=0.04
        Y_xs=0.45
        Y_ps=0.9
        m_x=0.014
        S_F=600
        Beta=0.00062

        dydt=[0,0,0,0]

        dydt[0] = mu_x * y[2] * y[0] / (K_x * y[0] + y[2]) - y[0] * \
                   (F - Beta * y[3]) / y[3]
        dydt[1] = mu_p * y[2] * y[0] / (K_p + y[2] + y[2] ** 2 / K_I) - y[1] * \
                   (F - Beta * y[3]) / y[3]
        dydt[2] = -mu_x * y[2] * y[0] / (Y_xs * (K_x * y[0] + y[2])) - mu_p * \
                   y[2] * y[0] / (Y_ps * (K_p + y[2] + y[2] ** 2 / K_I)) - m_x * y[0] + S_F * \
                   F / y[3] - y[2] * (F - Beta * y[3]) / y[3]
        dydt[3] = F - Beta * y[3]
        return dydt

    sol = solve_ivp(ode_fun, [0, 192], [0.1,0,S0,100],dense_output=True,rtol=1e-6)

    # t = numpy.linspace(0, 192, 100)
    # z = sol.sol(t)
    # plt.plot(t, z.T)
    # plt.xlabel('t')
    # plt.legend(['X', 'P', 'S', 'V'], shadow=True)
    # plt.title('penicillin')
    # plt.show()

    y_tf=[0,0,0,0]
    for i in range(4):
        y_tf[i]=sol.y[i][-1]
    return y_tf

def get_simulation_data_plant():
    n=0
    data = pandas.DataFrame(columns=['F','S0','obj','con'])
    for F in numpy.linspace(start=0.1, stop=0.3, num=20):
        for S0 in numpy.linspace(start=10, stop=60, num=20):
            print(n)
            y_tf = plant_simulation(F=F, S0=S0)
            data = data.append({'F':F,'S0':S0,'obj':y_tf[1]*120,'con':y_tf[3]}, ignore_index=True)
            n+=1
    data.to_csv("data/plant_data_20.txt", sep='\t')

def get_simulation_data_test():
    n = 0
    data = pandas.DataFrame(columns=['F', 'S0', 'obj', 'con'])
    for F in numpy.linspace(start=0.1, stop=0.3, num=20):
        for S0 in numpy.linspace(start=10, stop=60, num=20):
            data = data.append({'F': F, 'S0': S0, 'obj': S0*S0*S0, 'con': F*F*F*1e5}, ignore_index=True)
            n += 1
    data.to_csv("data/test_data_20.txt", sep='\t')

def get_simulation_data_model():
    n=0
    data = pandas.DataFrame(columns=['F','S0','obj','con'])
    for F in numpy.linspace(start=0.1, stop=0.3, num=20):
        for S0 in numpy.linspace(start=10, stop=60, num=20):
            print(n)
            y_tf = model_simulation(F=F, S0=S0)
            data = data.append({'F':F,'S0':S0,'obj':y_tf[1]*120,'con':y_tf[3]}, ignore_index=True)
            n+=1
    data.to_csv("data/model_data_20.txt", sep='\t')

def plot_contour(data_file, x_list, y_list, filename):
    data_stored = pandas.read_csv(data_file, index_col=0, header=0, sep='\t')
    x_len = len(x_list)
    y_len = len(y_list)
    Z = numpy.zeros(shape=(x_len, y_len))

    for i in range(x_len):
        for j in range(y_len):
            Z[j, i] = data_stored.iloc[i*y_len+j].loc['obj']

    CS = plt.contour(x_list,y_list,Z,20)
    plt.clabel(CS, inline=True, fmt='%d')
    plt.savefig(filename)
    plt.close()

def plot_con_contour(data_file, x_list, y_list, filename):
    data_stored = pandas.read_csv(data_file, index_col=0, header=0, sep='\t')
    x_len = len(x_list)
    y_len = len(y_list)
    Z = numpy.zeros(shape=(x_len, y_len))

    for i in range(x_len):
        for j in range(y_len):
            Z[j, i] = data_stored.iloc[i*y_len+j].loc['con']

    CS = plt.contour(x_list,y_list,Z,20)
    plt.clabel(CS, inline=True, fmt='%d')
    plt.savefig(filename)
    plt.close()

def plot_hessian_minimal_eigen(data_file, x_list, y_list, filename):
    data_stored = pandas.read_csv(data_file, index_col=0, header=0, sep='\t')
    mesh_size=len(x_list)
    Z = numpy.zeros(shape=(mesh_size, mesh_size))
    for i in range(mesh_size):
        for j in range(mesh_size):
            Z[i, j] = data_stored.iloc[i * mesh_size + j].loc['obj']
    minimal_eigen_value=numpy.zeros(shape=(mesh_size, mesh_size))
    for i in range(mesh_size):
        if i == 0:
            i_index = [0, 1, 2]
        elif i == mesh_size-1:
            i_index = [mesh_size-3, mesh_size-2, mesh_size-1]
        else:
            i_index = [i - 1, i, i + 1]
        for j in range(mesh_size):
            if j == 0:
                j_index = [0, 1, 2]
            elif j == mesh_size - 1:
                j_index = [mesh_size - 3, mesh_size - 2, mesh_size - 1]
            else:
                j_index = [j-1, j, j+1]
            local_data=numpy.zeros(shape=(3, 3))
            for k1 in range(3):
                for k2 in range(3):
                    local_data[k1,k2]=Z[i_index[k1],j_index[k2]]
            hessian=numpy.zeros(shape=(2, 2))
            hessian[0,0]=(local_data[2,0]-local_data[1,0])-(local_data[1,0]-local_data[0,0])
            hessian[1,1]=(local_data[0,2]-local_data[0,1])-(local_data[0,1]-local_data[0,0])
            hessian[0,1]=(local_data[1,1]-local_data[1,0])-(local_data[0,1]-local_data[0,0])
            hessian[1,0]=(local_data[1,1]-local_data[1,0])-(local_data[0,1]-local_data[0,0])

            minimal_eigen_value[j, i] = max(numpy.linalg.eig(hessian)[0])

    minimal_eigen_value=ndimage.filters.gaussian_filter(minimal_eigen_value, [1,1], mode='constant')
    CS = plt.contour(x_list, y_list, minimal_eigen_value,10)
    plt.clabel(CS, inline=True, fmt='%d')
    plt.savefig(filename)
    plt.close()

def plot_hessian_along_y_minimal_eigen(data_file, x_list, y_list, filename):
    data_stored = pandas.read_csv(data_file, index_col=0, header=0, sep='\t')
    mesh_size=len(x_list)
    Z = numpy.zeros(shape=(mesh_size, mesh_size))
    for i in range(mesh_size):
        for j in range(mesh_size):
            Z[i, j] = data_stored.iloc[i * mesh_size + j].loc['obj']
    # Z = ndimage.filters.gaussian_filter(Z, [2, 2], mode='constant')
    minimal_eigen_value=numpy.zeros(shape=(mesh_size, mesh_size))
    for i in range(mesh_size):
        if i == 0:
            i_index = [0, 1, 2]
        elif i == mesh_size-1:
            i_index = [mesh_size-3, mesh_size-2, mesh_size-1]
        else:
            i_index = [i - 1, i, i + 1]
        for j in range(mesh_size):
            if j == 0:
                j_index = [0, 1, 2]
            elif j == mesh_size - 1:
                j_index = [mesh_size - 3, mesh_size - 2, mesh_size - 1]
            else:
                j_index = [j-1, j, j+1]
            local_data=numpy.zeros(shape=(3, 3))
            for k1 in range(3):
                for k2 in range(3):
                    local_data[k1,k2]=Z[i_index[k1],j_index[k2]]
            # hessian=(local_data[2,0]-local_data[1,0])-(local_data[1,0]-local_data[0,0])
            hessian=(local_data[0,2]-local_data[0,1])-(local_data[0,1]-local_data[0,0])

            minimal_eigen_value[j, i] = hessian

    minimal_eigen_value = ndimage.filters.gaussian_filter(minimal_eigen_value, [1, 1], mode='constant')
    CS = plt.contour(x_list, y_list, minimal_eigen_value,10)
    plt.clabel(CS, inline=True, fmt='%d')
    plt.savefig(filename)
    plt.close()


if __name__ == "__main__":
    get_simulation_data_plant()
    F = numpy.linspace(start=0.1, stop=0.3, num=20)
    S0 = numpy.linspace(start=10, stop=60, num=20)
    plot_contour("data/plant_data_20.txt", F, S0, "pic/plant_contour20.png")
    plot_con_contour("data/plant_data_20.txt", F, S0, "pic/plant_con_contour20.png")
    plot_hessian_minimal_eigen("data/plant_data_20.txt", F, S0, "pic/plant_hessian_min_eig20.png")
    plot_hessian_along_y_minimal_eigen("data/plant_data_20.txt", F, S0, "pic/plant_hes_alongy_min_eig20.png")

    get_simulation_data_model()
    # F = numpy.linspace(start=0.1, stop=0.3, num=20)
    # S0 = numpy.linspace(start=10, stop=60, num=20)
    plot_contour("data/model_data_20.txt", F, S0, "pic/model_contour20.png")
    plot_con_contour("data/model_data_20.txt", F, S0, "pic/model_con_contour20.png")
    plot_hessian_minimal_eigen("data/model_data_20.txt", F, S0, "pic/model_hessian_min_eig20.png")
    plot_hessian_along_y_minimal_eigen("data/model_data_20.txt", F, S0, "pic/model_hes_alongy_min_eig20.png")

    # get_simulation_data_test()
    # F = numpy.linspace(start=0.1, stop=0.3, num=20)
    # S0 = numpy.linspace(start=10, stop=60, num=20)
    # plot_contour("data/test_data_20.txt", F, S0, "pic/test_contour20.png")
    # plot_con_contour("data/test_data_20.txt", F, S0, "pic/test_con_contour20.png")
    # plot_hessian_minimal_eigen("data/test_data_20.txt", F, S0, "pic/test_hessian_min_eig20.png")
    # plot_hessian_along_y_minimal_eigen("data/test_data_20.txt", F, S0, "pic/test_hes_alongy_min_eig20.png")

