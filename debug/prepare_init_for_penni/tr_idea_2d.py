from scipy.optimize import minimize
import numpy
import pandas
import matplotlib.pyplot as plt

def standard_tr():
    def plant(x):
        return -numpy.exp(-x[0]*x[0])+2*x[1]*x[1]

    def model(x):
        return -numpy.exp(-(x[0]-1)*(x[0]-1))+x[1]*x[1]

    def updated_model(x, xk):
        return -numpy.exp(-(x[0]-1)*(x[0]-1))+x[1]*x[1] +\
               (x[0]-xk[0])*(2*xk[0]*numpy.exp(-xk[0]**2)-(2*(xk[0]-1)*numpy.exp(-(xk[0]-1)**2)))+\
               (x[1]-xk[1])*(4*xk[1]-2*xk[1])

    x_k=[(1,1),]
    plant_f=[plant(x_k[0])]
    model_f=[updated_model(x_k[0],x_k[0])]
    accepted_x_k=[x_k[0]]
    accepted_model_f=[model_f[0]]
    accepted_plant_f=[plant_f[0]]
    delta=[0.2]
    eta_history=[0,]
    data=pandas.DataFrame(index=range(40), columns=("x_k1",'accepted_x_k1',
                                                    "x_k2", 'accepted_x_k2',
                                                    'model_f','accepted_model_f','ref_model_f',
                                                    'plant_f','accepted_plant_f',
                                                    'eta','delta'))
    for i in range(40):
        # print(i)
        fun = lambda x: updated_model(x, accepted_x_k[i])
        cons=[]
        cons.append({'type': 'ineq', 'fun': (lambda x: delta[i]**2-(x[0] - accepted_x_k[i][0])**2-(x[1] - accepted_x_k[i][1])**2)})

        res = minimize(fun, x0=accepted_x_k[i], method='SLSQP', constraints=cons,
                       options={'ftol': 1e-10, 'disp': False})
        x_k.append(list(res.x))
        plant_f.append(plant(x_k[i+1]))
        model_f.append(updated_model(x_k[i+1],accepted_x_k[i]))
        ref_model_f=updated_model(accepted_x_k[i],accepted_x_k[i])
        data.loc[i+1, 'ref_model_f'] = ref_model_f
        eta=(plant_f[i+1]-accepted_plant_f[i])/(model_f[i+1]-ref_model_f)
        eta_history.append(eta)
        # print(eta)
        if eta > 0.9:
            delta.append(delta[i]*2)
            accepted_x_k.append(x_k[i+1])
            accepted_model_f.append(model_f[i+1])
            accepted_plant_f.append(plant_f[i+1])
        elif eta < 0.1:
            delta.append(delta[i]/2)
            accepted_x_k.append(accepted_x_k[i])
            accepted_model_f.append(accepted_model_f[i])
            accepted_plant_f.append(accepted_plant_f[i])
        else:
            delta.append(delta[i])
            accepted_x_k.append(x_k[i + 1])
            accepted_model_f.append(model_f[i + 1])
            accepted_plant_f.append(plant_f[i + 1])
        data.loc[i, 'x_k1'] = x_k[i][0]
        data.loc[i, 'accepted_x_k1'] = accepted_x_k[i][0]
        data.loc[i, 'x_k2'] = x_k[i][1]
        data.loc[i, 'accepted_x_k2'] = accepted_x_k[i][1]
        data.loc[i, 'model_f'] = model_f[i]
        data.loc[i, 'accepted_model_f'] = accepted_model_f[i]
        data.loc[i, 'plant_f'] = plant_f[i]
        data.loc[i, 'accepted_plant_f'] = accepted_plant_f[i]
        data.loc[i, 'eta'] = eta_history[i]
        data.loc[i, 'delta'] = delta[i]

    data.to_csv("data/test_2d_tr_standard.txt", sep='\t')

def sector_tr():
    def plant(x):
        return -numpy.exp(-x[0]*x[0])+2*x[1]*x[1]

    def get_plant_gradient(xk):
        return [2*xk[0]*numpy.exp(-xk[0]**2), 4*xk[1]]

    def model(x):
        return -numpy.exp(-(x[0]-1)*(x[0]-1))+x[1]*x[1]

    def updated_model(x, xk):
        return -numpy.exp(-(x[0]-1)*(x[0]-1))+x[1]*x[1] +\
               (x[0]-xk[0])*(2*xk[0]*numpy.exp(-xk[0]**2)-(2*(xk[0]-1)*numpy.exp(-(xk[0]-1)**2)))+\
               (x[1]-xk[1])*(4*xk[1]-2*xk[1])

    x_k=[(1,1),]
    plant_f=[plant(x_k[0])]
    model_f=[updated_model(x_k[0],x_k[0])]
    accepted_x_k=[x_k[0]]
    accepted_model_f=[model_f[0]]
    accepted_plant_f=[plant_f[0]]
    delta=[0.2]
    delta_gamma_factor=10
    eta_history=[0,]
    data=pandas.DataFrame(index=range(40), columns=("x_k1",'accepted_x_k1',
                                                    "x_k2", 'accepted_x_k2',
                                                    'model_f','accepted_model_f','ref_model_f',
                                                    'plant_f','accepted_plant_f',
                                                    'eta','delta','gamma'))
    for i in range(40):
        # print(i)
        fun = lambda x: updated_model(x, accepted_x_k[i])
        plant_gradient=get_plant_gradient(accepted_x_k[i])
        cons=[]
        cons.append({'type': 'ineq', 'fun': (lambda x: delta[i]**2-(x[0] - accepted_x_k[i][0])**2-(x[1] - accepted_x_k[i][1])**2)})
        cons.append({'type': 'ineq', 'fun': (lambda x: (x[0] - accepted_x_k[i][0])*plant_gradient[0]+\
                                                       (x[1] - accepted_x_k[i][1])*plant_gradient[1]-\
                                                       (1-max(delta_gamma_factor*delta[i],0.1))*numpy.sqrt(((x[0] - accepted_x_k[i][0])**2+(x[1] - accepted_x_k[i][1])**2)*\
                                                       (plant_gradient[0]**2+plant_gradient[1]**2)))})

        res = minimize(fun, x0=accepted_x_k[i], method='SLSQP', constraints=cons,
                       options={'ftol': 1e-10, 'disp': False})
        x_k.append(list(res.x))
        plant_f.append(plant(x_k[i+1]))
        model_f.append(updated_model(x_k[i+1],accepted_x_k[i]))
        ref_model_f=updated_model(accepted_x_k[i],accepted_x_k[i])
        data.loc[i+1, 'ref_model_f'] = ref_model_f
        eta=(plant_f[i+1]-accepted_plant_f[i])/(model_f[i+1]-ref_model_f)
        # if model_f[i+1]==ref_model_f:
        #     eta=0
        eta_history.append(eta)
        # print(eta)
        if eta > 0.9:
            delta.append(delta[i]*2)
            accepted_x_k.append(x_k[i+1])
            accepted_model_f.append(model_f[i+1])
            accepted_plant_f.append(plant_f[i+1])
        elif eta < 0.1:
            delta.append(delta[i]/2)
            accepted_x_k.append(accepted_x_k[i])
            accepted_model_f.append(accepted_model_f[i])
            accepted_plant_f.append(accepted_plant_f[i])
        else:
            delta.append(delta[i])
            accepted_x_k.append(x_k[i + 1])
            accepted_model_f.append(model_f[i + 1])
            accepted_plant_f.append(plant_f[i + 1])
        data.loc[i, 'x_k1'] = x_k[i][0]
        data.loc[i, 'accepted_x_k1'] = accepted_x_k[i][0]
        data.loc[i, 'x_k2'] = x_k[i][1]
        data.loc[i, 'accepted_x_k2'] = accepted_x_k[i][1]
        data.loc[i, 'model_f'] = model_f[i]
        data.loc[i, 'accepted_model_f'] = accepted_model_f[i]
        data.loc[i, 'plant_f'] = plant_f[i]
        data.loc[i, 'accepted_plant_f'] = accepted_plant_f[i]
        data.loc[i, 'eta'] = eta_history[i]
        data.loc[i, 'delta'] = delta[i]
        data.loc[i, 'gamma'] = (1-max(delta_gamma_factor*delta[i],0.1))

    data.to_csv("data/test_2d_tr_sector.txt", sep='\t')


def plot_comparison():
    standard_tr = pandas.read_csv("data/test_2d_tr_standard.txt", sep='\t', header=0, index_col=0)
    sector_tr = pandas.read_csv("data/test_2d_tr_sector.txt", sep='\t', header=0, index_col=0)
    x=range(40)

    plt.plot(x,sector_tr.loc[x,'x_k1'], label='sector_tr_x1', marker='o')
    plt.plot(x, sector_tr.loc[x, 'x_k2'], label='sector_tr_x2', marker='o')
    plt.plot(x, standard_tr.loc[x, 'x_k1'], label='standard_tr_x1', linestyle='dashed')
    plt.plot(x,standard_tr.loc[x,'x_k2'], label='standard_tr_x2', linestyle='dashed')

    plt.legend()
    plt.show()


if __name__=="__main__":
    standard_tr()
    sector_tr()
    plot_comparison()