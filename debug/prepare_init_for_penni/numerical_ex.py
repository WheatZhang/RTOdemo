from scipy.optimize import minimize
import numpy
import pandas
import matplotlib.pyplot as plt

def hello_world():
    x_k=[(2,),]
    delta=[(2,)]
    i=0
    def updated_model(x, xk):
        return -numpy.exp(-(x[0]-1)*(x[0]-1)) +\
               (x[0]-xk[0])*(-2*xk[0]*numpy.exp(-xk[0]**2)-(-2*(xk[0]-1)*numpy.exp(-(xk[0]-1)**2)))

    fun = lambda x: updated_model(x, x_k[i])
    # fun = lambda x: -numpy.exp(-(x[0]-1)*(x[0]-1))
    cons = []
    for con_i in range(1):
        cons.append({'type': 'ineq', 'fun': (lambda x: delta[i][con_i] ** 2 - (x[con_i] - x_k[i][con_i]) ** 2)})

    res = minimize(fun, x0=x_k[i], method='SLSQP', constraints=cons,
                   options={'ftol': 1e-10, 'disp': True})
    print(res.x[0])


def standard_tr():
    def plant(x):
        return -numpy.exp(-x[0]*x[0])

    def model(x):
        return -numpy.exp(-(x[0]-1)*(x[0]-1))

    def updated_model(x, xk):
        return -numpy.exp(-(x[0]-1)*(x[0]-1)) +\
               (x[0]-xk[0])*(2*xk[0]*numpy.exp(-xk[0]**2)-(2*(xk[0]-1)*numpy.exp(-(xk[0]-1)**2)))

    x_k=[(2,),]
    plant_f=[plant(x_k[0])]
    model_f=[updated_model(x_k[0],x_k[0])]
    accepted_x_k=[x_k[0]]
    accepted_model_f=[model_f[0]]
    accepted_plant_f=[plant_f[0]]
    delta=[(1,)]
    eta_history=[0,]
    data=pandas.DataFrame(index=range(40), columns=("x_k",'accepted_x_k',
                                                    'model_f','accepted_model_f','ref_model_f',
                                                    'plant_f','accepted_plant_f',
                                                    'eta','delta'))
    for i in range(40):
        # print(i)
        fun = lambda x: updated_model(x, accepted_x_k[i])
        cons=[]
        for con_i in range(len(x_k[0])):
            cons.append({'type': 'ineq', 'fun': (lambda x: delta[i][con_i]**2-(x[con_i] - accepted_x_k[i][con_i])**2)})

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
            delta.append(list(numpy.array(delta[i])*2))
            accepted_x_k.append(x_k[i+1])
            accepted_model_f.append(model_f[i+1])
            accepted_plant_f.append(plant_f[i+1])
        elif eta < 0.1:
            delta.append(list(numpy.array(delta[i])/2))
            accepted_x_k.append(accepted_x_k[i])
            accepted_model_f.append(accepted_model_f[i])
            accepted_plant_f.append(accepted_plant_f[i])
        else:
            delta.append(delta[i])
            accepted_x_k.append(x_k[i + 1])
            accepted_model_f.append(model_f[i + 1])
            accepted_plant_f.append(plant_f[i + 1])
        data.loc[i, 'x_k'] = x_k[i][0]
        data.loc[i, 'accepted_x_k'] = accepted_x_k[i][0]
        data.loc[i, 'model_f'] = model_f[i]
        data.loc[i, 'accepted_model_f'] = accepted_model_f[i]
        data.loc[i, 'plant_f'] = plant_f[i]
        data.loc[i, 'accepted_plant_f'] = accepted_plant_f[i]
        data.loc[i, 'eta'] = eta_history[i]
        data.loc[i, 'delta'] = delta[i][0]

    data.to_csv("data/test_tr_standard.txt", sep='\t')


def gradient_tr():
    def plant(x):
        return -numpy.exp(-x[0]*x[0])

    def model(x):
        return -numpy.exp(-(x[0]-1)*(x[0]-1))

    def updated_model(x, xk):
        return -numpy.exp(-(x[0]-1)*(x[0]-1)) +\
               (x[0]-xk[0])*(2*xk[0]*numpy.exp(-xk[0]**2)-(2*(xk[0]-1)*numpy.exp(-(xk[0]-1)**2)))

    x_k=[(2,),]
    plant_f=[plant(x_k[0])]
    model_f=[updated_model(x_k[0],x_k[0])]
    delta=[(1,)]
    accepted_x_k = [x_k[0]]
    accepted_model_f = [model_f[0]]
    accepted_plant_f = [plant_f[0]]
    eta_history=[0,]
    data = pandas.DataFrame(index=range(40), columns=("x_k", 'accepted_x_k',
                                                      'model_f', 'accepted_model_f','ref_model_f',
                                                      'plant_f', 'accepted_plant_f',
                                                      'eta', 'delta'))
    for i in range(40):
        # print(i)
        fun = lambda x: updated_model(x, accepted_x_k[i])
        cons=[]
        # print(delta[i][0]*numpy.exp(-x_k[i][0]*x_k[i][0]))
        cons.append({'type': 'ineq', 'fun': (lambda x: (delta[i][0]*numpy.exp(-accepted_x_k[i][0]*accepted_x_k[i][0]))**2-(x[0] - accepted_x_k[i][0])**2 )})

        res = minimize(fun, x0=accepted_x_k[i], method='SLSQP', constraints=cons,
                       options={'ftol': 1e-10, 'disp': False})
        x_k.append(list(res.x))
        plant_f.append(plant(x_k[i+1]))
        model_f.append(updated_model(x_k[i+1],accepted_x_k[i]))
        ref_model_f = updated_model(accepted_x_k[i], accepted_x_k[i])
        data.loc[i + 1, 'ref_model_f'] = ref_model_f
        eta = (plant_f[i + 1] - accepted_plant_f[i]) / (model_f[i + 1] - ref_model_f)
        # print(eta)
        eta_history.append(eta)
        if eta > 0.9 and eta < 2:
            delta.append(list(numpy.array(delta[i])*2))
            accepted_x_k.append(x_k[i+1])
            accepted_model_f.append(model_f[i+1])
            accepted_plant_f.append(plant_f[i+1])
        elif eta < 0.1:
            delta.append(list(numpy.array(delta[i])/2))
            accepted_x_k.append(accepted_x_k[i])
            accepted_model_f.append(accepted_model_f[i])
            accepted_plant_f.append(accepted_plant_f[i])
        else:
            delta.append(delta[i])
            accepted_x_k.append(x_k[i + 1])
            accepted_model_f.append(model_f[i + 1])
            accepted_plant_f.append(plant_f[i + 1])

        data.loc[i, 'x_k'] = x_k[i][0]
        data.loc[i, 'accepted_x_k'] = accepted_x_k[i][0]
        data.loc[i, 'model_f'] = model_f[i]
        data.loc[i, 'accepted_model_f'] = accepted_model_f[i]
        data.loc[i, 'plant_f'] = plant_f[i]
        data.loc[i, 'accepted_plant_f'] = accepted_plant_f[i]
        data.loc[i, 'eta'] = eta_history[i]
        data.loc[i, 'delta'] = delta[i][0]


    data.to_csv("data/test_tr_gradient.txt", sep='\t')

def plot_comparison():
    standard_tr = pandas.read_csv("data/test_tr_standard.txt", sep='\t', header=0, index_col=0)
    gradient_tr = pandas.read_csv("data/test_tr_gradient.txt", sep='\t', header=0, index_col=0)
    x=range(40)
    plt.plot(x,standard_tr.loc[x,'plant_f'], label='standard_tr')
    plt.plot(x,gradient_tr.loc[x,'plant_f'], label='gradient_tr')
    plt.legend()
    plt.show()


def opposite_curvature_example():
    def plant(x):
        return -x * x

    def model(x, x_k):
        return x * x + (-4 * x_k) * x + 2 * x_k * x_k

    x_k = [10]
    plant_f = [plant(x_k[0])]
    model_f = [model(x_k[0], x_k[0])]
    delta = [3]
    for i in range(20):
        if x_k[i] < 0:
            x_k.append(x_k[i] + delta[i])
        else:
            x_k.append(x_k[i] - delta[i])
        plant_f.append(plant(x_k[i + 1]))
        model_f.append(model(x_k[i], x_k[i + 1]))
        eta = (plant_f[i + 1] - plant_f[i]) / (model_f[i + 1] - model_f[i])
        if eta > 0.9:
            delta.append(delta[i] * 2)
        elif eta < 0.1:
            delta.append(delta[i] / 2)
            x_k[i + 1] = x_k[i]
        else:
            delta.append(delta[i])

    print(x_k)
    print(delta)

    kappa_grad = 1
    x_k = [10]
    plant_f = [plant(x_k[0])]
    model_f = [model(x_k[0], x_k[0])]
    delta = [1]
    for i in range(20):
        if x_k[i] < 0:
            x_k.append(x_k[i] + delta[i] * kappa_grad * numpy.abs(2 * x_k[i]))
        else:
            x_k.append(x_k[i] - delta[i] * kappa_grad * numpy.abs(2 * x_k[i]))
        plant_f.append(plant(x_k[i + 1]))
        model_f.append(model(x_k[i], x_k[i + 1]))
        eta = (plant_f[i + 1] - plant_f[i]) / (model_f[i + 1] - model_f[i])
        if eta > 0.9:
            delta.append(delta[i] * 2)
        elif eta < 0.1:
            delta.append(delta[i] / 2)
            x_k[i + 1] = x_k[i]
        else:
            delta.append(delta[i])

    print(x_k)
    print(delta)

if __name__=="__main__":
    # hello_world()
    standard_tr()
    gradient_tr()
    plot_comparison()