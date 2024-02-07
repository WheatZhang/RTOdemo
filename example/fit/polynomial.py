from rtolib.util.fit_poly_krig import PolynomialFitter, MonomialBaseGroup,\
    GaussianRBFNetworkFunc,GaussianRBFNetworkFuncFitter, KrigingFunction,\
    KrigingFitter
from rtolib.util.cvx_quadr import get_hypercube_sampling
import numpy as np


def ex_polynomial_fitting():
    max_order = 2
    dimension = 3
    way = 'Compound'

    fitter = PolynomialFitter(dimension, MonomialBaseGroup.fitting_base_library(dimension,way,max_order))

    def my_func(x):
        return x[0]**2+x[1]+x[0]*x[2]
    n_data = 100
    lb = []
    ub = []
    for i in range(3):
        lb.append(0)
        ub.append(i+1)
    X = get_hypercube_sampling(dimension, n_data, lb, ub, seed=1)
    Y = np.zeros((n_data,))
    for i in range(n_data):
        Y[i] = my_func(X[i])
    polynomial_func = fitter.fit(X, Y, method='minimize')
    fitted_Y = np.zeros((n_data,))
    for i in range(n_data):
        fitted_Y[i] = polynomial_func.evaluate(X[i])
    print(max(abs(Y-fitted_Y)))
    print(max(Y)-min(Y))
    polynomial_func.set_var_name(['x0','x1','x2'],'y')
    print(polynomial_func.generate_text(name='text_func'))


def ex_kriging_fitting():
    max_order = 2
    dimension = 3
    way = 'Compound'
    base_group = MonomialBaseGroup.fitting_base_library(dimension, way, max_order)
    bandwidth = np.array([0.5,0.5,0.5])

    fitter = KrigingFitter(dimension, base_group, bandwidth)

    def my_func(x):
        return x[0]**2+x[1]+x[0]*x[2]+np.sin(x[1]*x[2]-x[0])
    n_data = 100
    lb = []
    ub = []
    for i in range(3):
        lb.append(0)
        ub.append(i+1)
    X = get_hypercube_sampling(dimension, n_data, lb, ub, seed=1)
    Y = np.zeros((n_data,))
    for i in range(n_data):
        Y[i] = my_func(X[i])
    rbfn_func = fitter.fit(X, Y)
    fitted_Y = np.zeros((n_data,))
    for i in range(n_data):
        fitted_Y[i] = rbfn_func.evaluate(X[i])
    print(max(abs(Y-fitted_Y)))
    print(max(Y)-min(Y))
    rbfn_func.set_var_name(['v0','v1','v2'],'y')
    string = rbfn_func.generate_text(name='text_func')
    with open("ex_kriging.py",'w') as fp:
        fp.write("import numpy as np\n\n")
        fp.write(string)

    from ex_kriging import text_func

    fitted_Y2 = np.zeros((n_data,))
    for i in range(n_data):
        fitted_Y2[i] = text_func(*X[i])
    print(max(abs(Y - fitted_Y2)))
    print(max(Y) - min(Y))



if __name__ == "__main__":
    # ex_polynomial_fitting()
    ex_kriging_fitting()