import scipy.stats.qmc
import numpy as np

def is_value_dict_same(dict1, dict2):
    for k in dict1:
        if abs(dict1[k]-dict2[k]) > 1e-6:
            return False
    return True

def distance_of_two_dicts(dict1, dict2):
    norm = 0
    for k in dict1.keys():
        norm+=(dict1[k]-dict2[k])**2
    norm = np.sqrt(norm)
    return norm


def save_iteration_data_in_dict(history_data, filename):
    headers = history_data[1].keys()
    rows = history_data.keys()
    with open(filename, 'w') as fp:
        fp.write("RTO_iter")
        for h in headers:
            fp.write('\t' + h)
        fp.write('\n')
        for r in rows:
            fp.write("%d" % r)
            for h in headers:
                try:
                    if isinstance(history_data[r][h], str):
                        fp.write('\t%s' % history_data[r][h])
                    else:
                        fp.write('\t%e' % history_data[r][h])
                except:
                    fp.write('\t')
            fp.write('\n')

def get_attr_function(x):
    return lambda m: getattr(m, x)

def modified_output_function(pyomo_model, cv, mvs):
    return lambda m: pyomo_model.output_variables[cv+ '_unmodified'].__call__(m) + getattr(m, cv + "_eps") + \
                     sum([(pyomo_model.input_variables[mv].__call__(m) - getattr(m, mv + "_base")) *\
                          getattr(m,mv + "_" + cv + "_lam") for mv in mvs])
    # return lambda m:pyomo_model.output_variables[cv].__call__(m)+getattr(m,cv+"_eps")+\
    #                 sum([(pyomo_model.input_variables[mv].__call__(m)-getattr(m,mv+"_base"))*getattr(m,mv+"_"+cv+"_lam") for mv in mvs])

def get_hypercube_sampling(dimension, n_data, lb, ub, seed=1):
    engine = scipy.stats.qmc.LatinHypercube(d=dimension, seed=seed)
    samples = engine.random(n_data)
    sample_X = np.zeros((n_data,dimension))
    for i, s in enumerate(samples):
        for j in range(dimension):
            sample_X[i,j] = s[j] * (ub[j] - lb[j]) + lb[j]
    return sample_X

def square_circle_mapping(x, radius, circle_scaling):
    '''
    x in [-1,1]^n
    :param x: ndarray
    :param radius:
    :param circle_scaling: ndarray or list
    :return:
    '''
    dimension = x.shape[0]
    x_abs = np.zeros((dimension,))
    for i in range(dimension):
        x_abs[i] = abs(x[i])
    max_element = np.max(x_abs)
    if max_element < 1e-5:
        k_shrink = 1
    else:
        for i in range(dimension):
            x_abs[i] *= 1/max_element
        k_shrink = np.linalg.norm(x_abs)
    ret = np.zeros((dimension,))
    for i in range(dimension):
        ret[i] = x[i]/k_shrink*radius*circle_scaling[i]
    return ret

def sort_samples(X):
    from python_tsp.distances import great_circle_distance_matrix
    from python_tsp.heuristics import solve_tsp_simulated_annealing
    distance_matrix = great_circle_distance_matrix(X)
    permutation, distance = solve_tsp_simulated_annealing(distance_matrix, max_processing_time=500)
    return X[permutation,:]
