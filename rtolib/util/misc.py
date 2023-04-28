import scipy.stats.qmc
import numpy as np

def is_value_dict_same(dict1, dict2):
    for k in dict1:
        if abs(dict1[k]-dict2[k]) > 1e-6:
            return False
    return True


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