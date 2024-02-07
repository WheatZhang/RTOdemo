#!/usr/bin/env python
#-*- coding:utf-8 -*-
import numpy as np
from scipy.optimize import minimize
import re
from sklearn.linear_model import LinearRegression, Ridge

def GetVarSymbol(nth):
    return "x"+str(nth)
def GetFullXAsList(n):
    s = "["
    for i in range(n):
        if i != 0:
            s += ","
        s+=GetVarSymbol(i)
    s += "]"
    return s

class ConvergenceError(Exception):
    def __init__(self, str):
        print(str)

class MonomialBase(object):
    '''
    定义了多项式拟合所需要的每一个单项式（基）
    与MonomialFunction的区别是不含系数
    '''
    def __init__(self, power_coeff):
        self.power_coeff = power_coeff
        self.dimension = len(power_coeff)
        if self.dimension > 8:
            raise ValueError("The dimension is too high.")
        for i in self.power_coeff:
            if i > 10:
                raise ValueError("The order is too high.")

    def  __eq__(self, baseB):
        if not isinstance(baseB, MonomialBase):
            return False
        if self.dimension != baseB.dimension:
            return False
        for index,item in enumerate(self.power_coeff):
            if baseB.power_coeff[index] != item:
                return False
        return True

    def evaluate(self, sample):
        if len(sample) != self.dimension:
            raise ValueError("The dimension of the base is incorrect.")
        result = 1
        for i in range(self.dimension):
            for j in range(self.power_coeff[i]):
                result *= sample[i]
        return result

    def generate_str(self):
        str = ""
        first_flag = True
        for i in range(self.dimension):
            if self.power_coeff[i] == 0:
                continue
            if first_flag:
                first_flag = False
            else:
                str += "*"
            if self.power_coeff[i] == 1:
                str += GetVarSymbol(i)
            else:
                str += GetVarSymbol(i)+"_%d"%self.power_coeff[i]
        if first_flag:
            return "1"
        return str

class MonomialBaseGroup(object):
    '''
    定义了多项式拟合的一组完整的基。
    self.bases中存储了若干MonomialBase对象，每个对应一个单项式基。
    '''
    def __init__(self, name, dimension):
        self.name = name
        self.bases = []
        self.dimension = dimension
        self.max_order = [0 for i in range(dimension)]

    def numberOfBases(self):
        return len(self.bases)

    def add_base(self, base):
        if base.dimension != self.dimension:
            raise ValueError("The dimension of the base is incorrect.")
        for i in range(self.dimension):
            if base.power_coeff[i] > self.max_order[i]:
                self.max_order[i] = base.power_coeff[i]
        self.bases.append(base)

    def evaluate(self, sample, coeff):
        result = 0
        for index, item in enumerate(self.bases):
            result += item.evaluate(sample) * coeff[index]
        return result

    def generate_power_series(self):
        str = ""
        for i in range(self.dimension):
            str += MonomialBaseGroup.generate_single_power_series(GetVarSymbol(i), self.max_order[i])
        return str

    def generate_expression(self, coeff):
        str = ""
        str += "\toutput = \\\n"
        n_items = len(self.bases)
        for index, item in enumerate(self.bases):
            if index != n_items-1:
                str += "\t    "+item.generate_str()+"*%.25e+\\\n"%coeff[index]
            else:
                str += "\t    "+item.generate_str()+"*%.25e\n"%coeff[index]
        return str

    def __contains__(self, e):
        for i in self.bases:
            if i == e:
                return True
        return False

    @staticmethod
    def generate_single_power_series(var_symbol, order):
        '''
        This method is used in generating the target py function text.
        Like:
            x1_2=x1*x1
            x1_3=x1_2*x1
            x1_4=x1_3*x1
        :param var_symbol:the name of the variable, String
        :param order:Int. In the above example, the order is 4.
        :return:
        '''
        if order <= 1:
            return ""
        str = ""
        if order == 2:
            str += "\t"+var_symbol+"_2 = "+var_symbol+" * "+var_symbol+"\n"
            return str
        str += "\t"+var_symbol + "_2 = " + var_symbol + "*" + var_symbol+"\n"
        for i in range(3,order+1):
            str +="\tSYMBOL_%d = SYMBOL_%d * SYMBOL\n" %(i,i-1)
        str = str.replace("SYMBOL",var_symbol)
        return str

    @staticmethod
    def _base_A_iterater(order, dimension):
        if dimension == 1:
            for i in range(0, order + 1):
                yield [i]
        else:
            for i in range(0, order + 1):
                prefix = [i]
                for suffix in MonomialBaseGroup._base_A_iterater(order, dimension - 1):
                    yield prefix + suffix

    @staticmethod
    def _base_B_iterater(order, dimension):
        if dimension == 1:
            for i in range(0, order+1):
                yield [i]
        else:
            for i in range(0, order+1):
                prefix = [i]
                for suffix in MonomialBaseGroup._base_B_iterater(order-i, dimension-1):
                    yield prefix+suffix

    @staticmethod
    def general_fitting_base_A(order, dimension):
        # This function used to get a certain set of polynomial base. This is integrated in the funciton fitting_base_library.
        # e.g. dimension=2(x,y), order=3:
        # 1,x, x2,x3,
        # y,xy,x2y,x3y,
        # y2, xy2, x2y2,x3y2,
        # y3,xy3,x2y3,x3y3
        fitting_base_group = MonomialBaseGroup('TypeADim%dOrder%d' % (dimension, order), dimension)
        for item in MonomialBaseGroup._base_A_iterater(order, dimension):
            fitting_base_group.add_base((MonomialBase(item)))
        return fitting_base_group

    @staticmethod
    def general_fitting_base_B(order, dimension):
        # This function used to get a certain set of polynomial base. This is integrated in the funciton fitting_base_library.
        # e.g. dimension=2(x,y), order=3:
        # 1,x, x2,x3,
        # y,xy,x2y,
        # y2, xy2,
        # y3
        fitting_base_group = MonomialBaseGroup('TypeBDim%dOrder%d'%(dimension, order), dimension)
        for item in MonomialBaseGroup._base_B_iterater(order, dimension):
            fitting_base_group.add_base((MonomialBase(item)))
        return fitting_base_group

    @staticmethod
    def fitting_base_library(dimension, way, order):
        '''
        There are two types of base included in this library, individual and compound.
        A two dimensional individual base is:
            1, x, y, x2, xy, y2
        A two dimensional compound base is:
            1, x, y, xy
        In a word, an compound base consist of all the product of powers of independant variables that the order
        of monomial is less or equal to the disired order.
        A compound base consist of all the product of powers of independant variables that the maximum order of
        each independant variable in the polynomial is no greater than the disired order.
        Usually, the compound mode have relatively better performance.
        :param dimension:
        :param way: "Individual" or 'Compound'
        :param order:
        :return:
        '''
        if way == "Individual":
            return MonomialBaseGroup.general_fitting_base_A(order, dimension)
        elif way == 'Compound':
            return MonomialBaseGroup.general_fitting_base_B(order, dimension)
        else:
            raise ValueError("Unknown way.")


class ApproximatonFunction(object):
    def __init__(self, dimension):
        self.dimension = dimension

    def set_var_name(self, indep_var_name, dep_var_name):
        if self.dimension != len(indep_var_name):
            raise ValueError("The dimension of the indep_var_name is incorrect.")
        for name in indep_var_name:
            if re.match(r"x\d+", name) or re.match(r"x\d+_\d+", name) or name == 'y':
                raise NameError("Illegal independent variable name: %s"%name)
        self.indep_var_name = indep_var_name
        self.dep_var_name = dep_var_name

    def evaluate(self, x):
        return NotImplementedError

    def generate_text(self):
        return NotImplementedError

class PolynomialFunction(ApproximatonFunction):
    def __init__(self, dimension, base_group, coeff, max_value, min_value):
        self.dimension = dimension
        self.base_group = base_group
        self.coeff = coeff
        self.max_value = max_value
        self.min_value = min_value

    def evaluate(self, x):
        data = []
        for index,item in enumerate(x):
            data.append((item - self.min_value[index]) / (self.max_value[index] - self.min_value[index]))
        raw_output = self.base_group.evaluate(data, self.coeff)
        return raw_output * (self.max_value[-1] - self.min_value[-1]) + self.min_value[-1]

    def generate_text(self, name=None, flag_function_def=True):
        str = ""
        if flag_function_def:
            str += "def " + name + "("
            first_flag = True
            for i in range(self.dimension):
                if first_flag:
                    str += self.indep_var_name[i]
                    first_flag = False
                else:
                    str += ',' + self.indep_var_name[i]
            str += "):\n"
        for i in range(self.dimension):
            str += "\t" + GetVarSymbol(i) + " = (" + self.indep_var_name[i] + "-%.25e)/%.25e\n" % (
                self.min_value[i], self.max_value[i] - self.min_value[i])
        str += self.base_group.generate_power_series()
        str += self.base_group.generate_expression(self.coeff)
        str += "\toutput = output*%.25e+%.25e\n" % (
            self.max_value[-1] - self.min_value[-1], self.min_value[-1])
        if flag_function_def:
            str += "\t" + self.dep_var_name + " = output\n"
            str += "\treturn " + self.dep_var_name+"\n\n"
        return str


class PolynomialFitter(object):
    def __init__(self, dimension, base_group):
        assert isinstance(base_group, MonomialBaseGroup)
        if dimension != base_group.dimension:
            return ValueError("The dimension of the base is incorrect.")
        self.base_group = base_group
        self.dimension = dimension
        self.X_center = None
        self.y_center = None
        self.original_X_center = None

    def set_center(self, x, y):
        '''
        :param x: A list or an ndarray
        :param y: A float
        :return:
        '''
        self.original_X_center = np.array(x)
        self.original_y_center = y

    def cost_function(self, params, X, y):
        error = y - X.dot(params)
        return error.dot(error)

    def center_constraints(self, params):
        error = self.y_center - self.X_center.dot(params)
        return error

    def fit(self, data_X, data_Y, init_coeff=None, method='minimize'):
        if method == 'sklearn-ols' and self.original_X_center is not None:
            raise Exception("sklearn-ols method can not handle centered ")
        if data_X.shape[1]!=self.dimension:
            return ValueError("The dimension of the data is incorrect.")
        min_value = np.zeros(shape=(self.dimension+1,))
        max_value = np.zeros(shape=(self.dimension+1,))
        for i in range(self.dimension):
            min_value[i] = np.min(data_X[:, i])
            max_value[i] = np.max(data_X[:, i])
            if max_value[i] == min_value[i]:
                raise ValueError("input %d"%i + " has only one value in all samples.")
        min_value[self.dimension] = np.min(data_Y[:])
        max_value[self.dimension] = np.max(data_Y[:])
        data_X_scaled = np.zeros(shape=data_X.shape)
        data_Y_scaled = np.zeros(shape=data_Y.shape)
        n_data = data_X_scaled.shape[0]
        # scale raw data
        for i in range(self.dimension):
            data_X_scaled[:, i] = (data_X[:, i] - min_value[i]) / (
                            max_value[i] - min_value[i])
        data_Y_scaled = (data_Y-min_value[self.dimension]) / (
                            max_value[self.dimension] - min_value[self.dimension])
        # adapt data to the base group
        temp = []
        for n in range(n_data):
            temp.append(self.base_group.bases[0].evaluate(data_X_scaled[n, :]))
        data_for_fit = np.transpose([temp])
        for index, item in enumerate(self.base_group.bases):
            if index == 0:
                continue
            temp = []
            for n in range(n_data):
                temp.append(item.evaluate(data_X_scaled[n, :]))
            data_for_fit = np.hstack((data_for_fit,np.transpose([temp])))
        # scale and preceed the center if needed
        if self.original_X_center is not None:
            original_X_center = np.zeros(shape=(1,self.original_X_center.shape[0]))
            for index, item in enumerate(self.original_X_center):
                original_X_center[0,index] = (self.original_X_center[index] - min_value[index]) / (
                            max_value[index] - min_value[index])
            self.X_center = np.zeros(shape=(self.base_group.numberOfBases(),))
            for i in range(self.X_center.shape[0]):
                self.X_center[i] = self.base_group.bases[i].evaluate(original_X_center[0,self.indep_var_index])
            self.y_center = (self.original_y_center - min_value[self.dimension]) / (
                            max_value[self.dimension] - min_value[self.dimension])
            # print(self.X_center)
            # print(self.y_center)
        if method=='minimize':
            # get initial values for the fitting parameters
            if init_coeff == None:
                x0 = np.zeros(shape=(self.base_group.numberOfBases(),))
            else:
                x0 = init_coeff
            # start fitting
            if self.X_center is not None:
                # if set center constraints
                output = minimize(self.cost_function, x0, args=(data_for_fit, data_Y_scaled),\
                                  constraints={'type': 'eq', 'fun':self.center_constraints},options={'disp': False})
            else:
                output = minimize(self.cost_function, x0, args=(data_for_fit, data_Y_scaled),options={'disp': False})
            coeff = output.x
        elif method == 'sklearn-ols':
            ls_fitter = LinearRegression(fit_intercept=False)
            ls_fitter.fit(data_for_fit, data_Y_scaled)
            coeff = ls_fitter.coef_

        # return the result
        return PolynomialFunction(self.dimension, self.base_group,
                                  coeff, max_value, min_value)


class GaussianKernel():
    def __init__(self, dimension, bandwidth):
        '''

        :param dimension:
        :param bandwidth: list of size 'dimension'
        '''
        self.dimension = dimension
        self.bandwidth = bandwidth

    def evaluate_single(self, input_value):
        if input_value.shape[0] != self.dimension:
            raise ValueError("The dimension is not correct")
        ret = 1
        for d in range(self.dimension):
            ret *= np.exp(-(input_value[d]**2)/2/self.bandwidth[d]/self.bandwidth[d])
        return ret

    def evaluate(self, input_value1, input_value2):
        return self.evaluate_single(input_value1-input_value2)


class GaussianRBFNetworkFunc(ApproximatonFunction):
    def __init__(self, dimension, bandwidth, x_points, coeff, offset=0):
        '''

        :param dimension:
        :param bandwidth: list of size 'dimension'
        '''
        self.dimension = dimension
        self.bandwidth = bandwidth
        self.x_points = x_points
        self.data_size = self.x_points.shape[0]
        self.kernel = GaussianKernel(dimension, bandwidth)
        self.coeff = coeff
        self.offset = offset

    def evaluate(self, x):
        '''

        :param data_point: tuple (input,) like (1,3)
        :return:
        '''
        if x.shape[0] != self.dimension:
            raise ValueError("The dimension is not correct")
        ret = 0
        for d in range(self.data_size):
            ret += self.kernel.evaluate(x,self.x_points[d])*self.coeff[d]
        return ret + self.offset

    def generate_text(self, name=None, flag_function_def=True):
        str = ""
        if flag_function_def:
            str += "def " + name + "("
            first_flag = True
            for i in range(self.dimension):
                if first_flag:
                    str += self.indep_var_name[i]
                    first_flag = False
                else:
                    str += ',' + self.indep_var_name[i]
            str += "):\n"
        name_list="["
        name_list+=",".join(self.indep_var_name)
        name_list+="]"
        str += "\tfull_x_list = "+name_list+"\n"
        str += "\tdata_x = np.array("+np.array2string(self.x_points,separator=',', formatter={'float_kind':lambda x: "%.25e" % x})+")\n"
        str += "\tcoeff = np.array(" + np.array2string(self.coeff,separator=',', formatter={'float_kind': lambda x: "%.25e" % x}) + ")\n"
        str += "\tbandwidth = np.array(" + np.array2string(self.bandwidth,separator=',',
                                                       formatter={'float_kind': lambda x: "%.25e" % x}) + ")\n"
        str +=\
'''
\toutput = 0
\tfor i in range(%d):
\t\ttemp = 0
\t\tfor j in range(%d):
\t\t\ttemp += ((full_x_list[j]-data_x[i,j])/bandwidth[j])**2
\t\ttemp *= -0.5
\t\toutput += np.exp(temp)*coeff[i]
'''%(self.data_size, self.dimension)
        str += "\toutput += %.25e\n"%self.offset
        if flag_function_def:
            str += "\t" + self.dep_var_name + " = output\n"
            str += "\treturn " + self.dep_var_name + "\n\n"
        return str

class GaussianRBFNetworkFuncFitter():
    def __init__(self, dimension, bandwidth):
        self.dimension = dimension
        self.bandwidth = bandwidth
        self.kernel = GaussianKernel(dimension, bandwidth)

    def fit(self, sample_X, sample_Y, point_X=None, fit_intercept=False, alpha=0.0):
        '''

        :param sample_X:
        :param sample_Y:
        :param point_X:
        :param fit_intercept:bool
        :param alpha: [0,inf)
        :return:
        '''
        if point_X is None:
            point_X = sample_X
        n_p = point_X.shape[0]
        n_s = sample_X.shape[0]
        if n_p > n_s:
            raise ValueError("too few sample points")

        X = np.zeros((n_s, n_p))
        for i in range(n_s):
            for j in range(n_p):
                X[i,j] = self.kernel.evaluate(sample_X[i],point_X[j])
        if alpha == 0.0:
            model = LinearRegression(fit_intercept=fit_intercept)
        else:
            model = Ridge(alpha=1.0, fit_intercept=fit_intercept)

        model.fit(X, sample_Y)
        coeff = np.array(model.coef_)
        offset = model.intercept_
        return GaussianRBFNetworkFunc(self.dimension, self.bandwidth, point_X, coeff, offset)


class KrigingFunction(ApproximatonFunction):
    def __init__(self, polynomial, rbf_network):
        assert polynomial.dimension == rbf_network.dimension
        self.dimension = polynomial.dimension
        self.polynomial = polynomial
        self.rbf_network = rbf_network

    def evaluate(self, x):
        return self.polynomial.evaluate(x)+self.rbf_network.evaluate(x)

    def generate_text(self, name=None, flag_function_def=True):
        str = ""
        if flag_function_def:
            str += "def " + name + "("
            first_flag = True
            for i in range(self.dimension):
                if first_flag:
                    str += self.indep_var_name[i]
                    first_flag = False
                else:
                    str += ',' + self.indep_var_name[i]
            str += "):\n"
        str += self.polynomial.generate_text(flag_function_def=False)
        str += "\t"+self.dep_var_name+" = output\n"
        str += self.rbf_network.generate_text(flag_function_def=False)
        str += "\t" + self.dep_var_name + " += output\n"
        if flag_function_def:
            str += "\treturn " + self.dep_var_name + "\n\n"
        return str

    def set_var_name(self, indep_var_name, dep_var_name):
        self.polynomial.set_var_name(indep_var_name, dep_var_name)
        self.rbf_network.set_var_name(indep_var_name, dep_var_name)
        super().set_var_name(indep_var_name, dep_var_name)

class KrigingFitter():
    def __init__(self, dimension, base_group, bandwidth):
        self.polynomial_fitter = PolynomialFitter(dimension, base_group)
        self.rbf_network_fitter = GaussianRBFNetworkFuncFitter(dimension, bandwidth)
    def fit(self, sample_X, sample_Y, point_X=None, fit_intercept=False, alpha=0.0):
        polynomial = self.polynomial_fitter.fit(sample_X, sample_Y, init_coeff=None, method='sklearn-ols')
        polynomial_Y = np.zeros(shape=sample_Y.shape)
        for i in range(polynomial_Y.shape[0]):
            polynomial_Y[i] = polynomial.evaluate(sample_X[i])
        rbf_network = self.rbf_network_fitter.fit(sample_X, sample_Y-polynomial_Y, point_X, fit_intercept, alpha)
        return KrigingFunction(polynomial, rbf_network)


if __name__=="__main__":
    pass