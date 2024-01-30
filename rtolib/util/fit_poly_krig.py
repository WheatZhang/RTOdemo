#!/usr/bin/env python
#-*- coding:utf-8 -*-
import numpy as np
from scipy.optimize import minimize
import re
from sklearn.linear_model import LinearRegression

def GetVarSymbol(nth):
    return "x"+str(nth)

class ConvergenceError(Exception):
    def __init__(self, str):
        print(str)

class MonomialItem(object):
    def __init__(self, power_coeff, coeff):
        self.power_coeff = power_coeff
        self.coeff = coeff

    def calc_derivative(self, direction):
        if self.power_coeff[direction] == 0:
            return MonomialItem([0,0,0], 0)
        else:
            coeff=self.coeff*self.power_coeff[direction]
            power_coeff = [0,0,0]
            for i in range(3):
                power_coeff[i] = self.power_coeff[i]
            power_coeff[direction] -= 1
            return MonomialItem(power_coeff, coeff)

class fittingBase(object):
    '''
    定义了多项式拟合所需要的每一个单项式（基）
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
        if not isinstance(baseB, fittingBase):
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

class fittingBaseGroup(object):
    '''
    定义了多项式拟合的一组完整的基。
    self.bases中存储了若干fittingBase对象，每个对应一个单项式基。
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
            str += fittingBaseGroup.generate_single_power_series(GetVarSymbol(i), self.max_order[i])
        return str

    def generate_expression(self, coeff):
        str = ""
        str += "\toutput = \\\n"
        n_items = len(self.bases)
        for index, item in enumerate(self.bases):
            if index != n_items-1:
                str += "\t    "+item.generate_str()+"*%.8e+\\\n"%coeff[index]
            else:
                str += "\t    "+item.generate_str()+"*%.8e\n"%coeff[index]
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
                for suffix in fittingBaseGroup._base_A_iterater(order, dimension - 1):
                    yield prefix + suffix

    @staticmethod
    def _base_B_iterater(order, dimension):
        if dimension == 1:
            for i in range(0, order+1):
                yield [i]
        else:
            for i in range(0, order+1):
                prefix = [i]
                for suffix in fittingBaseGroup._base_B_iterater(order-i, dimension-1):
                    yield prefix+suffix

    @staticmethod
    def general_fitting_base_A(order, dimension):
        # This function used to get a certain set of polynomial base. This is integrated in the funciton fitting_base_library.
        # e.g. dimension=2(x,y), order=3:
        # 1,x, x2,x3,
        # y,xy,x2y,x3y,
        # y2, xy2, x2y2,x3y2,
        # y3,xy3,x2y3,x3y3
        fitting_base_group = fittingBaseGroup('TypeADim%dOrder%d' % (dimension, order), dimension)
        for item in fittingBaseGroup._base_A_iterater(order, dimension):
            fitting_base_group.add_base((fittingBase(item)))
        return fitting_base_group

    @staticmethod
    def general_fitting_base_B(order, dimension):
        # This function used to get a certain set of polynomial base. This is integrated in the funciton fitting_base_library.
        # e.g. dimension=2(x,y), order=3:
        # 1,x, x2,x3,
        # y,xy,x2y,
        # y2, xy2,
        # y3
        fitting_base_group = fittingBaseGroup('TypeBDim%dOrder%d'%(dimension, order), dimension)
        for item in fittingBaseGroup._base_B_iterater(order, dimension):
            fitting_base_group.add_base((fittingBase(item)))
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
            return fittingBaseGroup.general_fitting_base_A(order, dimension)
        elif way == 'Compound':
            return fittingBaseGroup.general_fitting_base_B(order, dimension)
        else:
            raise ValueError("Unknown way.")

class GeneralFitter(object):
    def __init__(self, name, data=None, indep_var_index=None, dep_var_index=None, base_group=None, sign_derivative=False):
        if not sign_derivative:
            if data is None:
                return
            if not isinstance(data, np.ndarray):
                raise TypeError("The data for fitting process must be an ndarray.")
            self.name = name
            self.base_group = base_group
            self.data = data
            self.dimension = len(indep_var_index)
            if self.dimension != self.base_group.dimension:
                raise ValueError("The dimension of the base is incorrect.")
            self.indep_var_index =indep_var_index
            self.dep_var_index = dep_var_index
            self.X_center = None
            self.y_center = None
            self.original_X_center = None
        else:
            self.name=name

    def set_center(self, x, y):
        '''
        :param x: A list or an ndarray
        :param y: A float
        :return:
        '''
        self.original_X_center = np.array(x)
        self.original_y_center = y

    def set_var_name(self, indep_var_name, dep_var_name):
        if self.dimension != len(indep_var_name):
            raise ValueError("The dimension of the indep_var_name is incorrect.")
        for name in indep_var_name:
            if re.match(r"x\d+_\d+", name) or name == 'y':
                raise NameError("Illegal independent variable name: %s"%name)
        self.indep_var_name = indep_var_name
        self.dep_var_name = dep_var_name

    def cost_function_OLD_VERSION(self, params, data):
        error = 0
        for n_data in range(data.shape[0]):
            error+=np.linalg.norm(data[n_data,self.dep_var_index] - \
                                  self.base_group.evaluate(data[n_data, self.indep_var_index],params))
        return error

    def fit_OLD_VERSION(self, fitterB=None):
        self.min_value = np.zeros(shape=(self.dimension+1,))
        self.max_value = np.zeros(shape=(self.dimension+1,))
        for i in range(self.dimension):
            self.min_value[i] = np.min(self.data[:, self.indep_var_index[i]])
            self.max_value[i] = np.max(self.data[:, self.indep_var_index[i]])
        self.min_value[self.dimension] = np.min(self.data[:, self.dep_var_index])
        self.max_value[self.dimension] = np.max(self.data[:, self.dep_var_index])
        data = np.zeros(shape=self.data.shape)
        for i in range(self.data.shape[0]):
            for index,item in enumerate(self.indep_var_index):
                data[i,item] = (self.data[i,item]-self.min_value[index])/(self.max_value[index]-self.min_value[index])
            data[i, self.dep_var_index] = (self.data[i, self.dep_var_index] - self.min_value[-1]) / (self.max_value[-1] - self.min_value[-1])
        if fitterB == None:
            x0 = np.zeros(shape=(self.base_group.numberOfBases(),))
        else:
            x0 = np.zeros(shape=(self.base_group.numberOfBases(),))
            for index, item in enumerate(fitterB.base_group.bases):
                for index2, item2 in enumerate(self.base_group.bases):
                    if item == item2:
                        self.coeff[item2] = fitterB.coeff[item]
                        break
        output = minimize(self.cost_function, x0, args=(data))
        self.coeff = output.x

    def cost_function(self, params, X, y):
        error = y - X.dot(params)
        return error.dot(error)

    def center_constraints(self, params):
        error = self.y_center - self.X_center.dot(params)
        return error

    def fit(self, fitterB=None, method='minimize'):
        if method == 'sklearn-ols' and self.original_X_center is not None:
            raise Exception("sklearn-ols method can not handle centered ")
        self.min_value = np.zeros(shape=(self.dimension+1,))
        self.max_value = np.zeros(shape=(self.dimension+1,))
        for i in range(self.dimension):
            self.min_value[i] = np.min(self.data[:, self.indep_var_index[i]])
            self.max_value[i] = np.max(self.data[:, self.indep_var_index[i]])
        self.min_value[self.dimension] = np.min(self.data[:, self.dep_var_index])
        self.max_value[self.dimension] = np.max(self.data[:, self.dep_var_index])
        data = np.zeros(shape=self.data.shape)
        # scale raw data
        for i in range(self.data.shape[0]):
            for index,item in enumerate(self.indep_var_index):
                if self.max_value[index] == self.min_value[index]:
                    raise ValueError(self.indep_var_name[index]+" has only one value in all samples.")
                data[i,item] = (self.data[i,item]-self.min_value[index])/(self.max_value[index]-self.min_value[index])
            if self.max_value[-1] == self.min_value[-1]:
                raise ValueError(self.dep_var_name+" has only one value in all samples.")
            data[i, self.dep_var_index] = (self.data[i, self.dep_var_index] - self.min_value[-1]) / (self.max_value[-1] - self.min_value[-1])
        # adapt data to the base group. This makes the optimization problem an SQP.
        temp = []
        for n_data in range(data.shape[0]):
            temp.append(self.base_group.bases[0].evaluate(data[n_data, self.indep_var_index]))
        data_for_fit = np.transpose([temp])
        for index, item in enumerate(self.base_group.bases):
            if index == 0:
                continue
            temp = []
            for n_data in range(data.shape[0]):
                temp.append(item.evaluate(data[n_data, self.indep_var_index]))
            data_for_fit = np.hstack((data_for_fit,np.transpose([temp])))
        # scale and preceed the center if needed
        if self.original_X_center is not None:
            original_X_center = np.zeros(shape=(1,self.original_X_center.shape[0]))
            for index, item in enumerate(self.original_X_center):
                original_X_center[0,index] = (self.original_X_center[index] - self.min_value[index]) / (
                            self.max_value[index] - self.min_value[index])
            self.X_center = np.zeros(shape=(self.base_group.numberOfBases(),))
            for i in range(self.X_center.shape[0]):
                self.X_center[i] = self.base_group.bases[i].evaluate(original_X_center[0,self.indep_var_index])
            self.y_center = (self.original_y_center - self.min_value[self.dimension]) / (
                            self.max_value[self.dimension] - self.min_value[self.dimension])
            # print(self.X_center)
            # print(self.y_center)
        if method=='minimize':
            # get initial values for the fitting parameters
            if fitterB == None:
                x0 = np.zeros(shape=(self.base_group.numberOfBases(),))
            else:
                x0 = np.zeros(shape=(self.base_group.numberOfBases(),))
                for index, item in enumerate(fitterB.base_group.bases):
                    for index2, item2 in enumerate(self.base_group.bases):
                        if item == item2:
                            x0[index2] = fitterB.coeff[index]
                            break
            # start fitting
            if self.X_center is not None:
                # if set center constraints
                output = minimize(self.cost_function, x0, args=(data_for_fit, data[:, self.dep_var_index]),\
                                  constraints={'type': 'eq', 'fun':self.center_constraints},options={'disp': False})
            else:
                output = minimize(self.cost_function, x0, args=(data_for_fit, data[:, self.dep_var_index]),options={'disp': False})
            self.coeff = output.x
        elif method == 'sklearn-ols':
            ls_fitter = LinearRegression(fit_intercept=False)
            ls_fitter.fit(data_for_fit, data[:, self.dep_var_index])
            self.coeff = ls_fitter.coef_

    def evaluate(self, sample):
        if not hasattr(self, "coeff"):
            raise AttributeError("Must do regression first.")
        data = []
        for index,item in enumerate(sample):
            data.append((item - self.min_value[index]) / (self.max_value[index] - self.min_value[index]))
        raw_output = self.base_group.evaluate(data, self.coeff)
        return raw_output * (self.max_value[-1] - self.min_value[-1]) + self.min_value[-1]

    def generate_result(self):
        if not hasattr(self, "coeff"):
            raise AttributeError("Must do regression first.")
        str = ""
        str += "def " + self.name + "("
        first_flag = True
        for i in range(self.dimension):
            if first_flag:
                str += self.indep_var_name[i]
                first_flag = False
            else:
                str += ',' + self.indep_var_name[i]
        str += "):\n"
        for i in range(self.dimension):
            str += "\t" + GetVarSymbol(i) + " = (" + self.indep_var_name[i] + "-%.8e)/%.8e\n" % (
            self.min_value[i], self.max_value[i] - self.min_value[i])
        str += self.base_group.generate_power_series()
        str += self.base_group.generate_expression(self.coeff)
        str += "\t" + self.dep_var_name + " = output*%.8e+%.8e\n" % (
        self.max_value[-1] - self.min_value[-1], self.min_value[-1])
        str += "\treturn " + self.dep_var_name+"\n\n"
        return str

    def get_derivative_fitter(self, partial_var_index):
        if partial_var_index>= self.dimension:
            raise IndexError("The partial derivative variable index beyond the scope.")
        derivative_fitter = GeneralFitter(name = self.name+"_p"+self.indep_var_name[partial_var_index], sign_derivative=True)
        derivative_fitter.indep_var_name = self.indep_var_name
        derivative_fitter.dep_var_name = self.dep_var_name
        derivative_fitter.dimension = self.dimension
        derivative_fitter.base_group = fittingBaseGroup("DerivativeFitterBase",self.dimension)
        derivative_fitter.coeff = []
        derivative_fitter.max_value = self.max_value.copy()
        derivative_fitter.min_value = self.min_value.copy()
        derivative_fitter.max_value[self.dimension] = 1
        derivative_fitter.min_value[self.dimension] = 0
        for i in range(len(self.base_group.bases)):
            if self.base_group.bases[i].power_coeff[partial_var_index] > 0:
                power_coeff = []
                for index, item in enumerate(self.base_group.bases[i].power_coeff):
                    if index == partial_var_index:
                        power_coeff.append(item-1)
                    else:
                        power_coeff.append(item)
                derivative_fitter.base_group.add_base(fittingBase(power_coeff))
                derivative_fitter.coeff.append(self.base_group.bases[i].power_coeff[partial_var_index]*self.coeff[i])
        for i in range(len(derivative_fitter.coeff)):
            derivative_fitter.coeff[i] *= (self.max_value[self.dimension]-self.min_value[self.dimension])/\
                                       (self.max_value[partial_var_index]-self.min_value[partial_var_index])
        return derivative_fitter

    def fit_and_write(self, fitterB=None, method = 'minimize'):
        self.fit(fitterB, method)
        return self.generate_result()

class BasePolynomialFitting(object):
    def __init__(self, thermo_name, method = 'minimize'):
        self.thermo_name = thermo_name
        self.X_center = None
        self.y_center = None
        self.nominal_value = None
        self.deviation = None
        self.method = method

    def set_data(self, df, names):
        self.csv_data = df
        self.name_in_df = names
        self.indep_var_names = names[:-1]
        self.dep_var_name = names[-1]
        self.dimension = len(self.indep_var_names)

    def set_range(self, nominal_value, deviation):
        self.nominal_value = nominal_value
        self.deviation = deviation

    def set_center(self, X_center, y_center):
        self.X_center = X_center
        self.y_center = y_center

    def get_indicator_func(self, df):
        result = (df[self.indep_var_names[0]]>(self.nominal_value[0]-self.deviation[0])) & (df[self.indep_var_names[0]]<(self.nominal_value[0]+self.deviation[0]))
        if self.dimension == 1:
            return result
        for i in range(1,len(self.indep_var_names)):
            result = result & (df[self.indep_var_names[i]]>(self.nominal_value[i]-self.deviation[i])) & (df[self.indep_var_names[i]]<(self.nominal_value[i]+self.deviation[i]))
        return result

class AdaptivePolynomialFitting(BasePolynomialFitting):
    '''
    choose a proper order adaptively, starting from a low order
    '''
    def set_fitting_param(self, way, max_order, desired_error):
        self.max_order = max_order
        self.way = way
        self.target_error = desired_error

    def generate_derivative_func(self, indexOrName):
        if not hasattr(self, "last_fitter"):
            raise AttributeError("The fitting should be done first.")
        if isinstance(indexOrName, str):
            index = self.indep_var_names.index(indexOrName)
        else:
            index = indexOrName
        der_fitter = self.last_fitter.get_derivative_fitter(index)
        return der_fitter.generate_result()

    def generate_funcion(self):
        if self.nominal_value is not None:
            selected_data = self.csv_data[self.get_indicator_func(self.csv_data)]
        else:
            selected_data = self.csv_data
        selected_data = selected_data.sample(frac=1.0)  # 全部打乱
        cut_idx = int(round(0.5 * selected_data.shape[0]))
        df_test, df_train = selected_data.iloc[:cut_idx], selected_data.iloc[cut_idx:]
        if selected_data.shape[0] < self.max_order*self.dimension:
            raise ValueError("Too few samples.At least %d."%(self.max_order*self.dimension))
        # print(selected_data.shape[0])
        # print(df_test.shape[0])
        for name in self.name_in_df:
            if name not in df_train.columns:
                raise KeyError("The dataset has no key called %s." % name)
        train_data = np.array(df_train.loc[:,self.name_in_df])
        test_data = np.array(df_test.loc[:,self.name_in_df])
        indep_var_index = []
        for i in range(len(self.indep_var_names)):
            indep_var_index.append(i)
        dep_var_index = self.dimension

        last_base = fittingBaseGroup.fitting_base_library(self.dimension, self.way, 1)
        last_fitter = GeneralFitter(self.thermo_name, train_data, indep_var_index, dep_var_index, last_base)
        last_fitter.set_var_name(self.indep_var_names,self.dep_var_name)
        if self.X_center is not None:
            last_fitter.set_center(self.X_center, self.y_center)
        last_fitter.fit(method=self.method)
        last_predict = np.zeros(shape = (test_data.shape[0],1))
        for i in range(test_data.shape[0]):
            last_predict[i,0] = last_fitter.evaluate(test_data[i,:-1])
        last_error = last_predict-test_data[:,[-1]]
        max_last_error = max(abs(last_error))[0]
        flag = False
        order = 1
        if max_last_error>self.target_error:
            if self.max_order > 1:
                for i in range(2,self.max_order+1):
                    new_base = fittingBaseGroup.fitting_base_library(self.dimension, self.way, i)
                    new_fitter = GeneralFitter(self.thermo_name, train_data, indep_var_index, dep_var_index, new_base)
                    new_fitter.set_var_name(self.indep_var_names, self.dep_var_name)
                    if self.X_center is not None:
                        new_fitter.set_center(self.X_center, self.y_center)
                    new_fitter.fit(last_fitter,method=self.method)
                    new_predict = np.zeros(shape=(test_data.shape[0], 1))
                    for j in range(test_data.shape[0]):
                        new_predict[j, 0] = new_fitter.evaluate(test_data[j, :-1])
                    new_error = new_predict - test_data[:, [-1]]
                    max_new_error = max(abs(new_error))[0]
                    self.current_mean_error = np.mean(abs(new_error))
                    if max_new_error > max_last_error:
                        raise ConvergenceError("Mission impossible, the data should be more localized.")
                    else:
                        last_fitter = new_fitter
                        max_last_error = max_new_error
                        order = i
                        if max_last_error < self.target_error:
                            flag = True
                            break
            else:
                raise ConvergenceError("Mission impossible, the order should be bigger.")
        else:
            flag = True
        self.current_max_error = max_last_error
        self.order = order
        self.last_fitter = last_fitter
        if flag == False:
            raise ConvergenceError("Mission impossible, the order should be bigger.")
        print("The max error on test dataset is %e." % max_last_error)
        return last_fitter.generate_result()

class FixOrderPolynomialFitting(BasePolynomialFitting):
    '''
    choose a proper order adaptively, starting from a low order
    '''
    def set_fitting_param(self, way, order, desired_error):
        self.order = order
        self.way = way
        self.target_error = desired_error

    def generate_derivative_func(self, indexOrName):
        if not hasattr(self, "this_fitter"):
            raise AttributeError("The fitting should be done first.")
        if isinstance(indexOrName, str):
            index = self.indep_var_names.index(indexOrName)
        else:
            index = indexOrName
        der_fitter = self.this_fitter.get_derivative_fitter(index)
        return der_fitter.generate_result()

    def generate_funcion(self):
        '''
        Revised 20190917. High order fitting using low order fitting result as init values.
        :return:
        '''
        # the error here means absolute error
        if self.nominal_value is not None:
            selected_data = self.csv_data[self.get_indicator_func(self.csv_data)]
        else:
            selected_data = self.csv_data
        selected_data = selected_data.sample(frac=1.0)  # 全部打乱
        cut_idx = int(round(0.5 * selected_data.shape[0]))
        df_test, df_train = selected_data.iloc[:cut_idx], selected_data.iloc[cut_idx:]
        if selected_data.shape[0] < self.order*self.dimension:
            raise ValueError("Too few samples.At least %d, now it is %d."%(self.order*self.dimension,selected_data.shape[0]))
        # print(selected_data.shape[0])
        # print(df_test.shape[0])
        for name in self.name_in_df:
            if name not in df_train.columns:
                raise KeyError("The dataset has no key called %s." % name)
        train_data = np.array(df_train.loc[:,self.name_in_df])
        test_data = np.array(df_test.loc[:,self.name_in_df])
        indep_var_index = []
        for i in range(len(self.indep_var_names)):
            indep_var_index.append(i)
        dep_var_index = self.dimension

        for fitting_order in range(1, self.order +1):
            base = fittingBaseGroup.fitting_base_library(self.dimension, self.way, fitting_order)
            this_fitter = GeneralFitter(self.thermo_name, train_data, indep_var_index, dep_var_index, base)
            this_fitter.set_var_name(self.indep_var_names, self.dep_var_name)
            if self.X_center is not None:
                this_fitter.set_center(self.X_center, self.y_center)
            if fitting_order == 1:
                this_fitter.fit(method=self.method)
                last_fitter = this_fitter
            else:
                this_fitter.fit(last_fitter,method=self.method)
        this_predict = np.zeros(shape = (test_data.shape[0],1))
        for i in range(test_data.shape[0]):
            this_predict[i,0] = this_fitter.evaluate(test_data[i,:-1])
        this_error = this_predict-test_data[:,[-1]]
        max_error = max(abs(this_error))[0]
        print("The max error on test dataset is %e."%max_error)
        self.current_mean_error = np.mean(abs(this_error))
        if max_error > self.target_error:
            print("Target precision not reached.")
        self.current_max_error = max_error
        self.this_fitter = this_fitter
        return this_fitter.generate_result()

if __name__=="__main__":
    for item in fittingBaseGroup._base_A_iterater(3,4):
        print(item)