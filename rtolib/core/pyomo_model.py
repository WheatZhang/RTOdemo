#!/usr/bin/env python
#-*- coding:utf-8 -*-
from rtolib.core.basic import Simulator, Optimizer, ParameterEstimator, ProblemDescription
from pyomo.environ import *
from pyomo.opt import SolverStatus, TerminationCondition
import rtolib.util.init_value as init_value
from rtolib.util.misc import get_attr_function, modified_output_function
import copy
from enum import Enum

__name__ = ['PyomoModelSolvingStatus','PyomoModel',
           'ModifierType','PyomoModelWithModifiers']

class ModifierType(Enum):
    RTO = 1
    OUTPUT = 2
    BOTH = 3


class PyomoModel():
    def __init__(self):
        '''
        output_variables: key:string, value: a function that takes a pyomo ConcreteModel as an input, and
                                            returns a pyomo expr
        input_variables: same type as output_variables. input_variables include MV and SPEC
        parameters: same type as output_variables
        objective: a string, supposed to be minimized
        default_value: key:string  value: double, insistent with the initial value file
        parameter_scaling_factors: key:string  value: double
        '''
        self.output_variables = {}
        self.parameters={}
        self.input_variables = {}
        self.default_value = {}
        self.parameter_scaling_factors = {}
        self.initial_value_file=""

    def build_body(self, model):
        raise NotImplementedError('This is an interface method.')

    def build_rto(self, model, cv_func):
        raise NotImplementedError('This is an interface method.')

    def build(self, model):
        self.build_body(model)
        self.build_rto(model, self.output_variables)

    def load_init(self, model):
        raise Exception('This is an interface method.')


class PyomoModelWithModifiers(PyomoModel):
    def __init__(self, pyomo_model, modifier_type, mvs, cvs):
        '''

        :param pyomo_model: PyomoModel
        :param mvs: list of strings
        :param cvs: list of strings
        :return:
        '''
        assert isinstance(pyomo_model, PyomoModel)
        self.modifier_mvs = mvs
        self.modifier_cvs = cvs
        self.modifier_type = modifier_type
        self.output_variables = {}
        for k, v in pyomo_model.output_variables.items():
            if k in self.modifier_cvs:
                self.output_variables[k + '_unmodified'] = pyomo_model.output_variables[k]
                self.output_variables[k] = get_attr_function(k + "_modified")
            else:
                self.output_variables[k] = pyomo_model.output_variables[k]
        self.parameters = copy.deepcopy(pyomo_model.parameters)
        self.input_variables = pyomo_model.input_variables
        self.default_value = copy.deepcopy(pyomo_model.default_value)
        self.parameter_scaling_factors = copy.deepcopy(pyomo_model.parameter_scaling_factors)
        self.initial_value_file = pyomo_model.initial_value_file
        self.base_pyomo_model = pyomo_model

        for name in self.modifier_names_iterator():
           self.parameters[name] = get_attr_function(name)
           self.default_value[name] = 0
           # TODO:需要更加科学的确定方法
           self.parameter_scaling_factors[name] = 1


    def modifier_names_iterator(self):
        for cv in self.modifier_cvs:
            yield cv + "_eps"
            for mv in self.modifier_mvs:
                yield mv+"_"+cv + "_lam"

    def build(self, model):
        self.base_pyomo_model.build_body(model)
        for name in self.modifier_names_iterator():
            setattr(model, name, Var(initialize=0))
        for mv in self.modifier_mvs:
            setattr(model, mv + "_base", Var(initialize=0))
            getattr(model, mv + "_base").fixed = True
        # TODO: need to test if there are bugs
        if self.modifier_type == ModifierType.RTO:
            self.base_pyomo_model.build_rto(model, self.output_variables)
            for cv in self.modifier_cvs:
                temp_func = modified_output_function(self, cv, self.modifier_mvs)
                setattr(model, cv + "_modified", Expression(rule=temp_func))
        elif self.modifier_type == ModifierType.OUTPUT:
            for cv in self.modifier_cvs:
                temp_func=modified_output_function(self, cv, self.modifier_mvs)
                setattr(model, cv+"_modified", Expression(rule=temp_func))
            self.base_pyomo_model.build_rto(model, self.output_variables)
        elif self.modifier_type == ModifierType.BOTH:
            self.base_pyomo_model.build_rto(model, self.base_pyomo_model.output_variables)
            for cv in self.modifier_cvs:
                temp_func=modified_output_function(self, cv, self.modifier_mvs)
                setattr(model, cv+"_modified", Expression(rule=temp_func))
        return model

    def load_init(self, model):
        self.base_pyomo_model.load_init(model)

    def set_base_point(self, model, base_point):
        for mv in self.modifier_mvs:
            getattr(model, mv + "_base").fix(base_point[mv])
