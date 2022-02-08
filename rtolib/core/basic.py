#!/usr/bin/env python
#-*- coding:utf-8 -*-

__name__ = ['Simulator', 'Optimizer', 'ParameterEstimator', 'ProblemDescription']

class Simulator():
    def simulate(self, input_values, param_values=None):
        return NotImplementedError()


class Optimizer():
    def optimize(self, param_values=None):
        return NotImplementedError()


class ParameterEstimator():
    def estimate_parameter(self, data_sets):
        return NotImplementedError()


class ProblemDescription():
    def __init__(self, symbol_list, bounds, scaling_factors, default_values):
        '''
        An example:
        symbol_list = {
                'MV': ("u1","u2"),# input
                'CV': ("obj"),# output
                'OBJ': "obj",# objective function
                'SPEC': None,  # specification
            }
        bounds = {
                "obj":(None, None),
                "u1":(None, None),
                "u2":(None, None),
            }
        scaling_factors = {
                "obj":1,
                "u1":1,
                "u2":1,
            }
        default_values = {
                "obj":1,
                "u1":2,
                "u2":3,
            }
        :param symbol_list:
        :param bounds:
        :param scaling_factors:
        :param default_values:
        '''
        self.bounds = bounds
        if scaling_factors is None:
            scaling_factors = {}
            for var_name in self.bounds.keys():
                if self.bounds[var_name][0] is not None and self.bounds[var_name][1] is not None:
                    scaling_factors[var_name] = self.bounds[var_name][1] - self.bounds[var_name][0]
                else:
                    scaling_factors[var_name] = 1
        self.scaling_factors = scaling_factors
        self.symbol_list = symbol_list
        self.default_values = default_values

