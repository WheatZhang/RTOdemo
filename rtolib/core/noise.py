#!/usr/bin/env python
#-*- coding:utf-8 -*-
import numpy

__name__ = ['NoiseGenerator']

class NoiseGenerator():
    def __init__(self, noise_level=None):
        '''

        :param noise_level: A dict, stardard deviation
        '''
        self.noise_level = noise_level
        self.history_noise = {}

    def get_noise(self,rto_iter,test_point_no,variable_name):
        flag = False
        for k in self.history_noise.keys():
            if k[0] == rto_iter and k[1] == test_point_no:
                flag = True
                break
        if flag:
            return self.history_noise[(rto_iter, test_point_no, variable_name)]
        else:
            for k,v in self.noise_level.items():
                noise_value = numpy.random.randn()*v
                if noise_value > v*2:
                    noise_value = v*2
                elif noise_value < -v*2:
                    noise_value = -2*v
                self.history_noise[(rto_iter,test_point_no,k)] = noise_value
            return self.history_noise[(rto_iter, test_point_no, variable_name)]

    def save_noise(self, filename):
        with open(filename, 'w') as fp:
            fp.write("rto_iter\ttest_point_no\tvariable_name\tnoise\n")
            for key, value in self.history_noise.items():
                fp.write("%d\t"%key[0])
                fp.write("%d\t" % key[1])
                fp.write("%s\t"%key[2])
                fp.write("%.6e\n"%value)

    def load_noise(self, filename):
        with open(filename, 'r') as f:
            line = f.readline()
            line = f.readline()
            while line:
                line_split = line.split(sep='\t')
                self.history_noise[(int(line_split[0]),int(line_split[1]),line_split[2])] = float(line_split[3])
                line = f.readline()