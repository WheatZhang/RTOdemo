#!/usr/bin/env python
#-*- coding:utf-8 -*-
from pyomo.environ import *
from pyomo.dae import *
import rtolib.util.init_value as ivt
import pyomo
from bisect import bisect_right, bisect_left
import copy


