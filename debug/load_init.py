#!/usr/bin/env python
#-*- coding:utf-8 -*-
import rtolib.model.HPC.PlantHPC as PlantHPC
import rtolib.util.init_value as init_value
from pyomo.environ import *

rto_plant = PlantHPC.RTO_Plant_HPC()
plant_simulator = ConcreteModel()
rto_plant.build(plant_simulator)
rto_plant.get_prosim_simulator(plant_simulator)
init_value.load_init_from_template(plant_simulator, r"F:\Research\RTOdemo\rtolib\model\HPC\init\NominalOpt_InitValue.txt",\
                                                       ignore_init_mismatch=True)
starting_point = {
        'TA_Flow':2000*0.14,
        'MA_Flow':2000*0.86,
    }
specifications = {
                'GAN_Flow': 700,
                'LN_Flow': 100,
                            }
solver = SolverFactory('ipopt', executable=r"../external/bin/ipopt.exe")
rto_plant.promoted_simulate(plant_simulator, starting_point, specifications, solver, tee=True)

init_value.to_template(plant_simulator, "perfect_init.txt")

plant_simulator2 = ConcreteModel()
rto_plant.build(plant_simulator2)
init_value.load_init_from_template(plant_simulator2, "perfect_init.txt",ignore_init_mismatch=True)
solver.solve(plant_simulator2, tee=True)