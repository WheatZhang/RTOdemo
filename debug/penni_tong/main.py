# coding=utf-8

from pyomo.environ import *
from pyomo.dae import *
from FBRPP_model import FBR


primitive_model = FBR()
process = ConcreteModel()
process.tf = Param(initialize=6)
process.t = ContinuousSet(bounds=(0, process.tf))
primitive_model.build(process, 'SIMU')
discretizer = TransformationFactory('dae.collocation')
discretizer.apply_to(process, nfe=1, ncp=3, scheme='LAGRANGE-RADAU')
solver = SolverFactory('ipopt', executable=r"../../external/bin/ipopt.exe")
primitive_model.initialize(process)
for i in process.t:
    process.F[i].fixed = True
    if i == process.t.first():
        process.S[i].fixed = True
        process.X[i].fixed = True
        process.P[i].fixed = True
        process.V[i].fixed = True
def update_initial(instance):
    for v in instance.component_objects(Var, active=True):
        for index in v:
            v[index] = value(v[instance.tf])
observations = ['X', 'P', 'S', 'V','F']
measurements = {}
for obs in observations:
    measurements[obs] = []
for simulation in range(32):
    solver.solve(process, tee=False)
    for obs in observations:
        for time in process.t:
            if time == 0:
                continue
            measurements[obs].append(value(getattr(process, obs)[time]))
    update_initial(process)
# for i in process.t:
#     if i == process.t.first():
#         process.F[i].fix(0.11)
#         process.S[i].fix(6.0)
#         process.X[i].fix(0.1)
#         process.P[i].fix(0.0)
#         process.V[i].fix(100.0)
#     process.F[i].fix(0.11)
# process.pprint()
for obs in observations:
    print(obs + ': ', measurements[obs])

