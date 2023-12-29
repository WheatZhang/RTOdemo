from rtolib.core.pyomo_model import PyomoModel
from pyomo.environ import *
import os

class ASU_Quadratic_Model(PyomoModel):
    def __init__(self):
        self.output_variables = {
            "cost":(lambda m: m.OBJ),
            "OxygenPrdtPurityCon":(lambda m: m.PurityCon1),
            "OxygenPrdtImpurityCon":(lambda m: m.PurityCon2),
            "NitrogenPrdtPurityCon":(lambda m: m.PurityCon3),
            "NitrogenPrdtImpurityCon":(lambda m: m.PurityCon4),
            "ArgonPrdtImpurityCon":(lambda m: m.PurityCon5),
            "CrudeArgonImpurityCon":(lambda m: m.PurityCon6),
            "LOX_Con":(lambda m: m.ProductCon1),
            "GAN_Con":(lambda m: m.ProductCon2),
            "GAR_Con":(lambda m: m.ProductCon3),
            "MainCoolerTempDiffCon":(lambda m: m.TempDiffCon1),
            "HeatLPCTempDiffCon":(lambda m: m.TempDiffCon2),
            "HPA_GOX_LB_Con":(lambda m: m.TurbineCon1),
            "HPA_GOX_UB_Con":(lambda m: m.TurbineCon2),
            "Feed_LB_Con":(lambda m: m.TurbineCon3),
            "Feed_UB_Con":(lambda m: m.TurbineCon4),
            "HPA_LB_Con":(lambda m: m.TurbineCon7),
            "HPA_UB_Con":(lambda m: m.TurbineCon8),
            "TA_LB_Con":(lambda m: m.TurbineCon9),
            "TA_UB_Con":(lambda m: m.TurbineCon10),
            "OverallHeatCon":(lambda m: m.OverallHeatBlnc),
        }
        self.noised_outputs = {
        }
        self.input_variables = {
            'FeedMFlow':(lambda m: m.FeedMFlow),
            'TAFeedFrac':(lambda m: m.FeedSplitterOut1Ratio), #TA Frac
            'MAFeedFrac':(lambda m: m.FeedSplitterOut2Ratio), #MA_Frac
            '40_LIN':(lambda m: m.HPCCondPrdtMFlow), #40-LIN
            '52_WN':(lambda m: m.LPCExt46MFlow), #52-WN
            '47_ARC':(lambda m: m.LPCExt15MFlow), #47-ARC
            'ASC_RefluxRatio':(lambda m: m.ASCCondRefRatio), #ASC_Reflux
            'LOX_Frac':(lambda m: m.OxSplitterOutARatio), #LOX_Frac
        }
        self.parameters = {
        }
        self.default_value={}
        self.parameter_scaling_factors={}
        self.output_noise = {
        }
        self.initial_value_file = os.path.join(os.path.dirname(__file__) + r"\GOX17000SimInitValue_model.txt")

    def build_body(self, test_model):
        test_model.FeedMFlow=Var(initialize=4100, bounds=(3800,4200))
        test_model.FeedSplitterOut1Ratio=Var(initialize=0.19, bounds=(0.15,0.25))
        test_model.FeedSplitterOut2Ratio=Var(initialize=0.5, bounds=(0.45,0.6))
        test_model.HPCCondPrdtMFlow=Var(initialize=1450, bounds=(1200,1500))
        test_model.LPCExt46MFlow=Var(initialize=1800, bounds=(1600,1900))
        test_model.LPCExt15MFlow=Var(initialize=600, bounds=(600,700))
        test_model.ASCCondRefRatio=Var(initialize=28, bounds=(25,30))
        test_model.OxSplitterOutARatio=Var(initialize=0.09, bounds=(0.07,0.1))
        test_model.FeedMFlow.fixed = True
        test_model.FeedSplitterOut1Ratio.fixed = True
        test_model.FeedSplitterOut2Ratio.fixed = True
        test_model.HPCCondPrdtMFlow.fixed = True
        test_model.LPCExt46MFlow.fixed = True
        test_model.LPCExt15MFlow.fixed = True
        test_model.ASCCondRefRatio.fixed = True
        test_model.OxSplitterOutARatio.fixed = True

        test_model.obj = Var(initialize=0)
        test_model.e = Var(initialize=0)

        def eq1(m):
            return m.e == 0.0*( \
            m.FeedMFlow+\
            m.FeedSplitterOut1Ratio+\
            m.FeedSplitterOut2Ratio+\
            m.HPCCondPrdtMFlow+\
            m.LPCExt46MFlow+\
            m.LPCExt15MFlow+\
            m.ASCCondRefRatio+\
            m.OxSplitterOutARatio
                    )
        test_model.eq1 = Constraint(rule=eq1)

        def _obj(m):
            return 0  #10*(m.Fb**2+m.Tr**2)
        test_model._obj = Expression(rule=_obj)

    def build_rto(self, test_model, cv_func):
        def PurityCon1(m):
            return m._obj+4.136111e-01 \
                + 7.846681e-06 * m.FeedMFlow \
                + -8.714713e-01 * m.FeedSplitterOut1Ratio \
                + -1.016297e+00 * m.FeedSplitterOut2Ratio \
                + -3.666371e-05 * m.HPCCondPrdtMFlow \
                + -5.269763e-06 * m.LPCExt46MFlow \
                + -4.918646e-05 * m.LPCExt15MFlow \
                + 1.931166e-09 * m.FeedMFlow * m.FeedMFlow \
                + -6.583403e-06 * m.FeedMFlow * m.FeedSplitterOut1Ratio \
                + -1.083568e-05 * m.FeedMFlow * m.FeedSplitterOut2Ratio \
                + -3.299821e-10 * m.FeedMFlow * m.HPCCondPrdtMFlow \
                + -1.258876e-09 * m.FeedMFlow * m.LPCExt46MFlow \
                + -1.836294e-09 * m.FeedMFlow * m.LPCExt15MFlow \
                + -6.583403e-06 * m.FeedSplitterOut1Ratio * m.FeedMFlow \
                + 4.853459e-01 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut1Ratio \
                + 5.604734e-01 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut2Ratio \
                + 1.676384e-05 * m.FeedSplitterOut1Ratio * m.HPCCondPrdtMFlow \
                + 5.460984e-06 * m.FeedSplitterOut1Ratio * m.LPCExt46MFlow \
                + 2.047440e-05 * m.FeedSplitterOut1Ratio * m.LPCExt15MFlow \
                + -1.083568e-05 * m.FeedSplitterOut2Ratio * m.FeedMFlow \
                + 5.604734e-01 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut1Ratio \
                + 6.890774e-01 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut2Ratio \
                + 2.502339e-05 * m.FeedSplitterOut2Ratio * m.HPCCondPrdtMFlow \
                + 1.131900e-05 * m.FeedSplitterOut2Ratio * m.LPCExt46MFlow \
                + 1.173508e-05 * m.FeedSplitterOut2Ratio * m.LPCExt15MFlow \
                + -3.299821e-10 * m.HPCCondPrdtMFlow * m.FeedMFlow \
                + 1.676384e-05 * m.HPCCondPrdtMFlow * m.FeedSplitterOut1Ratio \
                + 2.502339e-05 * m.HPCCondPrdtMFlow * m.FeedSplitterOut2Ratio \
                + 2.127717e-09 * m.HPCCondPrdtMFlow * m.HPCCondPrdtMFlow \
                + 9.045169e-10 * m.HPCCondPrdtMFlow * m.LPCExt46MFlow \
                + 1.977965e-09 * m.HPCCondPrdtMFlow * m.LPCExt15MFlow \
                + -1.258876e-09 * m.LPCExt46MFlow * m.FeedMFlow \
                + 5.460984e-06 * m.LPCExt46MFlow * m.FeedSplitterOut1Ratio \
                + 1.131900e-05 * m.LPCExt46MFlow * m.FeedSplitterOut2Ratio \
                + 9.045169e-10 * m.LPCExt46MFlow * m.HPCCondPrdtMFlow \
                + 1.568911e-09 * m.LPCExt46MFlow * m.LPCExt46MFlow \
                + 2.268358e-11 * m.LPCExt46MFlow * m.LPCExt15MFlow \
                + -1.836294e-09 * m.LPCExt15MFlow * m.FeedMFlow \
                + 2.047440e-05 * m.LPCExt15MFlow * m.FeedSplitterOut1Ratio \
                + 1.173508e-05 * m.LPCExt15MFlow * m.FeedSplitterOut2Ratio \
                + 1.977965e-09 * m.LPCExt15MFlow * m.HPCCondPrdtMFlow \
                + 2.268358e-11 * m.LPCExt15MFlow * m.LPCExt46MFlow \
                + 2.685155e-08 * m.LPCExt15MFlow * m.LPCExt15MFlow

        test_model.PurityCon1 = Expression(rule=PurityCon1)

        def PurityCon2(m):
            return -1.708457e-05 \
                + 1.096767e-09 * m.FeedMFlow \
                + -8.840838e-06 * m.FeedSplitterOut1Ratio \
                + -8.281442e-06 * m.FeedSplitterOut2Ratio \
                + -1.244825e-09 * m.HPCCondPrdtMFlow \
                + 1.322817e-09 * m.LPCExt46MFlow \
                + 1.770678e-09 * m.LPCExt15MFlow \
                + 1.005781e-13 * m.FeedMFlow * m.FeedMFlow \
                + 1.059882e-10 * m.FeedMFlow * m.FeedSplitterOut1Ratio \
                + 3.053729e-10 * m.FeedMFlow * m.FeedSplitterOut2Ratio \
                + 2.461017e-13 * m.FeedMFlow * m.HPCCondPrdtMFlow \
                + -4.772098e-14 * m.FeedMFlow * m.LPCExt46MFlow \
                + -5.127895e-14 * m.FeedMFlow * m.LPCExt15MFlow \
                + 1.059882e-10 * m.FeedSplitterOut1Ratio * m.FeedMFlow \
                + 1.672470e-05 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut1Ratio \
                + 3.286501e-06 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut2Ratio \
                + 6.916296e-10 * m.FeedSplitterOut1Ratio * m.HPCCondPrdtMFlow \
                + -6.629728e-10 * m.FeedSplitterOut1Ratio * m.LPCExt46MFlow \
                + -5.804718e-10 * m.FeedSplitterOut1Ratio * m.LPCExt15MFlow \
                + 3.053729e-10 * m.FeedSplitterOut2Ratio * m.FeedMFlow \
                + 3.286501e-06 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut1Ratio \
                + 9.237251e-06 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut2Ratio \
                + 3.738333e-10 * m.FeedSplitterOut2Ratio * m.HPCCondPrdtMFlow \
                + -1.454945e-10 * m.FeedSplitterOut2Ratio * m.LPCExt46MFlow \
                + 6.692029e-10 * m.FeedSplitterOut2Ratio * m.LPCExt15MFlow \
                + 2.461017e-13 * m.HPCCondPrdtMFlow * m.FeedMFlow \
                + 6.916296e-10 * m.HPCCondPrdtMFlow * m.FeedSplitterOut1Ratio \
                + 3.738333e-10 * m.HPCCondPrdtMFlow * m.FeedSplitterOut2Ratio \
                + 8.130929e-13 * m.HPCCondPrdtMFlow * m.HPCCondPrdtMFlow \
                + -9.516807e-14 * m.HPCCondPrdtMFlow * m.LPCExt46MFlow \
                + -8.179090e-14 * m.HPCCondPrdtMFlow * m.LPCExt15MFlow \
                + -4.772098e-14 * m.LPCExt46MFlow * m.FeedMFlow \
                + -6.629728e-10 * m.LPCExt46MFlow * m.FeedSplitterOut1Ratio \
                + -1.454945e-10 * m.LPCExt46MFlow * m.FeedSplitterOut2Ratio \
                + -9.516807e-14 * m.LPCExt46MFlow * m.HPCCondPrdtMFlow \
                + 1.909061e-13 * m.LPCExt46MFlow * m.LPCExt46MFlow \
                + -4.919431e-14 * m.LPCExt46MFlow * m.LPCExt15MFlow \
                + -5.127895e-14 * m.LPCExt15MFlow * m.FeedMFlow \
                + -5.804718e-10 * m.LPCExt15MFlow * m.FeedSplitterOut1Ratio \
                + 6.692029e-10 * m.LPCExt15MFlow * m.FeedSplitterOut2Ratio \
                + -8.179090e-14 * m.LPCExt15MFlow * m.HPCCondPrdtMFlow \
                + -4.919431e-14 * m.LPCExt15MFlow * m.LPCExt46MFlow \
                + 8.332001e-13 * m.LPCExt15MFlow * m.LPCExt15MFlow

        test_model.PurityCon2 = Expression(rule=PurityCon2)

        def PurityCon3(m):
            return 2.242935e-02 \
                + -1.053913e-05 * m.FeedMFlow \
                + -5.736427e-02 * m.FeedSplitterOut1Ratio \
                + -5.178457e-02 * m.FeedSplitterOut2Ratio \
                + 2.026795e-05 * m.HPCCondPrdtMFlow \
                + 6.116006e-06 * m.LPCExt46MFlow \
                + -6.505817e-06 * m.LPCExt15MFlow \
                + 1.344177e-09 * m.FeedMFlow * m.FeedMFlow \
                + 7.023972e-06 * m.FeedMFlow * m.FeedSplitterOut1Ratio \
                + 7.303014e-06 * m.FeedMFlow * m.FeedSplitterOut2Ratio \
                + -2.684490e-09 * m.FeedMFlow * m.HPCCondPrdtMFlow \
                + -8.313123e-10 * m.FeedMFlow * m.LPCExt46MFlow \
                + -5.779891e-12 * m.FeedMFlow * m.LPCExt15MFlow \
                + 7.023972e-06 * m.FeedSplitterOut1Ratio * m.FeedMFlow \
                + 5.018358e-02 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut1Ratio \
                + 3.772446e-02 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut2Ratio \
                + -1.371671e-05 * m.FeedSplitterOut1Ratio * m.HPCCondPrdtMFlow \
                + -5.384219e-06 * m.FeedSplitterOut1Ratio * m.LPCExt46MFlow \
                + -1.835455e-06 * m.FeedSplitterOut1Ratio * m.LPCExt15MFlow \
                + 7.303014e-06 * m.FeedSplitterOut2Ratio * m.FeedMFlow \
                + 3.772446e-02 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut1Ratio \
                + 4.315120e-02 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut2Ratio \
                + -1.599738e-05 * m.FeedSplitterOut2Ratio * m.HPCCondPrdtMFlow \
                + -5.184172e-06 * m.FeedSplitterOut2Ratio * m.LPCExt46MFlow \
                + -2.070738e-06 * m.FeedSplitterOut2Ratio * m.LPCExt15MFlow \
                + -2.684490e-09 * m.HPCCondPrdtMFlow * m.FeedMFlow \
                + -1.371671e-05 * m.HPCCondPrdtMFlow * m.FeedSplitterOut1Ratio \
                + -1.599738e-05 * m.HPCCondPrdtMFlow * m.FeedSplitterOut2Ratio \
                + 6.724575e-09 * m.HPCCondPrdtMFlow * m.HPCCondPrdtMFlow \
                + 1.348080e-09 * m.HPCCondPrdtMFlow * m.LPCExt46MFlow \
                + 2.932483e-10 * m.HPCCondPrdtMFlow * m.LPCExt15MFlow \
                + -8.313123e-10 * m.LPCExt46MFlow * m.FeedMFlow \
                + -5.384219e-06 * m.LPCExt46MFlow * m.FeedSplitterOut1Ratio \
                + -5.184172e-06 * m.LPCExt46MFlow * m.FeedSplitterOut2Ratio \
                + 1.348080e-09 * m.LPCExt46MFlow * m.HPCCondPrdtMFlow \
                + 1.212001e-09 * m.LPCExt46MFlow * m.LPCExt46MFlow \
                + 6.767309e-10 * m.LPCExt46MFlow * m.LPCExt15MFlow \
                + -5.779891e-12 * m.LPCExt15MFlow * m.FeedMFlow \
                + -1.835455e-06 * m.LPCExt15MFlow * m.FeedSplitterOut1Ratio \
                + -2.070738e-06 * m.LPCExt15MFlow * m.FeedSplitterOut2Ratio \
                + 2.932483e-10 * m.LPCExt15MFlow * m.HPCCondPrdtMFlow \
                + 6.767309e-10 * m.LPCExt15MFlow * m.LPCExt46MFlow \
                + 5.915605e-09 * m.LPCExt15MFlow * m.LPCExt15MFlow

        test_model.PurityCon3 = Expression(rule=PurityCon3)

        def PurityCon4(m):
            return 1.168076e+00 \
                + -3.537967e-04 * m.FeedMFlow \
                + -6.840375e-02 * m.FeedSplitterOut1Ratio \
                + -2.563277e-01 * m.FeedSplitterOut2Ratio \
                + -1.973433e-04 * m.HPCCondPrdtMFlow \
                + -1.770412e-04 * m.LPCExt46MFlow \
                + -7.153776e-04 * m.LPCExt15MFlow \
                + 3.369746e-08 * m.FeedMFlow * m.FeedMFlow \
                + -1.577247e-05 * m.FeedMFlow * m.FeedSplitterOut1Ratio \
                + 1.064079e-05 * m.FeedMFlow * m.FeedSplitterOut2Ratio \
                + 1.475841e-08 * m.FeedMFlow * m.HPCCondPrdtMFlow \
                + 1.299420e-08 * m.FeedMFlow * m.LPCExt46MFlow \
                + 3.838369e-08 * m.FeedMFlow * m.LPCExt15MFlow \
                + -1.577247e-05 * m.FeedSplitterOut1Ratio * m.FeedMFlow \
                + 5.760030e-01 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut1Ratio \
                + 7.175613e-03 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut2Ratio \
                + -2.162220e-06 * m.FeedSplitterOut1Ratio * m.HPCCondPrdtMFlow \
                + -1.948524e-05 * m.FeedSplitterOut1Ratio * m.LPCExt46MFlow \
                + 3.560678e-05 * m.FeedSplitterOut1Ratio * m.LPCExt15MFlow \
                + 1.064079e-05 * m.FeedSplitterOut2Ratio * m.FeedMFlow \
                + 7.175613e-03 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut1Ratio \
                + 1.272442e-01 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut2Ratio \
                + 8.932618e-06 * m.FeedSplitterOut2Ratio * m.HPCCondPrdtMFlow \
                + 1.986469e-05 * m.FeedSplitterOut2Ratio * m.LPCExt46MFlow \
                + -1.179999e-05 * m.FeedSplitterOut2Ratio * m.LPCExt15MFlow \
                + 1.475841e-08 * m.HPCCondPrdtMFlow * m.FeedMFlow \
                + -2.162220e-06 * m.HPCCondPrdtMFlow * m.FeedSplitterOut1Ratio \
                + 8.932618e-06 * m.HPCCondPrdtMFlow * m.FeedSplitterOut2Ratio \
                + 1.816849e-08 * m.HPCCondPrdtMFlow * m.HPCCondPrdtMFlow \
                + 9.230130e-09 * m.HPCCondPrdtMFlow * m.LPCExt46MFlow \
                + 2.812008e-08 * m.HPCCondPrdtMFlow * m.LPCExt15MFlow \
                + 1.299420e-08 * m.LPCExt46MFlow * m.FeedMFlow \
                + -1.948524e-05 * m.LPCExt46MFlow * m.FeedSplitterOut1Ratio \
                + 1.986469e-05 * m.LPCExt46MFlow * m.FeedSplitterOut2Ratio \
                + 9.230130e-09 * m.LPCExt46MFlow * m.HPCCondPrdtMFlow \
                + 9.565036e-09 * m.LPCExt46MFlow * m.LPCExt46MFlow \
                + 1.513800e-08 * m.LPCExt46MFlow * m.LPCExt15MFlow \
                + 3.838369e-08 * m.LPCExt15MFlow * m.FeedMFlow \
                + 3.560678e-05 * m.LPCExt15MFlow * m.FeedSplitterOut1Ratio \
                + -1.179999e-05 * m.LPCExt15MFlow * m.FeedSplitterOut2Ratio \
                + 2.812008e-08 * m.LPCExt15MFlow * m.HPCCondPrdtMFlow \
                + 1.513800e-08 * m.LPCExt15MFlow * m.LPCExt46MFlow \
                + 2.579645e-07 * m.LPCExt15MFlow * m.LPCExt15MFlow

        test_model.PurityCon4 = Expression(rule=PurityCon4)

        def PurityCon5(m):
            return 4.727750e+02 \
                + -8.887676e-02 * m.FeedMFlow \
                + -1.276279e+01 * m.FeedSplitterOut1Ratio \
                + -7.778652e+01 * m.FeedSplitterOut2Ratio \
                + -2.517770e-02 * m.HPCCondPrdtMFlow \
                + -3.557531e-02 * m.LPCExt46MFlow \
                + -3.344375e-01 * m.LPCExt15MFlow \
                + 8.879974e-06 * m.FeedMFlow * m.FeedMFlow \
                + 8.712841e-03 * m.FeedMFlow * m.FeedSplitterOut1Ratio \
                + 6.405738e-03 * m.FeedMFlow * m.FeedSplitterOut2Ratio \
                + -2.870353e-06 * m.FeedMFlow * m.HPCCondPrdtMFlow \
                + 3.862978e-06 * m.FeedMFlow * m.LPCExt46MFlow \
                + 3.874812e-06 * m.FeedMFlow * m.LPCExt15MFlow \
                + 8.712841e-03 * m.FeedSplitterOut1Ratio * m.FeedMFlow \
                + 3.753299e+02 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut1Ratio \
                + -1.702099e+02 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut2Ratio \
                + -1.630068e-02 * m.FeedSplitterOut1Ratio * m.HPCCondPrdtMFlow \
                + -1.061616e-02 * m.FeedSplitterOut1Ratio * m.LPCExt46MFlow \
                + -2.925291e-02 * m.FeedSplitterOut1Ratio * m.LPCExt15MFlow \
                + 6.405738e-03 * m.FeedSplitterOut2Ratio * m.FeedMFlow \
                + -1.702099e+02 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut1Ratio \
                + 1.107001e+02 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut2Ratio \
                + -1.918381e-02 * m.FeedSplitterOut2Ratio * m.HPCCondPrdtMFlow \
                + 7.626847e-03 * m.FeedSplitterOut2Ratio * m.LPCExt46MFlow \
                + 5.573367e-03 * m.FeedSplitterOut2Ratio * m.LPCExt15MFlow \
                + -2.870353e-06 * m.HPCCondPrdtMFlow * m.FeedMFlow \
                + -1.630068e-02 * m.HPCCondPrdtMFlow * m.FeedSplitterOut1Ratio \
                + -1.918381e-02 * m.HPCCondPrdtMFlow * m.FeedSplitterOut2Ratio \
                + 4.100577e-05 * m.HPCCondPrdtMFlow * m.HPCCondPrdtMFlow \
                + 6.186580e-06 * m.HPCCondPrdtMFlow * m.LPCExt46MFlow \
                + -4.027500e-06 * m.HPCCondPrdtMFlow * m.LPCExt15MFlow \
                + 3.862978e-06 * m.LPCExt46MFlow * m.FeedMFlow \
                + -1.061616e-02 * m.LPCExt46MFlow * m.FeedSplitterOut1Ratio \
                + 7.626847e-03 * m.LPCExt46MFlow * m.FeedSplitterOut2Ratio \
                + 6.186580e-06 * m.LPCExt46MFlow * m.HPCCondPrdtMFlow \
                + 5.040519e-06 * m.LPCExt46MFlow * m.LPCExt46MFlow \
                + -1.799556e-05 * m.LPCExt46MFlow * m.LPCExt15MFlow \
                + 3.874812e-06 * m.LPCExt15MFlow * m.FeedMFlow \
                + -2.925291e-02 * m.LPCExt15MFlow * m.FeedSplitterOut1Ratio \
                + 5.573367e-03 * m.LPCExt15MFlow * m.FeedSplitterOut2Ratio \
                + -4.027500e-06 * m.LPCExt15MFlow * m.HPCCondPrdtMFlow \
                + -1.799556e-05 * m.LPCExt15MFlow * m.LPCExt46MFlow \
                + 2.894634e-04 * m.LPCExt15MFlow * m.LPCExt15MFlow

        test_model.PurityCon5 = Expression(rule=PurityCon5)

        def PurityCon6(m):
            return 7.597834e-01 \
                + 2.241797e-05 * m.FeedMFlow \
                + -1.473651e+00 * m.FeedSplitterOut1Ratio \
                + -1.726238e+00 * m.FeedSplitterOut2Ratio \
                + -9.171998e-05 * m.HPCCondPrdtMFlow \
                + -4.111114e-05 * m.LPCExt46MFlow \
                + 5.723076e-05 * m.LPCExt15MFlow \
                + 1.114961e-08 * m.FeedMFlow * m.FeedMFlow \
                + -1.374683e-05 * m.FeedMFlow * m.FeedSplitterOut1Ratio \
                + -3.204405e-05 * m.FeedMFlow * m.FeedSplitterOut2Ratio \
                + -1.299118e-09 * m.FeedMFlow * m.HPCCondPrdtMFlow \
                + -6.446765e-09 * m.FeedMFlow * m.LPCExt46MFlow \
                + -2.952952e-08 * m.FeedMFlow * m.LPCExt15MFlow \
                + -1.374683e-05 * m.FeedSplitterOut1Ratio * m.FeedMFlow \
                + 7.242313e-01 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut1Ratio \
                + 8.106840e-01 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut2Ratio \
                + 4.790556e-05 * m.FeedSplitterOut1Ratio * m.HPCCondPrdtMFlow \
                + 2.782067e-05 * m.FeedSplitterOut1Ratio * m.LPCExt46MFlow \
                + 3.907479e-06 * m.FeedSplitterOut1Ratio * m.LPCExt15MFlow \
                + -3.204405e-05 * m.FeedSplitterOut2Ratio * m.FeedMFlow \
                + 8.106840e-01 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut1Ratio \
                + 1.180833e+00 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut2Ratio \
                + 6.898462e-05 * m.FeedSplitterOut2Ratio * m.HPCCondPrdtMFlow \
                + 4.694287e-05 * m.FeedSplitterOut2Ratio * m.LPCExt46MFlow \
                + 3.061085e-06 * m.FeedSplitterOut2Ratio * m.LPCExt15MFlow \
                + -1.299118e-09 * m.HPCCondPrdtMFlow * m.FeedMFlow \
                + 4.790556e-05 * m.HPCCondPrdtMFlow * m.FeedSplitterOut1Ratio \
                + 6.898462e-05 * m.HPCCondPrdtMFlow * m.FeedSplitterOut2Ratio \
                + 1.399655e-08 * m.HPCCondPrdtMFlow * m.HPCCondPrdtMFlow \
                + 6.449265e-09 * m.HPCCondPrdtMFlow * m.LPCExt46MFlow \
                + -1.268880e-08 * m.HPCCondPrdtMFlow * m.LPCExt15MFlow \
                + -6.446765e-09 * m.LPCExt46MFlow * m.FeedMFlow \
                + 2.782067e-05 * m.LPCExt46MFlow * m.FeedSplitterOut1Ratio \
                + 4.694287e-05 * m.LPCExt46MFlow * m.FeedSplitterOut2Ratio \
                + 6.449265e-09 * m.LPCExt46MFlow * m.HPCCondPrdtMFlow \
                + 2.729145e-08 * m.LPCExt46MFlow * m.LPCExt46MFlow \
                + -3.278529e-09 * m.LPCExt46MFlow * m.LPCExt15MFlow \
                + -2.952952e-08 * m.LPCExt15MFlow * m.FeedMFlow \
                + 3.907479e-06 * m.LPCExt15MFlow * m.FeedSplitterOut1Ratio \
                + 3.061085e-06 * m.LPCExt15MFlow * m.FeedSplitterOut2Ratio \
                + -1.268880e-08 * m.LPCExt15MFlow * m.HPCCondPrdtMFlow \
                + -3.278529e-09 * m.LPCExt15MFlow * m.LPCExt46MFlow \
                + 2.020865e-07 * m.LPCExt15MFlow * m.LPCExt15MFlow

        test_model.PurityCon6 = Expression(rule=PurityCon6)

        def ProductCon1(m):
            return 2.877770e+03 \
                + -6.836490e-01 * m.FeedMFlow \
                + 1.722631e+03 * m.FeedSplitterOut1Ratio \
                + 1.039570e+03 * m.FeedSplitterOut2Ratio \
                + 4.803167e-01 * m.HPCCondPrdtMFlow \
                + -8.162213e-01 * m.LPCExt46MFlow \
                + -2.567861e+00 * m.LPCExt15MFlow \
                + 4.111142e-05 * m.FeedMFlow * m.FeedMFlow \
                + 6.606313e-03 * m.FeedMFlow * m.FeedSplitterOut1Ratio \
                + 1.020992e-01 * m.FeedMFlow * m.FeedSplitterOut2Ratio \
                + -5.068371e-05 * m.FeedMFlow * m.HPCCondPrdtMFlow \
                + 8.969536e-06 * m.FeedMFlow * m.LPCExt46MFlow \
                + 4.537436e-05 * m.FeedMFlow * m.LPCExt15MFlow \
                + 6.606313e-03 * m.FeedSplitterOut1Ratio * m.FeedMFlow \
                + 1.529773e+03 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut1Ratio \
                + -8.187579e+02 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut2Ratio \
                + 6.953559e-02 * m.FeedSplitterOut1Ratio * m.HPCCondPrdtMFlow \
                + 1.977669e-02 * m.FeedSplitterOut1Ratio * m.LPCExt46MFlow \
                + 5.313162e-02 * m.FeedSplitterOut1Ratio * m.LPCExt15MFlow \
                + 1.020992e-01 * m.FeedSplitterOut2Ratio * m.FeedMFlow \
                + -8.187579e+02 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut1Ratio \
                + 1.149685e+03 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut2Ratio \
                + -9.356048e-02 * m.FeedSplitterOut2Ratio * m.HPCCondPrdtMFlow \
                + 1.539399e-01 * m.FeedSplitterOut2Ratio * m.LPCExt46MFlow \
                + -1.944238e-01 * m.FeedSplitterOut2Ratio * m.LPCExt15MFlow \
                + -5.068371e-05 * m.HPCCondPrdtMFlow * m.FeedMFlow \
                + 6.953559e-02 * m.HPCCondPrdtMFlow * m.FeedSplitterOut1Ratio \
                + -9.356048e-02 * m.HPCCondPrdtMFlow * m.FeedSplitterOut2Ratio \
                + 9.325048e-05 * m.HPCCondPrdtMFlow * m.HPCCondPrdtMFlow \
                + -6.793521e-06 * m.HPCCondPrdtMFlow * m.LPCExt46MFlow \
                + -5.686818e-05 * m.HPCCondPrdtMFlow * m.LPCExt15MFlow \
                + 8.969536e-06 * m.LPCExt46MFlow * m.FeedMFlow \
                + 1.977669e-02 * m.LPCExt46MFlow * m.FeedSplitterOut1Ratio \
                + 1.539399e-01 * m.LPCExt46MFlow * m.FeedSplitterOut2Ratio \
                + -6.793521e-06 * m.LPCExt46MFlow * m.HPCCondPrdtMFlow \
                + 1.638267e-04 * m.LPCExt46MFlow * m.LPCExt46MFlow \
                + 7.110248e-05 * m.LPCExt46MFlow * m.LPCExt15MFlow \
                + 4.537436e-05 * m.LPCExt15MFlow * m.FeedMFlow \
                + 5.313162e-02 * m.LPCExt15MFlow * m.FeedSplitterOut1Ratio \
                + -1.944238e-01 * m.LPCExt15MFlow * m.FeedSplitterOut2Ratio \
                + -5.686818e-05 * m.LPCExt15MFlow * m.HPCCondPrdtMFlow \
                + 7.110248e-05 * m.LPCExt15MFlow * m.LPCExt46MFlow \
                + 1.501238e-03 * m.LPCExt15MFlow * m.LPCExt15MFlow

        test_model.ProductCon1 = Expression(rule=ProductCon1)

        def ProductCon2(m):
            return 4.739451e+03 \
                + -1.613283e+00 * m.FeedMFlow \
                + -1.267970e+02 * m.FeedSplitterOut1Ratio \
                + -1.812799e+03 * m.FeedSplitterOut2Ratio \
                + 9.181696e-02 * m.HPCCondPrdtMFlow \
                + 9.093777e-01 * m.LPCExt46MFlow \
                + -4.419766e-01 * m.LPCExt15MFlow \
                + 1.463377e-04 * m.FeedMFlow * m.FeedMFlow \
                + -3.197490e-01 * m.FeedMFlow * m.FeedSplitterOut1Ratio \
                + -2.693997e-01 * m.FeedMFlow * m.FeedSplitterOut2Ratio \
                + -8.854321e-06 * m.FeedMFlow * m.HPCCondPrdtMFlow \
                + -4.023817e-07 * m.FeedMFlow * m.LPCExt46MFlow \
                + 1.037063e-05 * m.FeedMFlow * m.LPCExt15MFlow \
                + -3.197490e-01 * m.FeedSplitterOut1Ratio * m.FeedMFlow \
                + 1.362963e+03 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut1Ratio \
                + 2.077923e+02 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut2Ratio \
                + -7.084009e-02 * m.FeedSplitterOut1Ratio * m.HPCCondPrdtMFlow \
                + -8.906894e-02 * m.FeedSplitterOut1Ratio * m.LPCExt46MFlow \
                + -1.198111e-01 * m.FeedSplitterOut1Ratio * m.LPCExt15MFlow \
                + -2.693997e-01 * m.FeedSplitterOut2Ratio * m.FeedMFlow \
                + 2.077923e+02 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut1Ratio \
                + 7.918278e+02 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut2Ratio \
                + 6.428724e-02 * m.FeedSplitterOut2Ratio * m.HPCCondPrdtMFlow \
                + 4.776610e-02 * m.FeedSplitterOut2Ratio * m.LPCExt46MFlow \
                + -2.086383e-01 * m.FeedSplitterOut2Ratio * m.LPCExt15MFlow \
                + -8.854321e-06 * m.HPCCondPrdtMFlow * m.FeedMFlow \
                + -7.084009e-02 * m.HPCCondPrdtMFlow * m.FeedSplitterOut1Ratio \
                + 6.428724e-02 * m.HPCCondPrdtMFlow * m.FeedSplitterOut2Ratio \
                + 2.603571e-05 * m.HPCCondPrdtMFlow * m.HPCCondPrdtMFlow \
                + 1.934626e-05 * m.HPCCondPrdtMFlow * m.LPCExt46MFlow \
                + 5.939661e-06 * m.HPCCondPrdtMFlow * m.LPCExt15MFlow \
                + -4.023817e-07 * m.LPCExt46MFlow * m.FeedMFlow \
                + -8.906894e-02 * m.LPCExt46MFlow * m.FeedSplitterOut1Ratio \
                + 4.776610e-02 * m.LPCExt46MFlow * m.FeedSplitterOut2Ratio \
                + 1.934626e-05 * m.LPCExt46MFlow * m.HPCCondPrdtMFlow \
                + 1.781952e-05 * m.LPCExt46MFlow * m.LPCExt46MFlow \
                + 2.599063e-05 * m.LPCExt46MFlow * m.LPCExt15MFlow \
                + 1.037063e-05 * m.LPCExt15MFlow * m.FeedMFlow \
                + -1.198111e-01 * m.LPCExt15MFlow * m.FeedSplitterOut1Ratio \
                + -2.086383e-01 * m.LPCExt15MFlow * m.FeedSplitterOut2Ratio \
                + 5.939661e-06 * m.LPCExt15MFlow * m.HPCCondPrdtMFlow \
                + 2.599063e-05 * m.LPCExt15MFlow * m.LPCExt46MFlow \
                + 8.549584e-04 * m.LPCExt15MFlow * m.LPCExt15MFlow

        test_model.ProductCon2 = Expression(rule=ProductCon2)

        def ProductCon3(m):
            return 2.484087e+02 \
                + -1.480364e-02 * m.FeedMFlow \
                + 9.380115e+01 * m.FeedSplitterOut1Ratio \
                + -1.299723e+02 * m.FeedSplitterOut2Ratio \
                + 2.008139e-02 * m.HPCCondPrdtMFlow \
                + -1.790772e-02 * m.LPCExt46MFlow \
                + -1.722219e-01 * m.LPCExt15MFlow \
                + 2.087657e-06 * m.FeedMFlow * m.FeedMFlow \
                + -1.078872e-03 * m.FeedMFlow * m.FeedSplitterOut1Ratio \
                + 5.210133e-03 * m.FeedMFlow * m.FeedSplitterOut2Ratio \
                + -2.363780e-06 * m.FeedMFlow * m.HPCCondPrdtMFlow \
                + 4.107025e-07 * m.FeedMFlow * m.LPCExt46MFlow \
                + -1.179590e-06 * m.FeedMFlow * m.LPCExt15MFlow \
                + -1.078872e-03 * m.FeedSplitterOut1Ratio * m.FeedMFlow \
                + 3.034095e+02 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut1Ratio \
                + -1.528115e+02 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut2Ratio \
                + -1.607016e-02 * m.FeedSplitterOut1Ratio * m.HPCCondPrdtMFlow \
                + -3.311347e-03 * m.FeedSplitterOut1Ratio * m.LPCExt46MFlow \
                + -6.903744e-02 * m.FeedSplitterOut1Ratio * m.LPCExt15MFlow \
                + 5.210133e-03 * m.FeedSplitterOut2Ratio * m.FeedMFlow \
                + -1.528115e+02 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut1Ratio \
                + 1.178660e+02 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut2Ratio \
                + -1.344670e-02 * m.FeedSplitterOut2Ratio * m.HPCCondPrdtMFlow \
                + 1.136751e-02 * m.FeedSplitterOut2Ratio * m.LPCExt46MFlow \
                + 5.540295e-02 * m.FeedSplitterOut2Ratio * m.LPCExt15MFlow \
                + -2.363780e-06 * m.HPCCondPrdtMFlow * m.FeedMFlow \
                + -1.607016e-02 * m.HPCCondPrdtMFlow * m.FeedSplitterOut1Ratio \
                + -1.344670e-02 * m.HPCCondPrdtMFlow * m.FeedSplitterOut2Ratio \
                + 1.806180e-05 * m.HPCCondPrdtMFlow * m.HPCCondPrdtMFlow \
                + -8.931125e-07 * m.HPCCondPrdtMFlow * m.LPCExt46MFlow \
                + 2.495760e-06 * m.HPCCondPrdtMFlow * m.LPCExt15MFlow \
                + 4.107025e-07 * m.LPCExt46MFlow * m.FeedMFlow \
                + -3.311347e-03 * m.LPCExt46MFlow * m.FeedSplitterOut1Ratio \
                + 1.136751e-02 * m.LPCExt46MFlow * m.FeedSplitterOut2Ratio \
                + -8.931125e-07 * m.LPCExt46MFlow * m.HPCCondPrdtMFlow \
                + 8.762515e-06 * m.LPCExt46MFlow * m.LPCExt46MFlow \
                + 1.528807e-05 * m.LPCExt46MFlow * m.LPCExt15MFlow \
                + -1.179590e-06 * m.LPCExt15MFlow * m.FeedMFlow \
                + -6.903744e-02 * m.LPCExt15MFlow * m.FeedSplitterOut1Ratio \
                + 5.540295e-02 * m.LPCExt15MFlow * m.FeedSplitterOut2Ratio \
                + 2.495760e-06 * m.LPCExt15MFlow * m.HPCCondPrdtMFlow \
                + 1.528807e-05 * m.LPCExt15MFlow * m.LPCExt46MFlow \
                + 5.726421e-05 * m.LPCExt15MFlow * m.LPCExt15MFlow

        test_model.ProductCon3 = Expression(rule=ProductCon3)

        def TempDiffCon1(m):
            return 3.451987e+02 \
                + -6.796999e-02 * m.FeedMFlow \
                + -1.991949e+00 * m.FeedSplitterOut1Ratio \
                + 1.959784e+01 * m.FeedSplitterOut2Ratio \
                + 2.985118e-02 * m.HPCCondPrdtMFlow \
                + -4.575662e-02 * m.LPCExt46MFlow \
                + -2.702342e-01 * m.LPCExt15MFlow \
                + 9.424658e-06 * m.FeedMFlow * m.FeedMFlow \
                + -1.235292e-03 * m.FeedMFlow * m.FeedSplitterOut1Ratio \
                + -9.691227e-03 * m.FeedMFlow * m.FeedSplitterOut2Ratio \
                + -1.435008e-06 * m.FeedMFlow * m.HPCCondPrdtMFlow \
                + -2.124835e-06 * m.FeedMFlow * m.LPCExt46MFlow \
                + 1.429278e-05 * m.FeedMFlow * m.LPCExt15MFlow \
                + -1.235292e-03 * m.FeedSplitterOut1Ratio * m.FeedMFlow \
                + 1.388117e+02 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut1Ratio \
                + -7.606320e+01 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut2Ratio \
                + -1.626051e-03 * m.FeedSplitterOut1Ratio * m.HPCCondPrdtMFlow \
                + 7.206331e-03 * m.FeedSplitterOut1Ratio * m.LPCExt46MFlow \
                + -4.173991e-03 * m.FeedSplitterOut1Ratio * m.LPCExt15MFlow \
                + -9.691227e-03 * m.FeedSplitterOut2Ratio * m.FeedMFlow \
                + -7.606320e+01 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut1Ratio \
                + 7.846691e+01 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut2Ratio \
                + -1.290692e-02 * m.FeedSplitterOut2Ratio * m.HPCCondPrdtMFlow \
                + 8.527967e-03 * m.FeedSplitterOut2Ratio * m.LPCExt46MFlow \
                + 4.852268e-03 * m.FeedSplitterOut2Ratio * m.LPCExt15MFlow \
                + -1.435008e-06 * m.HPCCondPrdtMFlow * m.FeedMFlow \
                + -1.626051e-03 * m.HPCCondPrdtMFlow * m.FeedSplitterOut1Ratio \
                + -1.290692e-02 * m.HPCCondPrdtMFlow * m.FeedSplitterOut2Ratio \
                + 1.577950e-05 * m.HPCCondPrdtMFlow * m.HPCCondPrdtMFlow \
                + -2.721676e-06 * m.HPCCondPrdtMFlow * m.LPCExt46MFlow \
                + -1.191051e-05 * m.HPCCondPrdtMFlow * m.LPCExt15MFlow \
                + -2.124835e-06 * m.LPCExt46MFlow * m.FeedMFlow \
                + 7.206331e-03 * m.LPCExt46MFlow * m.FeedSplitterOut1Ratio \
                + 8.527967e-03 * m.LPCExt46MFlow * m.FeedSplitterOut2Ratio \
                + -2.721676e-06 * m.LPCExt46MFlow * m.HPCCondPrdtMFlow \
                + 1.514841e-05 * m.LPCExt46MFlow * m.LPCExt46MFlow \
                + 8.365387e-06 * m.LPCExt46MFlow * m.LPCExt15MFlow \
                + 1.429278e-05 * m.LPCExt15MFlow * m.FeedMFlow \
                + -4.173991e-03 * m.LPCExt15MFlow * m.FeedSplitterOut1Ratio \
                + 4.852268e-03 * m.LPCExt15MFlow * m.FeedSplitterOut2Ratio \
                + -1.191051e-05 * m.LPCExt15MFlow * m.HPCCondPrdtMFlow \
                + 8.365387e-06 * m.LPCExt15MFlow * m.LPCExt46MFlow \
                + 6.575398e-05 * m.LPCExt15MFlow * m.LPCExt15MFlow

        test_model.TempDiffCon1 = Expression(rule=TempDiffCon1)

        def TempDiffCon2(m):
            return 1.311966e+02 \
                + -5.138072e-03 * m.FeedMFlow \
                + -5.831702e+01 * m.FeedSplitterOut1Ratio \
                + -8.082675e+01 * m.FeedSplitterOut2Ratio \
                + -4.452770e-02 * m.HPCCondPrdtMFlow \
                + -2.507427e-02 * m.LPCExt46MFlow \
                + -4.592842e-02 * m.LPCExt15MFlow \
                + 4.924680e-06 * m.FeedMFlow * m.FeedMFlow \
                + -1.347210e-02 * m.FeedMFlow * m.FeedSplitterOut1Ratio \
                + -1.046151e-02 * m.FeedMFlow * m.FeedSplitterOut2Ratio \
                + -1.083377e-06 * m.FeedMFlow * m.HPCCondPrdtMFlow \
                + -1.327534e-06 * m.FeedMFlow * m.LPCExt46MFlow \
                + 8.384334e-06 * m.FeedMFlow * m.LPCExt15MFlow \
                + -1.347210e-02 * m.FeedSplitterOut1Ratio * m.FeedMFlow \
                + 1.330109e+02 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut1Ratio \
                + 2.297401e+01 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut2Ratio \
                + 1.103084e-02 * m.FeedSplitterOut1Ratio * m.HPCCondPrdtMFlow \
                + -7.007193e-03 * m.FeedSplitterOut1Ratio * m.LPCExt46MFlow \
                + -1.534999e-02 * m.FeedSplitterOut1Ratio * m.LPCExt15MFlow \
                + -1.046151e-02 * m.FeedSplitterOut2Ratio * m.FeedMFlow \
                + 2.297401e+01 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut1Ratio \
                + 4.974894e+01 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut2Ratio \
                + 1.124722e-02 * m.FeedSplitterOut2Ratio * m.HPCCondPrdtMFlow \
                + -5.762404e-03 * m.FeedSplitterOut2Ratio * m.LPCExt46MFlow \
                + -1.225246e-02 * m.FeedSplitterOut2Ratio * m.LPCExt15MFlow \
                + -1.083377e-06 * m.HPCCondPrdtMFlow * m.FeedMFlow \
                + 1.103084e-02 * m.HPCCondPrdtMFlow * m.FeedSplitterOut1Ratio \
                + 1.124722e-02 * m.HPCCondPrdtMFlow * m.FeedSplitterOut2Ratio \
                + 1.358046e-05 * m.HPCCondPrdtMFlow * m.HPCCondPrdtMFlow \
                + 3.623825e-06 * m.HPCCondPrdtMFlow * m.LPCExt46MFlow \
                + -1.448519e-05 * m.HPCCondPrdtMFlow * m.LPCExt15MFlow \
                + -1.327534e-06 * m.LPCExt46MFlow * m.FeedMFlow \
                + -7.007193e-03 * m.LPCExt46MFlow * m.FeedSplitterOut1Ratio \
                + -5.762404e-03 * m.LPCExt46MFlow * m.FeedSplitterOut2Ratio \
                + 3.623825e-06 * m.LPCExt46MFlow * m.HPCCondPrdtMFlow \
                + 1.073869e-05 * m.LPCExt46MFlow * m.LPCExt46MFlow \
                + -1.388330e-05 * m.LPCExt46MFlow * m.LPCExt15MFlow \
                + 8.384334e-06 * m.LPCExt15MFlow * m.FeedMFlow \
                + -1.534999e-02 * m.LPCExt15MFlow * m.FeedSplitterOut1Ratio \
                + -1.225246e-02 * m.LPCExt15MFlow * m.FeedSplitterOut2Ratio \
                + -1.448519e-05 * m.LPCExt15MFlow * m.HPCCondPrdtMFlow \
                + -1.388330e-05 * m.LPCExt15MFlow * m.LPCExt46MFlow \
                + 8.792963e-05 * m.LPCExt15MFlow * m.LPCExt15MFlow

        test_model.TempDiffCon2 = Expression(rule=TempDiffCon2)

        def TurbineCon1(m):
            return 1.929154e+03 \
                + -9.149507e-01 * m.FeedMFlow \
                + -2.365893e+01 * m.FeedSplitterOut1Ratio \
                + -2.337344e+02 * m.FeedSplitterOut2Ratio \
                + 3.405128e-02 * m.HPCCondPrdtMFlow \
                + -5.709730e-02 * m.LPCExt46MFlow \
                + -1.895487e-01 * m.LPCExt15MFlow \
                + 7.538829e-05 * m.FeedMFlow * m.FeedMFlow \
                + 1.423657e-01 * m.FeedMFlow * m.FeedSplitterOut1Ratio \
                + 9.080850e-02 * m.FeedMFlow * m.FeedSplitterOut2Ratio \
                + -9.962140e-06 * m.FeedMFlow * m.HPCCondPrdtMFlow \
                + 1.130049e-05 * m.FeedMFlow * m.LPCExt46MFlow \
                + -6.230546e-06 * m.FeedMFlow * m.LPCExt15MFlow \
                + 1.423657e-01 * m.FeedSplitterOut1Ratio * m.FeedMFlow \
                + 1.370110e+03 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut1Ratio \
                + -7.762174e+01 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut2Ratio \
                + -6.520978e-03 * m.FeedSplitterOut1Ratio * m.HPCCondPrdtMFlow \
                + -2.358483e-02 * m.FeedSplitterOut1Ratio * m.LPCExt46MFlow \
                + -3.584111e-01 * m.FeedSplitterOut1Ratio * m.LPCExt15MFlow \
                + 9.080850e-02 * m.FeedSplitterOut2Ratio * m.FeedMFlow \
                + -7.762174e+01 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut1Ratio \
                + 5.921761e+02 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut2Ratio \
                + 4.096751e-02 * m.FeedSplitterOut2Ratio * m.HPCCondPrdtMFlow \
                + 2.411873e-02 * m.FeedSplitterOut2Ratio * m.LPCExt46MFlow \
                + -8.931490e-02 * m.FeedSplitterOut2Ratio * m.LPCExt15MFlow \
                + -9.962140e-06 * m.HPCCondPrdtMFlow * m.FeedMFlow \
                + -6.520978e-03 * m.HPCCondPrdtMFlow * m.FeedSplitterOut1Ratio \
                + 4.096751e-02 * m.HPCCondPrdtMFlow * m.FeedSplitterOut2Ratio \
                + 8.615390e-05 * m.HPCCondPrdtMFlow * m.HPCCondPrdtMFlow \
                + -1.538989e-06 * m.HPCCondPrdtMFlow * m.LPCExt46MFlow \
                + 1.047859e-06 * m.HPCCondPrdtMFlow * m.LPCExt15MFlow \
                + 1.130049e-05 * m.LPCExt46MFlow * m.FeedMFlow \
                + -2.358483e-02 * m.LPCExt46MFlow * m.FeedSplitterOut1Ratio \
                + 2.411873e-02 * m.LPCExt46MFlow * m.FeedSplitterOut2Ratio \
                + -1.538989e-06 * m.LPCExt46MFlow * m.HPCCondPrdtMFlow \
                + 1.542954e-05 * m.LPCExt46MFlow * m.LPCExt46MFlow \
                + -3.426560e-05 * m.LPCExt46MFlow * m.LPCExt15MFlow \
                + -6.230546e-06 * m.LPCExt15MFlow * m.FeedMFlow \
                + -3.584111e-01 * m.LPCExt15MFlow * m.FeedSplitterOut1Ratio \
                + -8.931490e-02 * m.LPCExt15MFlow * m.FeedSplitterOut2Ratio \
                + 1.047859e-06 * m.LPCExt15MFlow * m.HPCCondPrdtMFlow \
                + -3.426560e-05 * m.LPCExt15MFlow * m.LPCExt46MFlow \
                + 3.931440e-04 * m.LPCExt15MFlow * m.LPCExt15MFlow

        test_model.TurbineCon1 = Expression(rule=TurbineCon1)

        def TurbineCon2(m):
            return 4.896889e+03 \
                + -6.615718e-01 * m.FeedMFlow \
                + -1.818634e+02 * m.FeedSplitterOut1Ratio \
                + -1.893035e+02 * m.FeedSplitterOut2Ratio \
                + 4.946921e-01 * m.HPCCondPrdtMFlow \
                + -8.514016e-01 * m.LPCExt46MFlow \
                + -4.943770e+00 * m.LPCExt15MFlow \
                + 9.150370e-05 * m.FeedMFlow * m.FeedMFlow \
                + -5.915217e-02 * m.FeedMFlow * m.FeedSplitterOut1Ratio \
                + -9.684638e-02 * m.FeedMFlow * m.FeedSplitterOut2Ratio \
                + -4.665040e-05 * m.FeedMFlow * m.HPCCondPrdtMFlow \
                + -1.388759e-05 * m.FeedMFlow * m.LPCExt46MFlow \
                + 2.823670e-04 * m.FeedMFlow * m.LPCExt15MFlow \
                + -5.915217e-02 * m.FeedSplitterOut1Ratio * m.FeedMFlow \
                + 2.566931e+03 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut1Ratio \
                + -1.087508e+03 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut2Ratio \
                + 4.492365e-03 * m.FeedSplitterOut1Ratio * m.HPCCondPrdtMFlow \
                + 7.558592e-02 * m.FeedSplitterOut1Ratio * m.LPCExt46MFlow \
                + -7.391571e-02 * m.FeedSplitterOut1Ratio * m.LPCExt15MFlow \
                + -9.684638e-02 * m.FeedSplitterOut2Ratio * m.FeedMFlow \
                + -1.087508e+03 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut1Ratio \
                + 1.487346e+03 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut2Ratio \
                + -1.743981e-01 * m.FeedSplitterOut2Ratio * m.HPCCondPrdtMFlow \
                + 2.307688e-01 * m.FeedSplitterOut2Ratio * m.LPCExt46MFlow \
                + -3.110198e-01 * m.FeedSplitterOut2Ratio * m.LPCExt15MFlow \
                + -4.665040e-05 * m.HPCCondPrdtMFlow * m.FeedMFlow \
                + 4.492365e-03 * m.HPCCondPrdtMFlow * m.FeedSplitterOut1Ratio \
                + -1.743981e-01 * m.HPCCondPrdtMFlow * m.FeedSplitterOut2Ratio \
                + 1.975548e-04 * m.HPCCondPrdtMFlow * m.HPCCondPrdtMFlow \
                + -3.647214e-05 * m.HPCCondPrdtMFlow * m.LPCExt46MFlow \
                + -8.615980e-05 * m.HPCCondPrdtMFlow * m.LPCExt15MFlow \
                + -1.388759e-05 * m.LPCExt46MFlow * m.FeedMFlow \
                + 7.558592e-02 * m.LPCExt46MFlow * m.FeedSplitterOut1Ratio \
                + 2.307688e-01 * m.LPCExt46MFlow * m.FeedSplitterOut2Ratio \
                + -3.647214e-05 * m.LPCExt46MFlow * m.HPCCondPrdtMFlow \
                + 2.436086e-04 * m.LPCExt46MFlow * m.LPCExt46MFlow \
                + 5.069895e-06 * m.LPCExt46MFlow * m.LPCExt15MFlow \
                + 2.823670e-04 * m.LPCExt15MFlow * m.FeedMFlow \
                + -7.391571e-02 * m.LPCExt15MFlow * m.FeedSplitterOut1Ratio \
                + -3.110198e-01 * m.LPCExt15MFlow * m.FeedSplitterOut2Ratio \
                + -8.615980e-05 * m.LPCExt15MFlow * m.HPCCondPrdtMFlow \
                + 5.069895e-06 * m.LPCExt15MFlow * m.LPCExt46MFlow \
                + 1.757581e-03 * m.LPCExt15MFlow * m.LPCExt15MFlow

        test_model.TurbineCon2 = Expression(rule=TurbineCon2)

        def TurbineCon3(m):
            return 2.676946e+03 \
                + -1.000001e+00 * m.FeedMFlow \
                + 1.724592e-02 * m.FeedSplitterOut1Ratio \
                + -1.492163e-03 * m.FeedSplitterOut2Ratio \
                + 3.132209e-07 * m.HPCCondPrdtMFlow \
                + -1.395711e-06 * m.LPCExt46MFlow \
                + -1.795876e-06 * m.LPCExt15MFlow \
                + 1.647452e-10 * m.FeedMFlow * m.FeedMFlow \
                + -2.032854e-06 * m.FeedMFlow * m.FeedSplitterOut1Ratio \
                + -3.878549e-07 * m.FeedMFlow * m.FeedSplitterOut2Ratio \
                + -7.685020e-11 * m.FeedMFlow * m.HPCCondPrdtMFlow \
                + -4.289785e-11 * m.FeedMFlow * m.LPCExt46MFlow \
                + 3.354718e-10 * m.FeedMFlow * m.LPCExt15MFlow \
                + -2.032854e-06 * m.FeedSplitterOut1Ratio * m.FeedMFlow \
                + 4.501317e-02 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut1Ratio \
                + -9.710598e-04 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut2Ratio \
                + 1.283487e-06 * m.FeedSplitterOut1Ratio * m.HPCCondPrdtMFlow \
                + -2.111883e-06 * m.FeedSplitterOut1Ratio * m.LPCExt46MFlow \
                + -5.551282e-06 * m.FeedSplitterOut1Ratio * m.LPCExt15MFlow \
                + -3.878549e-07 * m.FeedSplitterOut2Ratio * m.FeedMFlow \
                + -9.710598e-04 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut1Ratio \
                + 3.297396e-03 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut2Ratio \
                + 3.020174e-07 * m.FeedSplitterOut2Ratio * m.HPCCondPrdtMFlow \
                + 1.112866e-06 * m.FeedSplitterOut2Ratio * m.LPCExt46MFlow \
                + -1.484395e-07 * m.FeedSplitterOut2Ratio * m.LPCExt15MFlow \
                + -7.685020e-11 * m.HPCCondPrdtMFlow * m.FeedMFlow \
                + 1.283487e-06 * m.HPCCondPrdtMFlow * m.FeedSplitterOut1Ratio \
                + 3.020174e-07 * m.HPCCondPrdtMFlow * m.FeedSplitterOut2Ratio \
                + 1.309166e-10 * m.HPCCondPrdtMFlow * m.HPCCondPrdtMFlow \
                + 1.056200e-11 * m.HPCCondPrdtMFlow * m.LPCExt46MFlow \
                + 3.276746e-12 * m.HPCCondPrdtMFlow * m.LPCExt15MFlow \
                + -4.289785e-11 * m.LPCExt46MFlow * m.FeedMFlow \
                + -2.111883e-06 * m.LPCExt46MFlow * m.FeedSplitterOut1Ratio \
                + 1.112866e-06 * m.LPCExt46MFlow * m.FeedSplitterOut2Ratio \
                + 1.056200e-11 * m.LPCExt46MFlow * m.HPCCondPrdtMFlow \
                + 5.221009e-10 * m.LPCExt46MFlow * m.LPCExt46MFlow \
                + -1.853575e-11 * m.LPCExt46MFlow * m.LPCExt15MFlow \
                + 3.354718e-10 * m.LPCExt15MFlow * m.FeedMFlow \
                + -5.551282e-06 * m.LPCExt15MFlow * m.FeedSplitterOut1Ratio \
                + -1.484395e-07 * m.LPCExt15MFlow * m.FeedSplitterOut2Ratio \
                + 3.276746e-12 * m.LPCExt15MFlow * m.HPCCondPrdtMFlow \
                + -1.853575e-11 * m.LPCExt15MFlow * m.LPCExt46MFlow \
                + 1.413755e-09 * m.LPCExt15MFlow * m.LPCExt15MFlow

        test_model.TurbineCon3 = Expression(rule=TurbineCon3)

        def TurbineCon4(m):
            return -6.690816e+03 \
                + 9.997947e-01 * m.FeedMFlow \
                + -7.052821e-01 * m.FeedSplitterOut1Ratio \
                + -9.247743e-02 * m.FeedSplitterOut2Ratio \
                + -4.438051e-05 * m.HPCCondPrdtMFlow \
                + -3.485485e-04 * m.LPCExt46MFlow \
                + -1.426995e-03 * m.LPCExt15MFlow \
                + 3.072367e-08 * m.FeedMFlow * m.FeedMFlow \
                + 2.632164e-05 * m.FeedMFlow * m.FeedSplitterOut1Ratio \
                + -1.476965e-05 * m.FeedMFlow * m.FeedSplitterOut2Ratio \
                + -8.647477e-09 * m.FeedMFlow * m.HPCCondPrdtMFlow \
                + 7.990795e-09 * m.FeedMFlow * m.LPCExt46MFlow \
                + -4.808839e-09 * m.FeedMFlow * m.LPCExt15MFlow \
                + 2.632164e-05 * m.FeedSplitterOut1Ratio * m.FeedMFlow \
                + 1.173925e+00 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut1Ratio \
                + -1.581967e-01 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut2Ratio \
                + 1.658289e-05 * m.FeedSplitterOut1Ratio * m.HPCCondPrdtMFlow \
                + 1.543669e-05 * m.FeedSplitterOut1Ratio * m.LPCExt46MFlow \
                + -1.014905e-04 * m.FeedSplitterOut1Ratio * m.LPCExt15MFlow \
                + -1.476965e-05 * m.FeedSplitterOut2Ratio * m.FeedMFlow \
                + -1.581967e-01 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut1Ratio \
                + 2.251553e-01 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut2Ratio \
                + -2.887210e-05 * m.FeedSplitterOut2Ratio * m.HPCCondPrdtMFlow \
                + -2.030386e-05 * m.FeedSplitterOut2Ratio * m.LPCExt46MFlow \
                + 1.118116e-04 * m.FeedSplitterOut2Ratio * m.LPCExt15MFlow \
                + -8.647477e-09 * m.HPCCondPrdtMFlow * m.FeedMFlow \
                + 1.658289e-05 * m.HPCCondPrdtMFlow * m.FeedSplitterOut1Ratio \
                + -2.887210e-05 * m.HPCCondPrdtMFlow * m.FeedSplitterOut2Ratio \
                + 4.873651e-08 * m.HPCCondPrdtMFlow * m.HPCCondPrdtMFlow \
                + -2.953022e-09 * m.HPCCondPrdtMFlow * m.LPCExt46MFlow \
                + 6.548000e-08 * m.HPCCondPrdtMFlow * m.LPCExt15MFlow \
                + 7.990795e-09 * m.LPCExt46MFlow * m.FeedMFlow \
                + 1.543669e-05 * m.LPCExt46MFlow * m.FeedSplitterOut1Ratio \
                + -2.030386e-05 * m.LPCExt46MFlow * m.FeedSplitterOut2Ratio \
                + -2.953022e-09 * m.LPCExt46MFlow * m.HPCCondPrdtMFlow \
                + 5.219886e-08 * m.LPCExt46MFlow * m.LPCExt46MFlow \
                + 3.895181e-08 * m.LPCExt46MFlow * m.LPCExt15MFlow \
                + -4.808839e-09 * m.LPCExt15MFlow * m.FeedMFlow \
                + -1.014905e-04 * m.LPCExt15MFlow * m.FeedSplitterOut1Ratio \
                + 1.118116e-04 * m.LPCExt15MFlow * m.FeedSplitterOut2Ratio \
                + 6.548000e-08 * m.LPCExt15MFlow * m.HPCCondPrdtMFlow \
                + 3.895181e-08 * m.LPCExt15MFlow * m.LPCExt46MFlow \
                + 1.029680e-06 * m.LPCExt15MFlow * m.LPCExt15MFlow

        test_model.TurbineCon4 = Expression(rule=TurbineCon4)

        def TurbineCon7(m):
            return 6.491782e+03 \
                + -2.627825e+00 * m.FeedMFlow \
                + 2.609630e+02 * m.FeedSplitterOut1Ratio \
                + -1.152923e+03 * m.FeedSplitterOut2Ratio \
                + -1.679721e-01 * m.HPCCondPrdtMFlow \
                + -3.412584e-01 * m.LPCExt46MFlow \
                + -2.248197e+00 * m.LPCExt15MFlow \
                + 1.921956e-04 * m.FeedMFlow * m.FeedMFlow \
                + 3.608531e-01 * m.FeedMFlow * m.FeedSplitterOut1Ratio \
                + 4.290382e-01 * m.FeedMFlow * m.FeedSplitterOut2Ratio \
                + -5.041510e-07 * m.FeedMFlow * m.HPCCondPrdtMFlow \
                + 2.104403e-05 * m.FeedMFlow * m.LPCExt46MFlow \
                + 8.682139e-05 * m.FeedMFlow * m.LPCExt15MFlow \
                + 3.608531e-01 * m.FeedSplitterOut1Ratio * m.FeedMFlow \
                + 2.180077e+03 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut1Ratio \
                + -5.804482e+01 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut2Ratio \
                + -3.585611e-03 * m.FeedSplitterOut1Ratio * m.HPCCondPrdtMFlow \
                + -9.963875e-02 * m.FeedSplitterOut1Ratio * m.LPCExt46MFlow \
                + 2.841542e-01 * m.FeedSplitterOut1Ratio * m.LPCExt15MFlow \
                + 4.290382e-01 * m.FeedSplitterOut2Ratio * m.FeedMFlow \
                + -5.804482e+01 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut1Ratio \
                + 1.510505e+03 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut2Ratio \
                + 1.402917e-02 * m.FeedSplitterOut2Ratio * m.HPCCondPrdtMFlow \
                + 1.323196e-01 * m.FeedSplitterOut2Ratio * m.LPCExt46MFlow \
                + -8.954633e-03 * m.FeedSplitterOut2Ratio * m.LPCExt15MFlow \
                + -5.041510e-07 * m.HPCCondPrdtMFlow * m.FeedMFlow \
                + -3.585611e-03 * m.HPCCondPrdtMFlow * m.FeedSplitterOut1Ratio \
                + 1.402917e-02 * m.HPCCondPrdtMFlow * m.FeedSplitterOut2Ratio \
                + 1.118795e-04 * m.HPCCondPrdtMFlow * m.HPCCondPrdtMFlow \
                + -1.277515e-05 * m.HPCCondPrdtMFlow * m.LPCExt46MFlow \
                + 1.181741e-04 * m.HPCCondPrdtMFlow * m.LPCExt15MFlow \
                + 2.104403e-05 * m.LPCExt46MFlow * m.FeedMFlow \
                + -9.963875e-02 * m.LPCExt46MFlow * m.FeedSplitterOut1Ratio \
                + 1.323196e-01 * m.LPCExt46MFlow * m.FeedSplitterOut2Ratio \
                + -1.277515e-05 * m.LPCExt46MFlow * m.HPCCondPrdtMFlow \
                + 5.076646e-05 * m.LPCExt46MFlow * m.LPCExt46MFlow \
                + -1.347669e-04 * m.LPCExt46MFlow * m.LPCExt15MFlow \
                + 8.682139e-05 * m.LPCExt15MFlow * m.FeedMFlow \
                + 2.841542e-01 * m.LPCExt15MFlow * m.FeedSplitterOut1Ratio \
                + -8.954633e-03 * m.LPCExt15MFlow * m.FeedSplitterOut2Ratio \
                + 1.181741e-04 * m.LPCExt15MFlow * m.HPCCondPrdtMFlow \
                + -1.347669e-04 * m.LPCExt15MFlow * m.LPCExt46MFlow \
                + 8.790847e-04 * m.LPCExt15MFlow * m.LPCExt15MFlow

        test_model.TurbineCon7 = Expression(rule=TurbineCon7)

        def TurbineCon8(m):
            return 3.851954e+03 \
                + -1.407389e+00 * m.FeedMFlow \
                + -1.117096e+02 * m.FeedSplitterOut1Ratio \
                + -2.325308e+03 * m.FeedSplitterOut2Ratio \
                + -2.062197e-02 * m.HPCCondPrdtMFlow \
                + -1.444900e-01 * m.LPCExt46MFlow \
                + 4.758704e-02 * m.LPCExt15MFlow \
                + 2.795484e-04 * m.FeedMFlow * m.FeedMFlow \
                + -4.880579e-01 * m.FeedMFlow * m.FeedSplitterOut1Ratio \
                + -4.110307e-01 * m.FeedMFlow * m.FeedSplitterOut2Ratio \
                + 1.825708e-06 * m.FeedMFlow * m.HPCCondPrdtMFlow \
                + 2.471715e-06 * m.FeedMFlow * m.LPCExt46MFlow \
                + 3.139748e-05 * m.FeedMFlow * m.LPCExt15MFlow \
                + -4.880579e-01 * m.FeedSplitterOut1Ratio * m.FeedMFlow \
                + 1.821923e+03 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut1Ratio \
                + -6.429820e+01 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut2Ratio \
                + -7.710441e-02 * m.FeedSplitterOut1Ratio * m.HPCCondPrdtMFlow \
                + -1.328762e-01 * m.FeedSplitterOut1Ratio * m.LPCExt46MFlow \
                + 5.686396e-02 * m.FeedSplitterOut1Ratio * m.LPCExt15MFlow \
                + -4.110307e-01 * m.FeedSplitterOut2Ratio * m.FeedMFlow \
                + -6.429820e+01 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut1Ratio \
                + 1.265899e+03 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut2Ratio \
                + 7.688241e-02 * m.FeedSplitterOut2Ratio * m.HPCCondPrdtMFlow \
                + 1.044038e-01 * m.FeedSplitterOut2Ratio * m.LPCExt46MFlow \
                + -2.182886e-01 * m.FeedSplitterOut2Ratio * m.LPCExt15MFlow \
                + 1.825708e-06 * m.HPCCondPrdtMFlow * m.FeedMFlow \
                + -7.710441e-02 * m.HPCCondPrdtMFlow * m.FeedSplitterOut1Ratio \
                + 7.688241e-02 * m.HPCCondPrdtMFlow * m.FeedSplitterOut2Ratio \
                + 7.527786e-05 * m.HPCCondPrdtMFlow * m.HPCCondPrdtMFlow \
                + 1.446353e-05 * m.HPCCondPrdtMFlow * m.LPCExt46MFlow \
                + -3.642225e-05 * m.HPCCondPrdtMFlow * m.LPCExt15MFlow \
                + 2.471715e-06 * m.LPCExt46MFlow * m.FeedMFlow \
                + -1.328762e-01 * m.LPCExt46MFlow * m.FeedSplitterOut1Ratio \
                + 1.044038e-01 * m.LPCExt46MFlow * m.FeedSplitterOut2Ratio \
                + 1.446353e-05 * m.LPCExt46MFlow * m.HPCCondPrdtMFlow \
                + 1.804850e-05 * m.LPCExt46MFlow * m.LPCExt46MFlow \
                + -3.353461e-05 * m.LPCExt46MFlow * m.LPCExt15MFlow \
                + 3.139748e-05 * m.LPCExt15MFlow * m.FeedMFlow \
                + 5.686396e-02 * m.LPCExt15MFlow * m.FeedSplitterOut1Ratio \
                + -2.182886e-01 * m.LPCExt15MFlow * m.FeedSplitterOut2Ratio \
                + -3.642225e-05 * m.LPCExt15MFlow * m.HPCCondPrdtMFlow \
                + -3.353461e-05 * m.LPCExt15MFlow * m.LPCExt46MFlow \
                + 4.641023e-04 * m.LPCExt15MFlow * m.LPCExt15MFlow

        test_model.TurbineCon8 = Expression(rule=TurbineCon8)

        def TurbineCon9(m):
            return 3.208014e+03 \
                + -9.841714e-01 * m.FeedMFlow \
                + -7.558998e+02 * m.FeedSplitterOut1Ratio \
                + -5.291948e+02 * m.FeedSplitterOut2Ratio \
                + 2.295342e-02 * m.HPCCondPrdtMFlow \
                + -1.373395e-01 * m.LPCExt46MFlow \
                + 1.801168e-02 * m.LPCExt15MFlow \
                + 1.228876e-04 * m.FeedMFlow * m.FeedMFlow \
                + -4.787164e-01 * m.FeedMFlow * m.FeedSplitterOut1Ratio \
                + -1.639962e-02 * m.FeedMFlow * m.FeedSplitterOut2Ratio \
                + -1.429595e-05 * m.FeedMFlow * m.HPCCondPrdtMFlow \
                + -1.410488e-07 * m.FeedMFlow * m.LPCExt46MFlow \
                + 2.754627e-05 * m.FeedMFlow * m.LPCExt15MFlow \
                + -4.787164e-01 * m.FeedSplitterOut1Ratio * m.FeedMFlow \
                + 1.962893e+03 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut1Ratio \
                + -6.271229e+01 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut2Ratio \
                + 5.043181e-02 * m.FeedSplitterOut1Ratio * m.HPCCondPrdtMFlow \
                + -4.941073e-02 * m.FeedSplitterOut1Ratio * m.LPCExt46MFlow \
                + -1.117674e-01 * m.FeedSplitterOut1Ratio * m.LPCExt15MFlow \
                + -1.639962e-02 * m.FeedSplitterOut2Ratio * m.FeedMFlow \
                + -6.271229e+01 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut1Ratio \
                + 3.929170e+02 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut2Ratio \
                + 7.633996e-02 * m.FeedSplitterOut2Ratio * m.HPCCondPrdtMFlow \
                + 7.720145e-02 * m.FeedSplitterOut2Ratio * m.LPCExt46MFlow \
                + 2.494581e-02 * m.FeedSplitterOut2Ratio * m.LPCExt15MFlow \
                + -1.429595e-05 * m.HPCCondPrdtMFlow * m.FeedMFlow \
                + 5.043181e-02 * m.HPCCondPrdtMFlow * m.FeedSplitterOut1Ratio \
                + 7.633996e-02 * m.HPCCondPrdtMFlow * m.FeedSplitterOut2Ratio \
                + 3.484688e-05 * m.HPCCondPrdtMFlow * m.HPCCondPrdtMFlow \
                + 1.801919e-07 * m.HPCCondPrdtMFlow * m.LPCExt46MFlow \
                + 2.623980e-05 * m.HPCCondPrdtMFlow * m.LPCExt15MFlow \
                + -1.410488e-07 * m.LPCExt46MFlow * m.FeedMFlow \
                + -4.941073e-02 * m.LPCExt46MFlow * m.FeedSplitterOut1Ratio \
                + 7.720145e-02 * m.LPCExt46MFlow * m.FeedSplitterOut2Ratio \
                + 1.801919e-07 * m.LPCExt46MFlow * m.HPCCondPrdtMFlow \
                + 2.947846e-05 * m.LPCExt46MFlow * m.LPCExt46MFlow \
                + -4.174938e-06 * m.LPCExt46MFlow * m.LPCExt15MFlow \
                + 2.754627e-05 * m.LPCExt15MFlow * m.FeedMFlow \
                + -1.117674e-01 * m.LPCExt15MFlow * m.FeedSplitterOut1Ratio \
                + 2.494581e-02 * m.LPCExt15MFlow * m.FeedSplitterOut2Ratio \
                + 2.623980e-05 * m.LPCExt15MFlow * m.HPCCondPrdtMFlow \
                + -4.174938e-06 * m.LPCExt15MFlow * m.LPCExt46MFlow \
                + 9.961317e-05 * m.LPCExt15MFlow * m.LPCExt15MFlow

        test_model.TurbineCon9 = Expression(rule=TurbineCon9)

        def TurbineCon10(m):
            return 7.511160e+02 \
                + -9.990782e-01 * m.FeedMFlow \
                + -4.914530e+02 * m.FeedSplitterOut1Ratio \
                + -4.318461e+02 * m.FeedSplitterOut2Ratio \
                + -1.837085e-01 * m.HPCCondPrdtMFlow \
                + -1.954327e-01 * m.LPCExt46MFlow \
                + 6.761131e-02 * m.LPCExt15MFlow \
                + 1.271158e-04 * m.FeedMFlow * m.FeedMFlow \
                + 4.604325e-01 * m.FeedMFlow * m.FeedSplitterOut1Ratio \
                + -2.699428e-02 * m.FeedMFlow * m.FeedSplitterOut2Ratio \
                + -8.644853e-06 * m.FeedMFlow * m.HPCCondPrdtMFlow \
                + 3.257900e-06 * m.FeedMFlow * m.LPCExt46MFlow \
                + -4.146879e-06 * m.FeedMFlow * m.LPCExt15MFlow \
                + 4.604325e-01 * m.FeedSplitterOut1Ratio * m.FeedMFlow \
                + 1.850657e+03 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut1Ratio \
                + -4.588933e+01 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut2Ratio \
                + 6.050333e-02 * m.FeedSplitterOut1Ratio * m.HPCCondPrdtMFlow \
                + -1.447124e-02 * m.FeedSplitterOut1Ratio * m.LPCExt46MFlow \
                + -5.686985e-02 * m.FeedSplitterOut1Ratio * m.LPCExt15MFlow \
                + -2.699428e-02 * m.FeedSplitterOut2Ratio * m.FeedMFlow \
                + -4.588933e+01 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut1Ratio \
                + 4.172122e+02 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut2Ratio \
                + 8.315836e-02 * m.FeedSplitterOut2Ratio * m.HPCCondPrdtMFlow \
                + 8.431285e-02 * m.FeedSplitterOut2Ratio * m.LPCExt46MFlow \
                + -9.563644e-02 * m.FeedSplitterOut2Ratio * m.LPCExt15MFlow \
                + -8.644853e-06 * m.HPCCondPrdtMFlow * m.FeedMFlow \
                + 6.050333e-02 * m.HPCCondPrdtMFlow * m.FeedSplitterOut1Ratio \
                + 8.315836e-02 * m.HPCCondPrdtMFlow * m.FeedSplitterOut2Ratio \
                + 8.950745e-05 * m.HPCCondPrdtMFlow * m.HPCCondPrdtMFlow \
                + 8.110280e-06 * m.HPCCondPrdtMFlow * m.LPCExt46MFlow \
                + 2.096725e-05 * m.HPCCondPrdtMFlow * m.LPCExt15MFlow \
                + 3.257900e-06 * m.LPCExt46MFlow * m.FeedMFlow \
                + -1.447124e-02 * m.LPCExt46MFlow * m.FeedSplitterOut1Ratio \
                + 8.431285e-02 * m.LPCExt46MFlow * m.FeedSplitterOut2Ratio \
                + 8.110280e-06 * m.LPCExt46MFlow * m.HPCCondPrdtMFlow \
                + 3.851044e-05 * m.LPCExt46MFlow * m.LPCExt46MFlow \
                + -1.246923e-05 * m.LPCExt46MFlow * m.LPCExt15MFlow \
                + -4.146879e-06 * m.LPCExt15MFlow * m.FeedMFlow \
                + -5.686985e-02 * m.LPCExt15MFlow * m.FeedSplitterOut1Ratio \
                + -9.563644e-02 * m.LPCExt15MFlow * m.FeedSplitterOut2Ratio \
                + 2.096725e-05 * m.LPCExt15MFlow * m.HPCCondPrdtMFlow \
                + -1.246923e-05 * m.LPCExt15MFlow * m.LPCExt46MFlow \
                + 1.218963e-04 * m.LPCExt15MFlow * m.LPCExt15MFlow

        test_model.TurbineCon10 = Expression(rule=TurbineCon10)

        def OverallHeatBlnc(m):
            return 1.425871e+07 \
                + -3.266186e+03 * m.FeedMFlow \
                + -2.150964e+06 * m.FeedSplitterOut1Ratio \
                + -6.065051e+06 * m.FeedSplitterOut2Ratio \
                + 6.829814e+02 * m.HPCCondPrdtMFlow \
                + -6.835264e+02 * m.LPCExt46MFlow \
                + -8.192126e+03 * m.LPCExt15MFlow \
                + 2.649507e-01 * m.FeedMFlow * m.FeedMFlow \
                + 9.642475e+01 * m.FeedMFlow * m.FeedSplitterOut1Ratio \
                + 1.560305e+02 * m.FeedMFlow * m.FeedSplitterOut2Ratio \
                + -1.004180e-01 * m.FeedMFlow * m.HPCCondPrdtMFlow \
                + 4.842998e-02 * m.FeedMFlow * m.LPCExt46MFlow \
                + 5.567983e-01 * m.FeedMFlow * m.LPCExt15MFlow \
                + 9.642475e+01 * m.FeedSplitterOut1Ratio * m.FeedMFlow \
                + 7.753388e+06 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut1Ratio \
                + -5.833275e+04 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut2Ratio \
                + -3.448055e+02 * m.FeedSplitterOut1Ratio * m.HPCCondPrdtMFlow \
                + -3.060215e+02 * m.FeedSplitterOut1Ratio * m.LPCExt46MFlow \
                + -7.047711e+02 * m.FeedSplitterOut1Ratio * m.LPCExt15MFlow \
                + 1.560305e+02 * m.FeedSplitterOut2Ratio * m.FeedMFlow \
                + -5.833275e+04 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut1Ratio \
                + 4.234126e+06 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut2Ratio \
                + -1.384840e+02 * m.FeedSplitterOut2Ratio * m.HPCCondPrdtMFlow \
                + 2.068624e+02 * m.FeedSplitterOut2Ratio * m.LPCExt46MFlow \
                + 7.294073e+02 * m.FeedSplitterOut2Ratio * m.LPCExt15MFlow \
                + -1.004180e-01 * m.HPCCondPrdtMFlow * m.FeedMFlow \
                + -3.448055e+02 * m.HPCCondPrdtMFlow * m.FeedSplitterOut1Ratio \
                + -1.384840e+02 * m.HPCCondPrdtMFlow * m.FeedSplitterOut2Ratio \
                + 3.366985e-01 * m.HPCCondPrdtMFlow * m.HPCCondPrdtMFlow \
                + -4.424636e-02 * m.HPCCondPrdtMFlow * m.LPCExt46MFlow \
                + -2.030321e-02 * m.HPCCondPrdtMFlow * m.LPCExt15MFlow \
                + 4.842998e-02 * m.LPCExt46MFlow * m.FeedMFlow \
                + -3.060215e+02 * m.LPCExt46MFlow * m.FeedSplitterOut1Ratio \
                + 2.068624e+02 * m.LPCExt46MFlow * m.FeedSplitterOut2Ratio \
                + -4.424636e-02 * m.LPCExt46MFlow * m.HPCCondPrdtMFlow \
                + 9.954828e-02 * m.LPCExt46MFlow * m.LPCExt46MFlow \
                + -1.676427e-01 * m.LPCExt46MFlow * m.LPCExt15MFlow \
                + 5.567983e-01 * m.LPCExt15MFlow * m.FeedMFlow \
                + -7.047711e+02 * m.LPCExt15MFlow * m.FeedSplitterOut1Ratio \
                + 7.294073e+02 * m.LPCExt15MFlow * m.FeedSplitterOut2Ratio \
                + -2.030321e-02 * m.LPCExt15MFlow * m.HPCCondPrdtMFlow \
                + -1.676427e-01 * m.LPCExt15MFlow * m.LPCExt46MFlow \
                + 3.067898e+00 * m.LPCExt15MFlow * m.LPCExt15MFlow

        test_model.OverallHeatBlnc = Expression(rule=OverallHeatBlnc)

        def OBJ(m):
            return 2.490562e+04 \
                + -2.490993e+00 * m.FeedMFlow \
                + -6.915694e+03 * m.FeedSplitterOut1Ratio \
                + -1.039509e+04 * m.FeedSplitterOut2Ratio \
                + 1.684281e+00 * m.HPCCondPrdtMFlow \
                + -2.825967e+00 * m.LPCExt46MFlow \
                + -1.360548e+01 * m.LPCExt15MFlow \
                + 5.008749e-04 * m.FeedMFlow * m.FeedMFlow \
                + 1.188547e-02 * m.FeedMFlow * m.FeedSplitterOut1Ratio \
                + -5.129066e-01 * m.FeedMFlow * m.FeedSplitterOut2Ratio \
                + -7.370631e-05 * m.FeedMFlow * m.HPCCondPrdtMFlow \
                + 1.066254e-04 * m.FeedMFlow * m.LPCExt46MFlow \
                + 5.278739e-04 * m.FeedMFlow * m.LPCExt15MFlow \
                + 1.188547e-02 * m.FeedSplitterOut1Ratio * m.FeedMFlow \
                + 9.025538e+03 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut1Ratio \
                + -1.613370e+03 * m.FeedSplitterOut1Ratio * m.FeedSplitterOut2Ratio \
                + -5.675858e-01 * m.FeedSplitterOut1Ratio * m.HPCCondPrdtMFlow \
                + 9.963293e-02 * m.FeedSplitterOut1Ratio * m.LPCExt46MFlow \
                + 4.123504e-01 * m.FeedSplitterOut1Ratio * m.LPCExt15MFlow \
                + -5.129066e-01 * m.FeedSplitterOut2Ratio * m.FeedMFlow \
                + -1.613370e+03 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut1Ratio \
                + 6.550452e+03 * m.FeedSplitterOut2Ratio * m.FeedSplitterOut2Ratio \
                + -2.807381e-01 * m.FeedSplitterOut2Ratio * m.HPCCondPrdtMFlow \
                + 4.054410e-01 * m.FeedSplitterOut2Ratio * m.LPCExt46MFlow \
                + 1.378650e+00 * m.FeedSplitterOut2Ratio * m.LPCExt15MFlow \
                + -7.370631e-05 * m.HPCCondPrdtMFlow * m.FeedMFlow \
                + -5.675858e-01 * m.HPCCondPrdtMFlow * m.FeedSplitterOut1Ratio \
                + -2.807381e-01 * m.HPCCondPrdtMFlow * m.FeedSplitterOut2Ratio \
                + 1.641286e-04 * m.HPCCondPrdtMFlow * m.HPCCondPrdtMFlow \
                + -6.895480e-05 * m.HPCCondPrdtMFlow * m.LPCExt46MFlow \
                + 3.732744e-04 * m.HPCCondPrdtMFlow * m.LPCExt15MFlow \
                + 1.066254e-04 * m.LPCExt46MFlow * m.FeedMFlow \
                + 9.963293e-02 * m.LPCExt46MFlow * m.FeedSplitterOut1Ratio \
                + 4.054410e-01 * m.LPCExt46MFlow * m.FeedSplitterOut2Ratio \
                + -6.895480e-05 * m.LPCExt46MFlow * m.HPCCondPrdtMFlow \
                + 4.149032e-04 * m.LPCExt46MFlow * m.LPCExt46MFlow \
                + -1.471105e-04 * m.LPCExt46MFlow * m.LPCExt15MFlow \
                + 5.278739e-04 * m.LPCExt15MFlow * m.FeedMFlow \
                + 4.123504e-01 * m.LPCExt15MFlow * m.FeedSplitterOut1Ratio \
                + 1.378650e+00 * m.LPCExt15MFlow * m.FeedSplitterOut2Ratio \
                + 3.732744e-04 * m.LPCExt15MFlow * m.HPCCondPrdtMFlow \
                + -1.471105e-04 * m.LPCExt15MFlow * m.LPCExt46MFlow \
                + 7.478130e-03 * m.LPCExt15MFlow * m.LPCExt15MFlow

        test_model.OBJ = Expression(rule=OBJ)



