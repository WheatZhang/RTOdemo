from rtolib.core.pyomo_model import PyomoModel
from pyomo.environ import *
import os

class ASU_QuadObjLinearCon_Model(PyomoModel):
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
            return 0

        test_model.PurityCon1 = Expression(rule=PurityCon1)

        def PurityCon2(m):
            return 0

        test_model.PurityCon2 = Expression(rule=PurityCon2)

        def PurityCon3(m):
            return 0

        test_model.PurityCon3 = Expression(rule=PurityCon3)

        def PurityCon4(m):
            return 0

        test_model.PurityCon4 = Expression(rule=PurityCon4)

        def PurityCon5(m):
            return 0

        test_model.PurityCon5 = Expression(rule=PurityCon5)

        def PurityCon6(m):
            return 0

        test_model.PurityCon6 = Expression(rule=PurityCon6)

        def ProductCon1(m):
            return 0

        test_model.ProductCon1 = Expression(rule=ProductCon1)

        def ProductCon2(m):
            return 0

        test_model.ProductCon2 = Expression(rule=ProductCon2)

        def ProductCon3(m):
            return 0

        test_model.ProductCon3 = Expression(rule=ProductCon3)

        def TempDiffCon1(m):
            return 0

        test_model.TempDiffCon1 = Expression(rule=TempDiffCon1)

        def TempDiffCon2(m):
            return 0

        test_model.TempDiffCon2 = Expression(rule=TempDiffCon2)

        def TurbineCon1(m):
            return 0

        test_model.TurbineCon1 = Expression(rule=TurbineCon1)

        def TurbineCon2(m):
            return 0

        test_model.TurbineCon2 = Expression(rule=TurbineCon2)

        def TurbineCon3(m):
            return 0

        test_model.TurbineCon3 = Expression(rule=TurbineCon3)

        def TurbineCon4(m):
            return 0

        test_model.TurbineCon4 = Expression(rule=TurbineCon4)

        def TurbineCon7(m):
            return 0

        test_model.TurbineCon7 = Expression(rule=TurbineCon7)

        def TurbineCon8(m):
            return 0

        test_model.TurbineCon8 = Expression(rule=TurbineCon8)

        def TurbineCon9(m):
            return 0

        test_model.TurbineCon9 = Expression(rule=TurbineCon9)

        def TurbineCon10(m):
            return 0

        test_model.TurbineCon10 = Expression(rule=TurbineCon10)

        def OverallHeatBlnc(m):
            return 0

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



