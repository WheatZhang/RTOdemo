from rtolib.core import ProblemDescription


symbol_list = {
                'MV': ('FeedMFlow', 'TAFeedFrac', 'MAFeedFrac', '40_LIN',\
                       '52_WN', '47_ARC', 'ASC_RefluxRatio', 'LOX_Frac'),# input
                'CV': (),# output
                'OBJ': "cost",# objective function
                'CON': ("OxygenPrdtPurityCon","OxygenPrdtImpurityCon",\
                        "NitrogenPrdtPurityCon","NitrogenPrdtImpurityCon",\
                        "ArgonPrdtImpurityCon","CrudeArgonImpurityCon",\
                        "LOX_Con","GAN_Con","GAR_Con",\
                        "MainCoolerTempDiffCon", "HeatLPCTempDiffCon",\
                        "HPA_GOX_LB_Con","HPA_GOX_UB_Con",\
                        "Feed_LB_Con","Feed_UB_Con",\
                        "HPA_LB_Con", "HPA_UB_Con",\
                        "TA_LB_Con", "TA_UB_Con",\
                        "OverallHeatCon"),  # inequality constraint
                'SPEC': None,  # specification
            }
bounds = {
        'FeedMFlow':(3800,4200),
        'TAFeedFrac':(0.15,0.25),
        'MAFeedFrac':(0.45,0.6),
        '40_LIN':(1200,1500),
        '52_WN':(1600,1900),
        '47_ARC':(600,700),
        'ASC_RefluxRatio':(25,30),
        'LOX_Frac':(0.07,0.1),
        "cost":(None, None),
        "OxygenPrdtPurityCon":(0.99,1),
        "OxygenPrdtImpurityCon":(0,1e-4),
        "NitrogenPrdtPurityCon":(0.99,1),
        "NitrogenPrdtImpurityCon":(0,1e-4),
        "ArgonPrdtImpurityCon":(0,1e-4),
        "CrudeArgonImpurityCon":(0,1e-2),
        "LOX_Con":(600,900),
        "GAN_Con":(1000,1800),
        "GAR_Con":(10,40),
        "MainCoolerTempDiffCon":(0,10),
        "HeatLPCTempDiffCon":(0,10),
        "HPA_GOX_LB_Con":(-300,300),
        "HPA_GOX_UB_Con":(-300,300),
        "Feed_LB_Con":(-1000,1000),
        "Feed_UB_Con":(-1000,1000),
        "HPA_LB_Con":(-400,400),
        "HPA_UB_Con":(-400,400),
        "TA_LB_Con":(-400,400),
        "TA_UB_Con":(-400,400),
        "OverallHeatCon":(-2e7,2e7),
    }
scaling_factors = {
        'FeedMFlow':100,
        'TAFeedFrac':0.1,
        'MAFeedFrac':0.1,
        '40_LIN':100,
        '52_WN':100,
        '47_ARC':20,
        'ASC_RefluxRatio':1,
        'LOX_Frac':0.01,
        "cost":10,
        "OxygenPrdtPurityCon":1e-3,
        "OxygenPrdtImpurityCon":1e-3,
        "NitrogenPrdtPurityCon":1e-3,
        "NitrogenPrdtImpurityCon":1e-3,
        "ArgonPrdtImpurityCon":1e-3,
        "CrudeArgonImpurityCon":1e-2,
        "LOX_Con":100,
        "GAN_Con":100,
        "GAR_Con":10,
        "MainCoolerTempDiffCon":0.5,
        "HeatLPCTempDiffCon":0.5,
        "HPA_GOX_LB_Con":100,
        "HPA_GOX_UB_Con":100,
        "Feed_LB_Con":100,
        "Feed_UB_Con":100,
        "HPA_LB_Con":100,
        "HPA_UB_Con":100,
        "TA_LB_Con":100,
        "TA_UB_Con":100,
        "OverallHeatCon":1e7,
    }
default_values = {
        'FeedMFlow':4100,
        'TAFeedFrac':0.19,
        'MAFeedFrac':0.5,
        '40_LIN':1450,
        '52_WN':1800,
        '47_ARC':600,
        'ASC_RefluxRatio':28,
        'LOX_Frac':0.09,
    }
default_ASU_description=ProblemDescription(symbol_list, bounds, scaling_factors, default_values)
