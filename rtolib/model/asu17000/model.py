from rtolib.core.pyomo_model import PyomoModel
from pyomo.environ import *
import os

class ASU_Model(PyomoModel):
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
        binary_coeff = -0.015
        LOX = 17000

        test_model.Component = Set(initialize=['Nitrogen', 'Oxygen', 'Argon'])

        # ===================================
        #
        #      Process Parts
        #
        # ===================================
        test_model.PR_CrtcTemp = Param(test_model.Component,
                                       initialize={'Oxygen': 154.77, 'Argon': 150.71,
                                                   'Nitrogen': 126.15})  # Temp unit: K
        test_model.PR_CrtcPressure = Param(test_model.Component, initialize={'Oxygen': 5080000, 'Argon': 4864000,
                                                                             'Nitrogen': 3394000})  # Crtc pressure, pa
        test_model.PR_CompOmega = Param(test_model.Component,
                                        initialize={'Oxygen': 0.019, 'Argon': 0, 'Nitrogen': 0.045})
        test_model.PR_SatPresCorC = Param(test_model.Component, initialize={'Oxygen': -5.667, 'Nitrogen': -6.344})
        # binary_coeff = -0.015  # -0.01238
        dict = {
            ('Oxygen', 'Oxygen'): 0,
            ('Oxygen', 'Argon'): 0.0265,
            ('Oxygen', 'Nitrogen'): binary_coeff,
            ('Argon', 'Oxygen'): 0.0265,
            ('Argon', 'Argon'): 0,
            ('Argon', 'Nitrogen'): -0.004071,
            ('Nitrogen', 'Oxygen'): binary_coeff,
            ('Nitrogen', 'Argon'): -0.004071,
            ('Nitrogen', 'Nitrogen'): 0
        }
        test_model.PR_Kij = Param(test_model.Component, test_model.Component, initialize=dict)
        test_model.PR_RefTemp = Param(test_model.Component,
                                      initialize={'Oxygen': 298.15, 'Argon': 298.15, 'Nitrogen': 298.15})  # unit: K
        # Reference state enthalpy for ideal gas at reference temperature, cal/mol
        test_model.PR_DHFORM = Param(test_model.Component, initialize={'Oxygen': 0, 'Argon': 0, 'Nitrogen': 0})
        # Parameters to calculate the equation for the DIPPR ideal gas heat capacity model by Aly and Lee 1981,unit: cal/(mol*k)
        test_model.PR_CPIGDP1 = Param(test_model.Component,
                                      initialize={'Oxygen': 1 / 4180 * 29103, 'Argon': 1 / 4180 * 20786,
                                                  'Nitrogen': 1 / 4180 * 29105})
        test_model.PR_CPIGDP2 = Param(test_model.Component,
                                      initialize={'Oxygen': 1 / 4180 * 10040, 'Argon': 1 / 4180 * 0,
                                                  'Nitrogen': 1 / 4180 * 8614.9})
        test_model.PR_CPIGDP3 = Param(test_model.Component,
                                      initialize={'Oxygen': 2526.5, 'Argon': 0, 'Nitrogen': 1701.6})
        test_model.PR_CPIGDP4 = Param(test_model.Component,
                                      initialize={'Oxygen': 1 / 4180 * 9356, 'Argon': 1 / 4180 * 0,
                                                  'Nitrogen': 1 / 4180 * 103.47})
        test_model.PR_CPIGDP5 = Param(test_model.Component,
                                      initialize={'Oxygen': 1153.8, 'Argon': 0, 'Nitrogen': 909.79})
        test_model.PR_CPIGDP6 = Param(test_model.Component, initialize={'Oxygen': 50, 'Argon': -173.15, 'Nitrogen': 50})
        test_model.PR_CPIGDP7 = Param(test_model.Component,
                                      initialize={'Oxygen': 1500, 'Argon': 1226.85, 'Nitrogen': 1500})
        test_model.PR_RGas = Param(initialize=8.314)  # J/(mol*K)

        def PR_bi_Calc(m, comp):
            return 0.07780 * m.PR_RGas * m.PR_CrtcTemp[comp] / m.PR_CrtcPressure[comp]

        test_model.PR_bi = Param(test_model.Component, initialize=PR_bi_Calc)

        def PR_m_Calc(m, comp):
            return 0.37464 + 1.54226 * m.PR_CompOmega[comp] - 0.26992 * m.PR_CompOmega[comp] * m.PR_CompOmega[comp]

        test_model.PR_m = Param(test_model.Component, initialize=PR_m_Calc)
        test_model.SRK_CrtcTemp = Param(test_model.Component,
                                        initialize={'Oxygen': 154.58, 'Argon': 150.86,
                                                    'Nitrogen': 126.2})  # Temp unit: K
        test_model.SRK_CrtcPressure = Param(test_model.Component,
                                            initialize={'Oxygen': 5042995.9125, 'Argon': 4897999.8375,
                                                        'Nitrogen': 3400000.905})  # Crtc pressure, pa
        test_model.SRK_CompOmega = Param(test_model.Component,
                                         initialize={'Oxygen': 0.0221798, 'Argon': 0, 'Nitrogen': 0.0377215})
        dict = {
            ('Oxygen', 'Oxygen'): 0,
            ('Oxygen', 'Argon'): 0.01505,
            ('Oxygen', 'Nitrogen'): -0.0078,
            ('Argon', 'Oxygen'): 0.01505,
            ('Argon', 'Argon'): 0,
            ('Argon', 'Nitrogen'): -0.0016,
            ('Nitrogen', 'Oxygen'): -0.0078,
            ('Nitrogen', 'Argon'): -0.0016,
            ('Nitrogen', 'Nitrogen'): 0
        }
        test_model.SRK_Kij = Param(test_model.Component, test_model.Component, initialize=dict)
        test_model.SRK_RefTemp = Param(test_model.Component,
                                       initialize={'Oxygen': 298.15, 'Argon': 298.15, 'Nitrogen': 298.15})  # ??λK
        # Reference state enthalpy for ideal gas at reference temperature, cal/mol
        test_model.SRK_DHFORM = Param(test_model.Component, initialize={'Oxygen': 0, 'Argon': 0, 'Nitrogen': 0})
        # Parameters to calculate the equation for the DIPPR ideal gas heat capacity model by Aly and Lee 1981,unit: cal/(mol*k)
        test_model.SRK_CPIGDP1 = Param(test_model.Component,
                                       initialize={'Oxygen': 1 / 4180 * 29103, 'Argon': 1 / 4180 * 20786,
                                                   'Nitrogen': 1 / 4180 * 29105})
        test_model.SRK_CPIGDP2 = Param(test_model.Component,
                                       initialize={'Oxygen': 1 / 4180 * 10040, 'Argon': 1 / 4180 * 0,
                                                   'Nitrogen': 1 / 4180 * 8614.9})
        test_model.SRK_CPIGDP3 = Param(test_model.Component,
                                       initialize={'Oxygen': 2526.5, 'Argon': 0, 'Nitrogen': 1701.6})
        test_model.SRK_CPIGDP4 = Param(test_model.Component,
                                       initialize={'Oxygen': 1 / 4180 * 9356, 'Argon': 1 / 4180 * 0,
                                                   'Nitrogen': 1 / 4180 * 103.47})
        test_model.SRK_CPIGDP5 = Param(test_model.Component,
                                       initialize={'Oxygen': 1153.8, 'Argon': 0, 'Nitrogen': 909.79})
        test_model.SRK_CPIGDP6 = Param(test_model.Component,
                                       initialize={'Oxygen': 50, 'Argon': -173.15, 'Nitrogen': 50})
        test_model.SRK_CPIGDP7 = Param(test_model.Component,
                                       initialize={'Oxygen': 1500, 'Argon': 1226.85, 'Nitrogen': 1500})
        test_model.SRK_RGas = Param(initialize=8.314)  # J/(mol*K)

        def SRK_bi_Calc(m, comp):
            return 0.08664 * m.SRK_RGas * m.SRK_CrtcTemp[comp] / m.SRK_CrtcPressure[comp]

        test_model.SRK_bi = Param(test_model.Component, initialize=SRK_bi_Calc)

        def SRK_ai_Calc(m, comp):
            return 0.42747 * m.SRK_RGas * m.SRK_RGas * m.SRK_CrtcTemp[comp] * m.SRK_CrtcTemp[comp] / m.SRK_CrtcPressure[
                comp]

        test_model.SRK_ai = Param(test_model.Component, initialize=SRK_ai_Calc)

        def SRK_m_Calc(m, comp):
            return 0.48508 + 1.55171 * m.SRK_CompOmega[comp] - 0.15613 * m.SRK_CompOmega[comp] * m.SRK_CompOmega[comp]

        test_model.SRK_m = Param(test_model.Component, initialize=SRK_m_Calc)

        # ===================================
        #
        #      Parameters & Sets
        #
        # ===================================
        # --------Feed---------
        test_model.FeedMFlow = Var(initialize=4030.759821, bounds=(0, 8061.5196428571435))
        test_model.FeedMFrac = Param(test_model.Component,
                                     initialize={'Nitrogen': 0.78118, 'Oxygen': 0.2095, 'Argon': 0.00932})
        test_model.FeedPressure = Param(initialize=100.000000)
        test_model.FeedTemp = Param(initialize=39.850000)
        test_model.FeedMEtlp = Param(initialize=424.928400)
        # --------FeedSplitter---------
        test_model.FeedSplitterOut1Ratio = Var(initialize=0.184000, bounds=(0, 1))
        test_model.FeedSplitterOut2Ratio = Var(initialize=0.517000, bounds=(0, 1))
        # --------TaCooler---------
        test_model.TaCoolerPressure = Param(initialize=566.000000)
        test_model.TaCoolerTemp = Param(initialize=-173.370000)
        # --------HpaCooler---------
        test_model.HpaCoolerPressure = Param(initialize=564.000000)
        test_model.HpaCoolerVF = Param(initialize=0.080121)
        # --------MaCooler---------
        test_model.MaCoolerTemp = Param(initialize=-168.550000)
        test_model.MaCoolerPressure = Param(initialize=566.000000)
        # --------HPCZeroReboiled---------
        test_model.HPCZeroReboiledNull = Param(initialize=0.000000)
        # --------HPC---------
        test_model.HPCTrays = RangeSet(1, 41)
        test_model.HPCTopPressure = Param(initialize=539.219512)
        test_model.HPCBtmPressure = Param(initialize=564.000000)
        # --------HPCSump---------
        test_model.HPCSumpLiqRho = Param(initialize=27.000000)
        test_model.HPCSumpSumpCSArea = Param(initialize=6.000000)
        # --------HPCCond---------
        test_model.HPCCondPressure = Param(initialize=538.600000)
        test_model.HPCCondPrdtMFlow = Var(initialize=1381.050000, bounds=(690.5250000000001, 2762.1000000000004))
        # --------CoolLin---------
        test_model.CoolLinTemp = Param(initialize=-191.150000)
        test_model.CoolLinDeltaP = Param(initialize=0.000000)
        # --------QTrtl---------
        test_model.QTrtlTemp = Param(initialize=-175.000000)
        test_model.QTrtlHeat = Param(initialize=0.000000)
        # --------Throttle---------
        test_model.ThrottlePressure = Param(initialize=140.000000)
        # --------HeatLPC---------
        test_model.HeatLPCDeltaP = Param(initialize=0.000000)
        # --------LPCZeroReflux---------
        test_model.LPCZeroRefluxNull = Param(initialize=0.000000)
        # --------LPC---------
        test_model.LPCTrays = RangeSet(1, 52)
        test_model.LPCExt46MFlow = Var(initialize=1750.840179, bounds=(875.4200892857143, 3501.6803571428572))
        test_model.LPCExt15MFlow = Var(initialize=647.429911, bounds=(323.7149553571429, 1294.8598214285716))
        test_model.LPCTopPressure = Param(initialize=135.600000)
        test_model.LPCBtmPressure = Param(initialize=140.896154)
        # --------LPCReboiler---------
        test_model.LPCReboilerPressure = Param(initialize=141.000000)
        test_model.LPCReboilerLiqRho = Param(initialize=30282.000000)
        test_model.LPCReboilerRebCSArea = Param(initialize=6.000000)
        # --------LPCReboiler VLE & Enthalpy---------
        # --------ZeroReboiledASC---------
        test_model.ZeroReboiledASCNull = Param(initialize=0.000000)
        # --------ASC---------
        test_model.ASCTrays = RangeSet(1, 190)
        test_model.ASCTopPressure = Param(initialize=124.476842)
        test_model.ASCBtmPressure = Param(initialize=139.000000)
        # --------ASCSump---------
        test_model.ASCSumpLiqRho = Param(initialize=34.522040)
        test_model.ASCSumpSumpCSArea = Param(initialize=6.000000)
        # --------ASCCond---------
        test_model.ASCCondPressure = Param(initialize=124.400000)
        test_model.ASCCondRefRatio = Var(initialize=27.966300, bounds=(1, 100))
        # --------HeatGan---------
        test_model.HeatGanDeltaP = Param(initialize=0.000000)
        # --------RecoGan---------
        test_model.RecoGanTemp = Param(initialize=32.850000)
        test_model.RecoGanPressure = Param(initialize=135.600000)
        # --------RecoGar---------
        test_model.RecoGarTemp = Param(initialize=32.850000)
        test_model.RecoGarPressure = Param(initialize=124.400000)
        # --------RecoWn---------
        test_model.RecoWnTemp = Param(initialize=32.850000)
        test_model.RecoWnPressure = Param(initialize=124.400000)
        # --------OxSplitter---------
        test_model.OxSplitterOutARatio = Var(initialize=0.898074, bounds=(0, 1))
        # --------RecoGox---------
        test_model.RecoGoxTemp = Param(initialize=32.850000)
        test_model.RecoGoxPressure = Param(initialize=141.000000)

        # ===================================
        #
        #         Variables
        #
        # ===================================
        # --------FeedSplitter---------
        test_model.FeedSplitterOut1MFlow = Var(within=NonNegativeReals)
        test_model.FeedSplitterOut2MFlow = Var(within=NonNegativeReals)
        test_model.FeedSplitterOut3MFlow = Var(within=NonNegativeReals)
        # --------TaCooler---------
        test_model.TaCoolerVapMFlow = Var(within=NonNegativeReals)
        test_model.TaCoolerLiqMFlow = Var(within=NonNegativeReals)
        test_model.TaCoolerMHeat = Var(within=NonNegativeReals)
        test_model.TaCoolerLiqMFrac = Var(test_model.Component, bounds=(0, 1))
        test_model.TaCoolerVapMFrac = Var(test_model.Component, bounds=(0, 1))
        test_model.TaCoolerSelfSqrtTr = Var(test_model.Component, initialize=0.8)
        test_model.TaCoolerSelfSqrt_ai = Var(test_model.Component, test_model.Component, initialize=0.17,
                                             bounds=(0.05, 0.4))
        test_model.TaCoolerSelfZL = Var(initialize=0.023966, bounds=(0, 0.05))
        test_model.TaCoolerSelfZV = Var(initialize=0.857334, bounds=(0.2, 3))
        test_model.TaCoolerSelfZLBL = Var(initialize=-4.878610, bounds=(-10, -1))
        test_model.TaCoolerSelfZLdivBL = Var(initialize=1.30898, bounds=(0.8, 2))
        test_model.TaCoolerSelfZVBV = Var(initialize=-0.1732, bounds=(-0.3, 0))
        test_model.TaCoolerSelfZVdivBV = Var(initialize=0.053666900, bounds=(0, 0.2))
        test_model.TaCoolerSelfVLE_K = Var(test_model.Component, bounds=(0, 8))
        # --------HpaCooler---------
        test_model.HpaCoolerTemp = Var()
        test_model.HpaCoolerVapMFlow = Var(within=NonNegativeReals)
        test_model.HpaCoolerLiqMFlow = Var(within=NonNegativeReals)
        test_model.HpaCoolerMHeat = Var(within=NonNegativeReals)
        test_model.HpaCoolerLiqMFrac = Var(test_model.Component, bounds=(0, 1))
        test_model.HpaCoolerVapMFrac = Var(test_model.Component, bounds=(0, 1))
        test_model.HpaCoolerSelfSqrtTr = Var(test_model.Component, initialize=0.8)
        test_model.HpaCoolerSelfSqrt_ai = Var(test_model.Component, test_model.Component, initialize=0.17,
                                              bounds=(0.05, 0.4))
        test_model.HpaCoolerSelfZL = Var(initialize=0.023966, bounds=(0, 0.05))
        test_model.HpaCoolerSelfZV = Var(initialize=0.857334, bounds=(0.2, 3))
        test_model.HpaCoolerSelfZLBL = Var(initialize=-4.878610, bounds=(-10, -1))
        test_model.HpaCoolerSelfZLdivBL = Var(initialize=1.30898, bounds=(0.8, 2))
        test_model.HpaCoolerSelfZVBV = Var(initialize=-0.1732, bounds=(-0.3, 0))
        test_model.HpaCoolerSelfZVdivBV = Var(initialize=0.053666900, bounds=(0, 0.2))
        test_model.HpaCoolerSelfVLE_K = Var(test_model.Component, bounds=(0, 8))
        # --------MaCooler---------
        test_model.MaCoolerMHeat = Var(within=NonNegativeReals)
        test_model.MaCoolerSelfSqrtTr = Var(test_model.Component, initialize=0.8)
        test_model.MaCoolerSelfSqrt_ai = Var(test_model.Component, test_model.Component, initialize=0.17,
                                             bounds=(0.05, 0.4))
        test_model.MaCoolerSelfZV = Var(initialize=0.866923293, bounds=(0.2, 3))
        test_model.MaCoolerSelfZVdivBV = Var(initialize=0.050441702, bounds=(0, 0.2))
        # --------AirMixer---------
        test_model.AirMixerOutMFlow = Var(within=NonNegativeReals)
        test_model.AirMixerOutMFrac = Var(test_model.Component, bounds=(0, 1))
        test_model.AirMixerOutMEtlp = Var()
        # --------HPC---------
        test_model.HPCVapLvMFlow = Var(test_model.HPCTrays, within=NonNegativeReals)
        test_model.HPCLiqLvMFlow = Var(test_model.HPCTrays, within=NonNegativeReals)
        test_model.HPCTrayTemp = Var(test_model.HPCTrays)
        test_model.HPCTrayPressure = Var(test_model.HPCTrays)
        test_model.HPCLiqLvMFrac = Var(test_model.HPCTrays, test_model.Component, bounds=(0, 1))
        test_model.HPCVapLvMFrac = Var(test_model.HPCTrays, test_model.Component, bounds=(0, 1))
        test_model.HPCAllTraysSqrtTr = Var(test_model.HPCTrays, test_model.Component, initialize=0.8)
        test_model.HPCAllTraysSqrt_ai = Var(test_model.HPCTrays, test_model.Component, test_model.Component,
                                            initialize=0.17, bounds=(0.05, 0.4))
        test_model.HPCAllTraysZL = Var(test_model.HPCTrays, initialize=0.023966, bounds=(0, 0.05))
        test_model.HPCAllTraysZV = Var(test_model.HPCTrays, initialize=0.857334, bounds=(0.2, 3))
        test_model.HPCAllTraysZLBL = Var(test_model.HPCTrays, initialize=-4.878610, bounds=(-10, -1))
        test_model.HPCAllTraysZLdivBL = Var(test_model.HPCTrays, initialize=1.30898, bounds=(0.8, 2))
        test_model.HPCAllTraysZVBV = Var(test_model.HPCTrays, initialize=-0.1732, bounds=(-0.3, 0))
        test_model.HPCAllTraysZVdivBV = Var(test_model.HPCTrays, initialize=0.053666900, bounds=(0, 0.2))
        test_model.HPCAllTraysVLE_K = Var(test_model.HPCTrays, test_model.Component, bounds=(0, 8))
        # --------HPCSump---------
        test_model.HPCSumpOutMFlow = Var(within=NonNegativeReals)
        test_model.HPCSumpHldpMFrac = Var(test_model.Component, bounds=(0, 1))
        # --------HPCCond---------
        test_model.HPCCondRefMFlow = Var(within=NonNegativeReals)
        test_model.HPCCondMHeatOut = Var(within=NonNegativeReals)
        test_model.HPCCondOutletSqrtTr = Var(test_model.Component, initialize=0.8)
        test_model.HPCCondOutletSqrt_ai = Var(test_model.Component, test_model.Component, initialize=0.17,
                                              bounds=(0.05, 0.4))
        test_model.HPCCondOutletZL = Var(initialize=0.023916, bounds=(0, 0.05))
        test_model.HPCCondOutletZLdivBL = Var(initialize=1.3089, bounds=(0.8, 2))
        # --------CoolLin---------
        test_model.CoolLinPressure = Var()
        test_model.CoolLinMHeat = Var(within=NonNegativeReals)
        test_model.CoolLinSelfSqrtTr = Var(test_model.Component, initialize=0.8)
        test_model.CoolLinSelfSqrt_ai = Var(test_model.Component, test_model.Component, initialize=0.17,
                                            bounds=(0.05, 0.4))
        test_model.CoolLinSelfZL = Var(initialize=0.023916, bounds=(0, 0.05))
        test_model.CoolLinSelfZLdivBL = Var(initialize=1.3089, bounds=(0.8, 2))
        # --------Throttle---------
        test_model.ThrottleTemp = Var()
        test_model.ThrottleVapMFlow = Var(within=NonNegativeReals)
        test_model.ThrottleLiqMFlow = Var(within=NonNegativeReals)
        test_model.ThrottleLiqMFrac = Var(test_model.Component, bounds=(0, 1))
        test_model.ThrottleVapMFrac = Var(test_model.Component, bounds=(0, 1))
        test_model.ThrottleSelfSqrtTr = Var(test_model.Component, initialize=0.8)
        test_model.ThrottleSelfSqrt_ai = Var(test_model.Component, test_model.Component, initialize=0.17,
                                             bounds=(0.05, 0.4))
        test_model.ThrottleSelfZL = Var(initialize=0.023966, bounds=(0, 0.05))
        test_model.ThrottleSelfZV = Var(initialize=0.857334, bounds=(0.2, 3))
        test_model.ThrottleSelfZLBL = Var(initialize=-4.878610, bounds=(-10, -1))
        test_model.ThrottleSelfZLdivBL = Var(initialize=1.30898, bounds=(0.8, 2))
        test_model.ThrottleSelfZVBV = Var(initialize=-0.1732, bounds=(-0.3, 0))
        test_model.ThrottleSelfZVdivBV = Var(initialize=0.053666900, bounds=(0, 0.2))
        test_model.ThrottleSelfVLE_K = Var(test_model.Component, bounds=(0, 8))
        # --------HeatLPC---------
        test_model.HeatLPCTemp = Var()
        test_model.HeatLPCPressure = Var()
        test_model.HeatLPCVapMFlow = Var(within=NonNegativeReals)
        test_model.HeatLPCLiqMFlow = Var(within=NonNegativeReals)
        test_model.HeatLPCLiqMFrac = Var(test_model.Component, bounds=(0, 1))
        test_model.HeatLPCVapMFrac = Var(test_model.Component, bounds=(0, 1))
        test_model.HeatLPCSelfSqrtTr = Var(test_model.Component, initialize=0.8)
        test_model.HeatLPCSelfSqrt_ai = Var(test_model.Component, test_model.Component, initialize=0.17,
                                            bounds=(0.05, 0.4))
        test_model.HeatLPCSelfZL = Var(initialize=0.023966, bounds=(0, 0.05))
        test_model.HeatLPCSelfZV = Var(initialize=0.857334, bounds=(0.2, 3))
        test_model.HeatLPCSelfZLBL = Var(initialize=-4.878610, bounds=(-10, -1))
        test_model.HeatLPCSelfZLdivBL = Var(initialize=1.30898, bounds=(0.8, 2))
        test_model.HeatLPCSelfZVBV = Var(initialize=-0.1732, bounds=(-0.3, 0))
        test_model.HeatLPCSelfZVdivBV = Var(initialize=0.053666900, bounds=(0, 0.2))
        test_model.HeatLPCSelfVLE_K = Var(test_model.Component, bounds=(0, 8))
        # --------LPC---------
        test_model.LPCVapLvMFlow = Var(test_model.LPCTrays, within=NonNegativeReals)
        test_model.LPCLiqLvMFlow = Var(test_model.LPCTrays, within=NonNegativeReals)
        test_model.LPCTrayTemp = Var(test_model.LPCTrays)
        test_model.LPCTrayPressure = Var(test_model.LPCTrays)
        test_model.LPCLiqLvMFrac = Var(test_model.LPCTrays, test_model.Component, bounds=(0, 1))
        test_model.LPCVapLvMFrac = Var(test_model.LPCTrays, test_model.Component, bounds=(0, 1))
        test_model.LPCAllTraysSqrtTr = Var(test_model.LPCTrays, test_model.Component, initialize=0.8)
        test_model.LPCAllTraysSqrt_ai = Var(test_model.LPCTrays, test_model.Component, test_model.Component,
                                            initialize=0.17, bounds=(0.05, 0.4))
        test_model.LPCAllTraysZL = Var(test_model.LPCTrays, initialize=0.023966, bounds=(0, 0.05))
        test_model.LPCAllTraysZV = Var(test_model.LPCTrays, initialize=0.857334, bounds=(0.2, 3))
        test_model.LPCAllTraysZLBL = Var(test_model.LPCTrays, initialize=-4.878610, bounds=(-10, -1))
        test_model.LPCAllTraysZLdivBL = Var(test_model.LPCTrays, initialize=1.30898, bounds=(0.8, 2))
        test_model.LPCAllTraysZVBV = Var(test_model.LPCTrays, initialize=-0.1732, bounds=(-0.3, 0))
        test_model.LPCAllTraysZVdivBV = Var(test_model.LPCTrays, initialize=0.053666900, bounds=(0, 0.2))
        test_model.LPCAllTraysVLE_K = Var(test_model.LPCTrays, test_model.Component, bounds=(0, 8))
        # --------LPCReboiler---------
        test_model.LPCReboilerLiqLvMFlow = Var(within=NonNegativeReals)
        test_model.LPCReboilerVapLvMFlow = Var(within=NonNegativeReals)
        test_model.LPCReboilerTemp = Var()
        test_model.LPCReboilerLiqLvMFrac = Var(test_model.Component, bounds=(0, 1))
        test_model.LPCReboilerVapLvMFrac = Var(test_model.Component, bounds=(0, 1))
        test_model.LPCReboilerHdlpSqrtTr = Var(test_model.Component, initialize=0.8)
        test_model.LPCReboilerHdlpSqrt_ai = Var(test_model.Component, test_model.Component, initialize=0.17,
                                                bounds=(0.05, 0.4))
        test_model.LPCReboilerHdlpZL = Var(initialize=0.023966, bounds=(0, 0.05))
        test_model.LPCReboilerHdlpZV = Var(initialize=0.857334, bounds=(0.2, 3))
        test_model.LPCReboilerHdlpZLBL = Var(initialize=-4.878610, bounds=(-10, -1))
        test_model.LPCReboilerHdlpZLdivBL = Var(initialize=1.30898, bounds=(0.8, 2))
        test_model.LPCReboilerHdlpZVBV = Var(initialize=-0.1732, bounds=(-0.3, 0))
        test_model.LPCReboilerHdlpZVdivBV = Var(initialize=0.053666900, bounds=(0, 0.2))
        test_model.LPCReboilerHdlpVLE_K = Var(test_model.Component, bounds=(0, 8))
        # --------ASC---------
        test_model.ASCVapLvMFlow = Var(test_model.ASCTrays, within=NonNegativeReals)
        test_model.ASCLiqLvMFlow = Var(test_model.ASCTrays, within=NonNegativeReals)
        test_model.ASCTrayTemp = Var(test_model.ASCTrays)
        test_model.ASCTrayPressure = Var(test_model.ASCTrays)
        test_model.ASCLiqLvMFrac = Var(test_model.ASCTrays, test_model.Component, bounds=(0, 1))
        test_model.ASCVapLvMFrac = Var(test_model.ASCTrays, test_model.Component, bounds=(0, 1))
        test_model.ASCAllTraysSqrtTr = Var(test_model.ASCTrays, test_model.Component, initialize=0.8)
        test_model.ASCAllTraysSqrt_ai = Var(test_model.ASCTrays, test_model.Component, test_model.Component,
                                            initialize=0.17, bounds=(0.05, 0.4))
        test_model.ASCAllTraysZL = Var(test_model.ASCTrays, initialize=0.0044, bounds=(0, 0.05))
        test_model.ASCAllTraysZV = Var(test_model.ASCTrays, initialize=0.961788, bounds=(0.2, 3))
        test_model.ASCAllTraysZLB = Var(test_model.ASCTrays, initialize=-7.427797, bounds=(-10, -4))
        test_model.ASCAllTraysBZL = Var(test_model.ASCTrays, initialize=0.625059, bounds=(0.1, 2))
        test_model.ASCAllTraysZVB = Var(test_model.ASCTrays, initialize=-0.04261712, bounds=(-0.07, -0.01))
        test_model.ASCAllTraysBZV = Var(test_model.ASCTrays, initialize=0.004027, bounds=(0.001, 0.007))
        test_model.ASCAllTraysVLE_K = Var(test_model.ASCTrays, test_model.Component, bounds=(0, 8))
        test_model.ASCAllTrayssqrtaiTr = Var(test_model.ASCTrays, test_model.Component, initialize=0.29,
                                             bounds=(0.1, 0.6))
        test_model.ASCAllTrayssqrtam = Var(test_model.ASCTrays, initialize=0.411244, bounds=(0.2, 0.7))
        # --------ASCSump---------
        test_model.ASCSumpOutMFlow = Var(within=NonNegativeReals)
        test_model.ASCSumpHldpMFrac = Var(test_model.Component, bounds=(0, 1))
        # --------ASCCond---------
        test_model.ASCCondVapMFlow = Var(within=NonNegativeReals)
        test_model.ASCCondLiqMFlow = Var(within=NonNegativeReals)
        test_model.ASCCondTemp = Var()
        test_model.ASCCondMHeat = Var(within=NonNegativeReals)
        test_model.ASCCondLiqMFrac = Var(test_model.Component, bounds=(0, 1))
        test_model.ASCCondVapMFrac = Var(test_model.Component, bounds=(0, 1))
        test_model.ASCCondSelfSqrtTr = Var(test_model.Component, initialize=0.8)
        test_model.ASCCondSelfSqrt_ai = Var(test_model.Component, test_model.Component, initialize=0.17,
                                            bounds=(0.05, 0.4))
        test_model.ASCCondSelfZL = Var(initialize=0.0044, bounds=(0, 0.05))
        test_model.ASCCondSelfZV = Var(initialize=0.961788, bounds=(0.2, 3))
        test_model.ASCCondSelfZLB = Var(initialize=-7.427797, bounds=(-10, -4))
        test_model.ASCCondSelfBZL = Var(initialize=0.625059, bounds=(0.1, 2))
        test_model.ASCCondSelfZVB = Var(initialize=-0.04261712, bounds=(-0.07, -0.01))
        test_model.ASCCondSelfBZV = Var(initialize=0.004027, bounds=(0.001, 0.007))
        test_model.ASCCondSelfVLE_K = Var(test_model.Component, bounds=(0, 8))
        test_model.ASCCondSelfsqrtaiTr = Var(test_model.Component, initialize=0.29, bounds=(0.1, 0.6))
        test_model.ASCCondSelfsqrtam = Var(initialize=0.411244, bounds=(0.2, 0.7))
        # --------HeatGan---------
        test_model.HeatGanTemp = Var()
        test_model.HeatGanPressure = Var()
        test_model.HeatGanSelfSqrtTr = Var(test_model.Component, initialize=0.8)
        test_model.HeatGanSelfSqrt_ai = Var(test_model.Component, test_model.Component, initialize=0.17,
                                            bounds=(0.05, 0.4))
        test_model.HeatGanSelfZV = Var(initialize=0.866923293, bounds=(0.2, 3))
        test_model.HeatGanSelfZVdivBV = Var(initialize=0.050441702, bounds=(0, 0.2))
        # --------RecoGan---------
        test_model.RecoGanMHeat = Var(within=NonNegativeReals)
        test_model.RecoGanSelfSqrtTr = Var(test_model.Component, initialize=0.8)
        test_model.RecoGanSelfSqrt_ai = Var(test_model.Component, test_model.Component, initialize=0.17,
                                            bounds=(0.05, 0.4))
        test_model.RecoGanSelfZV = Var(initialize=0.866923293, bounds=(0.2, 3))
        test_model.RecoGanSelfZVdivBV = Var(initialize=0.050441702, bounds=(0, 0.2))
        # --------RecoGar---------
        test_model.RecoGarMHeat = Var(within=NonNegativeReals)
        test_model.RecoGarSelfSqrtTr = Var(test_model.Component, initialize=0.8)
        test_model.RecoGarSelfSqrt_ai = Var(test_model.Component, test_model.Component, initialize=0.17,
                                            bounds=(0.05, 0.4))
        test_model.RecoGarSelfZV = Var(initialize=0.866923293, bounds=(0.2, 3))
        test_model.RecoGarSelfZVdivBV = Var(initialize=0.050441702, bounds=(0, 0.2))
        # --------RecoWn---------
        test_model.RecoWnMHeat = Var(within=NonNegativeReals)
        test_model.RecoWnSelfSqrtTr = Var(test_model.Component, initialize=0.8)
        test_model.RecoWnSelfSqrt_ai = Var(test_model.Component, test_model.Component, initialize=0.17,
                                           bounds=(0.05, 0.4))
        test_model.RecoWnSelfZV = Var(initialize=0.866923293, bounds=(0.2, 3))
        test_model.RecoWnSelfZVdivBV = Var(initialize=0.050441702, bounds=(0, 0.2))
        # --------OxSplitter---------
        test_model.OxSplitterOutAMFlow = Var(within=NonNegativeReals)
        test_model.OxSplitterOutBMFlow = Var(within=NonNegativeReals)
        # --------RecoGox---------
        test_model.RecoGoxMHeat = Var(within=NonNegativeReals)
        test_model.RecoGoxSelfSqrtTr = Var(test_model.Component, initialize=0.8)
        test_model.RecoGoxSelfSqrt_ai = Var(test_model.Component, test_model.Component, initialize=0.17,
                                            bounds=(0.05, 0.4))
        test_model.RecoGoxSelfZV = Var(initialize=0.866923293, bounds=(0.2, 3))
        test_model.RecoGoxSelfZVdivBV = Var(initialize=0.050441702, bounds=(0, 0.2))

        # ===================================
        #
        #         Expressions
        #
        # ===================================
        # --------TaCooler Enthalpy---------
        def TaCoolerSelfTr_Calc(m, comp):
            return (m.TaCoolerTemp + 273.15) / m.PR_CrtcTemp[comp]

        test_model.TaCoolerSelfTr = Expression(test_model.Component, rule=TaCoolerSelfTr_Calc)

        def TaCoolerSelfalpha_Calc(m, comp):
            return (1 + m.PR_m[comp] * (1 - m.TaCoolerSelfSqrtTr[comp])) ** 2

        test_model.TaCoolerSelfalpha = Expression(test_model.Component, rule=TaCoolerSelfalpha_Calc)

        def TaCoolerSelfai_Calc(m, comp):
            return 0.45724 * m.TaCoolerSelfalpha[comp] * (m.PR_RGas ** 2) * (m.PR_CrtcTemp[comp] ** 2) / \
                m.PR_CrtcPressure[
                    comp]

        test_model.TaCoolerSelfai = Expression(test_model.Component, rule=TaCoolerSelfai_Calc)

        def TaCoolerSelfaij_Calc(m, comp1, comp2):
            return m.TaCoolerSelfSqrt_ai[comp1, comp2] * (1 - m.PR_Kij[comp1, comp2])

        test_model.TaCoolerSelfaij = Expression(test_model.Component, test_model.Component, rule=TaCoolerSelfaij_Calc)

        def TaCoolerSelfaL_Calc(m):
            return sum(
                [sum([m.TaCoolerSelfaij[c1, c2] * m.TaCoolerLiqMFrac[c1] * m.TaCoolerLiqMFrac[c2] for c1 in
                      m.Component])
                 for c2 in m.Component])

        test_model.TaCoolerSelfaL = Expression(rule=TaCoolerSelfaL_Calc)

        def TaCoolerSelfbL_Calc(m):
            return sum([m.PR_bi[c] * m.TaCoolerLiqMFrac[c] for c in m.Component])

        test_model.TaCoolerSelfbL = Expression(rule=TaCoolerSelfbL_Calc)

        def TaCoolerSelfAL_Calc(m):
            return m.TaCoolerSelfaL * m.TaCoolerPressure * 1000 / ((m.PR_RGas * (m.TaCoolerTemp + 273.15)) ** 2)

        test_model.TaCoolerSelfAL = Expression(rule=TaCoolerSelfAL_Calc)

        def TaCoolerSelfBL_Calc(m):
            return m.TaCoolerSelfbL * m.TaCoolerPressure * 1000 / (m.PR_RGas * (m.TaCoolerTemp + 273.15))

        test_model.TaCoolerSelfBL = Expression(rule=TaCoolerSelfBL_Calc)

        def TaCoolerSelfa1_Calc(m):
            return m.TaCoolerSelfBL - 1

        test_model.TaCoolerSelfa1 = Expression(rule=TaCoolerSelfa1_Calc)

        def TaCoolerSelfa2_Calc(m):
            return m.TaCoolerSelfAL - 3 * m.TaCoolerSelfBL ** 2 - 2 * m.TaCoolerSelfBL

        test_model.TaCoolerSelfa2 = Expression(rule=TaCoolerSelfa2_Calc)

        def TaCoolerSelfa3_Calc(m):
            return m.TaCoolerSelfBL ** 2 + m.TaCoolerSelfBL ** 3 - m.TaCoolerSelfAL * m.TaCoolerSelfBL

        test_model.TaCoolerSelfa3 = Expression(rule=TaCoolerSelfa3_Calc)

        def TaCoolerSelfaV_Calc(m):
            return sum(
                [sum([m.TaCoolerSelfaij[c1, c2] * m.TaCoolerVapMFrac[c1] * m.TaCoolerVapMFrac[c2] for c1 in
                      m.Component])
                 for c2 in m.Component])

        test_model.TaCoolerSelfaV = Expression(rule=TaCoolerSelfaV_Calc)

        def TaCoolerSelfbV_Calc(m):
            return sum([m.PR_bi[c] * m.TaCoolerVapMFrac[c] for c in m.Component])

        test_model.TaCoolerSelfbV = Expression(rule=TaCoolerSelfbV_Calc)

        def TaCoolerSelfAV_Calc(m):
            return m.TaCoolerSelfaV * m.TaCoolerPressure * 1000 / ((m.PR_RGas * (m.TaCoolerTemp + 273.15)) ** 2)

        test_model.TaCoolerSelfAV = Expression(rule=TaCoolerSelfAV_Calc)

        def TaCoolerSelfBV_Calc(m):
            return m.TaCoolerSelfbV * m.TaCoolerPressure * 1000 / (m.PR_RGas * (m.TaCoolerTemp + 273.15))

        test_model.TaCoolerSelfBV = Expression(rule=TaCoolerSelfBV_Calc)

        def TaCoolerSelfb1_Calc(m):
            return m.TaCoolerSelfBV - 1

        test_model.TaCoolerSelfb1 = Expression(rule=TaCoolerSelfb1_Calc)

        def TaCoolerSelfb2_Calc(m):
            return m.TaCoolerSelfAV - 3 * m.TaCoolerSelfBV ** 2 - 2 * m.TaCoolerSelfBV

        test_model.TaCoolerSelfb2 = Expression(rule=TaCoolerSelfb2_Calc)

        def TaCoolerSelfb3_Calc(m):
            return m.TaCoolerSelfBV ** 2 + m.TaCoolerSelfBV ** 3 - m.TaCoolerSelfAV * m.TaCoolerSelfBV

        test_model.TaCoolerSelfb3 = Expression(rule=TaCoolerSelfb3_Calc)

        def TaCoolerSelfSL_Calc(m, comp):
            return sum([m.TaCoolerSelfaij[c, comp] * m.TaCoolerLiqMFrac[c] for c in m.Component])

        test_model.TaCoolerSelfSL = Expression(test_model.Component, rule=TaCoolerSelfSL_Calc)

        def TaCoolerSelfSV_Calc(m, comp):
            return sum([m.TaCoolerSelfaij[c, comp] * m.TaCoolerVapMFrac[c] for c in m.Component])

        test_model.TaCoolerSelfSV = Expression(test_model.Component, rule=TaCoolerSelfSV_Calc)

        def TaCoolerSelfPhiL_Calc(m, comp):
            return exp(m.PR_bi[comp] / m.TaCoolerSelfbL * (m.TaCoolerSelfZL - 1) - \
                       m.TaCoolerSelfZLBL - m.TaCoolerSelfAL / m.TaCoolerSelfBL / 2 / sqrt(2) * (2 * \
                                                                                                 m.TaCoolerSelfSL[
                                                                                                     comp] / m.TaCoolerSelfaL -
                                                                                                 m.PR_bi[
                                                                                                     comp] / m.TaCoolerSelfbL) * \
                       m.TaCoolerSelfZLdivBL)

        test_model.TaCoolerSelfPhiL = Expression(test_model.Component, rule=TaCoolerSelfPhiL_Calc)

        def TaCoolerSelfPhiV_Calc(m, comp):
            return exp(m.PR_bi[comp] / m.TaCoolerSelfbV * (m.TaCoolerSelfZV - 1) - \
                       m.TaCoolerSelfZVBV - m.TaCoolerSelfAV / m.TaCoolerSelfBV / 2 / sqrt(2) * (2 * \
                                                                                                 m.TaCoolerSelfSV[
                                                                                                     comp] / m.TaCoolerSelfaV -
                                                                                                 m.PR_bi[
                                                                                                     comp] / m.TaCoolerSelfbV) * \
                       m.TaCoolerSelfZVdivBV)

        test_model.TaCoolerSelfPhiV = Expression(test_model.Component, rule=TaCoolerSelfPhiV_Calc)

        def TaCoolerSelfH0_Calc(m, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.TaCoolerTemp + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (
                        m.PR_DHFORM[comp] + m.PR_CPIGDP1[comp] * ((m.TaCoolerTemp + 273.15) - m.PR_RefTemp[comp]) + 2 *
                        m.PR_CPIGDP2[comp] * \
                        m.PR_CPIGDP3[comp] * (-exp(2 / (m.TaCoolerTemp + 273.15) * m.PR_CPIGDP3[comp]) + \
                                              exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp])) / (
                                exp(2 / (m.TaCoolerTemp + 273.15) * m.PR_CPIGDP3[comp]) - 1) / (
                                exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp]) - 1) - \
                        2 * m.PR_CPIGDP4[comp] * m.PR_CPIGDP5[comp] * (
                                exp(2 / (m.TaCoolerTemp + 273.15) * m.PR_CPIGDP5[comp]) - exp(
                            2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp])) / \
                        (exp(2 / (m.TaCoolerTemp + 273.15) * m.PR_CPIGDP5[comp]) + 1) / (
                                exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp]) + 1))

        test_model.TaCoolerSelfH0 = Expression(test_model.Component, rule=TaCoolerSelfH0_Calc)

        def TaCoolerSelfM_Calc(m, comp1, comp2):
            return m.PR_m[comp1] * m.TaCoolerSelfSqrtTr[comp1] / (2 * sqrt(m.TaCoolerSelfalpha[comp1]))

        test_model.TaCoolerSelfM = Expression(test_model.Component, test_model.Component,
                                              rule=TaCoolerSelfM_Calc)  # Calculate the enthalpy of liquid phase mixture %unit: J/mol

        def TaCoolerSelfHV_Calc(m):
            return sum([m.TaCoolerVapMFrac[c] * m.TaCoolerSelfH0[c] for c in m.Component]) + \
                m.PR_RGas * (m.TaCoolerTemp + 273.15) * (
                        (m.TaCoolerSelfZV - 1) - 1 / (2 ** 1.5) * (1 / m.TaCoolerSelfBV) * m.TaCoolerSelfZVdivBV * \
                        sum([sum([m.TaCoolerVapMFrac[c1] * m.TaCoolerVapMFrac[c2] * m.TaCoolerSelfAV * (
                                1 + m.TaCoolerSelfM[c2, c1] + \
                                m.TaCoolerSelfM[c1, c2]) for c1 in m.Component]) for c2 in m.Component]))

        test_model.TaCoolerVapMEtlp = Expression(rule=TaCoolerSelfHV_Calc)

        # Calculate the enthalpy of liquid phase mixture %unit: J/mol
        def TaCoolerSelfHL_Calc(m):
            return sum([m.TaCoolerLiqMFrac[c] * m.TaCoolerSelfH0[c] for c in m.Component]) + \
                m.PR_RGas * (m.TaCoolerTemp + 273.15) * (
                        (m.TaCoolerSelfZL - 1) - 1 / (2 ** 1.5) * (1 / m.TaCoolerSelfBL) * m.TaCoolerSelfZLdivBL * \
                        sum([sum([m.TaCoolerLiqMFrac[c1] * m.TaCoolerLiqMFrac[c2] * m.TaCoolerSelfAL * (
                                1 + m.TaCoolerSelfM[c2, c1] + \
                                m.TaCoolerSelfM[c1, c2]) for c1 in m.Component]) for c2 in m.Component]))

        test_model.TaCoolerLiqMEtlp = Expression(rule=TaCoolerSelfHL_Calc)

        def TaCoolerMEtlp(m):
            return (
                    m.TaCoolerLiqMEtlp * m.TaCoolerLiqMFlow + m.TaCoolerVapMEtlp * m.TaCoolerVapMFlow) / m.FeedSplitterOut1MFlow

        test_model.TaCoolerMEtlp = Expression(rule=TaCoolerMEtlp)

        # --------HpaCooler Enthalpy---------
        def HpaCoolerSelfTr_Calc(m, comp):
            return (m.HpaCoolerTemp + 273.15) / m.PR_CrtcTemp[comp]

        test_model.HpaCoolerSelfTr = Expression(test_model.Component, rule=HpaCoolerSelfTr_Calc)

        def HpaCoolerSelfalpha_Calc(m, comp):
            return (1 + m.PR_m[comp] * (1 - m.HpaCoolerSelfSqrtTr[comp])) ** 2

        test_model.HpaCoolerSelfalpha = Expression(test_model.Component, rule=HpaCoolerSelfalpha_Calc)

        def HpaCoolerSelfai_Calc(m, comp):
            return 0.45724 * m.HpaCoolerSelfalpha[comp] * (m.PR_RGas ** 2) * (m.PR_CrtcTemp[comp] ** 2) / \
                m.PR_CrtcPressure[
                    comp]

        test_model.HpaCoolerSelfai = Expression(test_model.Component, rule=HpaCoolerSelfai_Calc)

        def HpaCoolerSelfaij_Calc(m, comp1, comp2):
            return m.HpaCoolerSelfSqrt_ai[comp1, comp2] * (1 - m.PR_Kij[comp1, comp2])

        test_model.HpaCoolerSelfaij = Expression(test_model.Component, test_model.Component, rule=HpaCoolerSelfaij_Calc)

        def HpaCoolerSelfaL_Calc(m):
            return sum(
                [sum([m.HpaCoolerSelfaij[c1, c2] * m.HpaCoolerLiqMFrac[c1] * m.HpaCoolerLiqMFrac[c2] for c1 in
                      m.Component])
                 for c2 in m.Component])

        test_model.HpaCoolerSelfaL = Expression(rule=HpaCoolerSelfaL_Calc)

        def HpaCoolerSelfbL_Calc(m):
            return sum([m.PR_bi[c] * m.HpaCoolerLiqMFrac[c] for c in m.Component])

        test_model.HpaCoolerSelfbL = Expression(rule=HpaCoolerSelfbL_Calc)

        def HpaCoolerSelfAL_Calc(m):
            return m.HpaCoolerSelfaL * m.HpaCoolerPressure * 1000 / ((m.PR_RGas * (m.HpaCoolerTemp + 273.15)) ** 2)

        test_model.HpaCoolerSelfAL = Expression(rule=HpaCoolerSelfAL_Calc)

        def HpaCoolerSelfBL_Calc(m):
            return m.HpaCoolerSelfbL * m.HpaCoolerPressure * 1000 / (m.PR_RGas * (m.HpaCoolerTemp + 273.15))

        test_model.HpaCoolerSelfBL = Expression(rule=HpaCoolerSelfBL_Calc)

        def HpaCoolerSelfa1_Calc(m):
            return m.HpaCoolerSelfBL - 1

        test_model.HpaCoolerSelfa1 = Expression(rule=HpaCoolerSelfa1_Calc)

        def HpaCoolerSelfa2_Calc(m):
            return m.HpaCoolerSelfAL - 3 * m.HpaCoolerSelfBL ** 2 - 2 * m.HpaCoolerSelfBL

        test_model.HpaCoolerSelfa2 = Expression(rule=HpaCoolerSelfa2_Calc)

        def HpaCoolerSelfa3_Calc(m):
            return m.HpaCoolerSelfBL ** 2 + m.HpaCoolerSelfBL ** 3 - m.HpaCoolerSelfAL * m.HpaCoolerSelfBL

        test_model.HpaCoolerSelfa3 = Expression(rule=HpaCoolerSelfa3_Calc)

        def HpaCoolerSelfaV_Calc(m):
            return sum(
                [sum([m.HpaCoolerSelfaij[c1, c2] * m.HpaCoolerVapMFrac[c1] * m.HpaCoolerVapMFrac[c2] for c1 in
                      m.Component])
                 for c2 in m.Component])

        test_model.HpaCoolerSelfaV = Expression(rule=HpaCoolerSelfaV_Calc)

        def HpaCoolerSelfbV_Calc(m):
            return sum([m.PR_bi[c] * m.HpaCoolerVapMFrac[c] for c in m.Component])

        test_model.HpaCoolerSelfbV = Expression(rule=HpaCoolerSelfbV_Calc)

        def HpaCoolerSelfAV_Calc(m):
            return m.HpaCoolerSelfaV * m.HpaCoolerPressure * 1000 / ((m.PR_RGas * (m.HpaCoolerTemp + 273.15)) ** 2)

        test_model.HpaCoolerSelfAV = Expression(rule=HpaCoolerSelfAV_Calc)

        def HpaCoolerSelfBV_Calc(m):
            return m.HpaCoolerSelfbV * m.HpaCoolerPressure * 1000 / (m.PR_RGas * (m.HpaCoolerTemp + 273.15))

        test_model.HpaCoolerSelfBV = Expression(rule=HpaCoolerSelfBV_Calc)

        def HpaCoolerSelfb1_Calc(m):
            return m.HpaCoolerSelfBV - 1

        test_model.HpaCoolerSelfb1 = Expression(rule=HpaCoolerSelfb1_Calc)

        def HpaCoolerSelfb2_Calc(m):
            return m.HpaCoolerSelfAV - 3 * m.HpaCoolerSelfBV ** 2 - 2 * m.HpaCoolerSelfBV

        test_model.HpaCoolerSelfb2 = Expression(rule=HpaCoolerSelfb2_Calc)

        def HpaCoolerSelfb3_Calc(m):
            return m.HpaCoolerSelfBV ** 2 + m.HpaCoolerSelfBV ** 3 - m.HpaCoolerSelfAV * m.HpaCoolerSelfBV

        test_model.HpaCoolerSelfb3 = Expression(rule=HpaCoolerSelfb3_Calc)

        def HpaCoolerSelfSL_Calc(m, comp):
            return sum([m.HpaCoolerSelfaij[c, comp] * m.HpaCoolerLiqMFrac[c] for c in m.Component])

        test_model.HpaCoolerSelfSL = Expression(test_model.Component, rule=HpaCoolerSelfSL_Calc)

        def HpaCoolerSelfSV_Calc(m, comp):
            return sum([m.HpaCoolerSelfaij[c, comp] * m.HpaCoolerVapMFrac[c] for c in m.Component])

        test_model.HpaCoolerSelfSV = Expression(test_model.Component, rule=HpaCoolerSelfSV_Calc)

        def HpaCoolerSelfPhiL_Calc(m, comp):
            return exp(m.PR_bi[comp] / m.HpaCoolerSelfbL * (m.HpaCoolerSelfZL - 1) - \
                       m.HpaCoolerSelfZLBL - m.HpaCoolerSelfAL / m.HpaCoolerSelfBL / 2 / sqrt(2) * (2 * \
                                                                                                    m.HpaCoolerSelfSL[
                                                                                                        comp] / m.HpaCoolerSelfaL -
                                                                                                    m.PR_bi[
                                                                                                        comp] / m.HpaCoolerSelfbL) * \
                       m.HpaCoolerSelfZLdivBL)

        test_model.HpaCoolerSelfPhiL = Expression(test_model.Component, rule=HpaCoolerSelfPhiL_Calc)

        def HpaCoolerSelfPhiV_Calc(m, comp):
            return exp(m.PR_bi[comp] / m.HpaCoolerSelfbV * (m.HpaCoolerSelfZV - 1) - \
                       m.HpaCoolerSelfZVBV - m.HpaCoolerSelfAV / m.HpaCoolerSelfBV / 2 / sqrt(2) * (2 * \
                                                                                                    m.HpaCoolerSelfSV[
                                                                                                        comp] / m.HpaCoolerSelfaV -
                                                                                                    m.PR_bi[
                                                                                                        comp] / m.HpaCoolerSelfbV) * \
                       m.HpaCoolerSelfZVdivBV)

        test_model.HpaCoolerSelfPhiV = Expression(test_model.Component, rule=HpaCoolerSelfPhiV_Calc)

        def HpaCoolerSelfH0_Calc(m, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.HpaCoolerTemp + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (
                        m.PR_DHFORM[comp] + m.PR_CPIGDP1[comp] * ((m.HpaCoolerTemp + 273.15) - m.PR_RefTemp[comp]) + 2 *
                        m.PR_CPIGDP2[comp] * \
                        m.PR_CPIGDP3[comp] * (-exp(2 / (m.HpaCoolerTemp + 273.15) * m.PR_CPIGDP3[comp]) + \
                                              exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp])) / (
                                exp(2 / (m.HpaCoolerTemp + 273.15) * m.PR_CPIGDP3[comp]) - 1) / (
                                exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp]) - 1) - \
                        2 * m.PR_CPIGDP4[comp] * m.PR_CPIGDP5[comp] * (
                                exp(2 / (m.HpaCoolerTemp + 273.15) * m.PR_CPIGDP5[comp]) - exp(
                            2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp])) / \
                        (exp(2 / (m.HpaCoolerTemp + 273.15) * m.PR_CPIGDP5[comp]) + 1) / (
                                exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp]) + 1))

        test_model.HpaCoolerSelfH0 = Expression(test_model.Component, rule=HpaCoolerSelfH0_Calc)

        def HpaCoolerSelfM_Calc(m, comp1, comp2):
            return m.PR_m[comp1] * m.HpaCoolerSelfSqrtTr[comp1] / (2 * sqrt(m.HpaCoolerSelfalpha[comp1]))

        test_model.HpaCoolerSelfM = Expression(test_model.Component, test_model.Component,
                                               rule=HpaCoolerSelfM_Calc)  # Calculate the enthalpy of liquid phase mixture %unit: J/mol

        def HpaCoolerSelfHV_Calc(m):
            return sum([m.HpaCoolerVapMFrac[c] * m.HpaCoolerSelfH0[c] for c in m.Component]) + \
                m.PR_RGas * (m.HpaCoolerTemp + 273.15) * (
                        (m.HpaCoolerSelfZV - 1) - 1 / (2 ** 1.5) * (1 / m.HpaCoolerSelfBV) * m.HpaCoolerSelfZVdivBV * \
                        sum([sum([m.HpaCoolerVapMFrac[c1] * m.HpaCoolerVapMFrac[c2] * m.HpaCoolerSelfAV * (
                                1 + m.HpaCoolerSelfM[c2, c1] + \
                                m.HpaCoolerSelfM[c1, c2]) for c1 in m.Component]) for c2 in m.Component]))

        test_model.HpaCoolerVapMEtlp = Expression(rule=HpaCoolerSelfHV_Calc)

        # Calculate the enthalpy of liquid phase mixture %unit: J/mol
        def HpaCoolerSelfHL_Calc(m):
            return sum([m.HpaCoolerLiqMFrac[c] * m.HpaCoolerSelfH0[c] for c in m.Component]) + \
                m.PR_RGas * (m.HpaCoolerTemp + 273.15) * (
                        (m.HpaCoolerSelfZL - 1) - 1 / (2 ** 1.5) * (1 / m.HpaCoolerSelfBL) * m.HpaCoolerSelfZLdivBL * \
                        sum([sum([m.HpaCoolerLiqMFrac[c1] * m.HpaCoolerLiqMFrac[c2] * m.HpaCoolerSelfAL * (
                                1 + m.HpaCoolerSelfM[c2, c1] + \
                                m.HpaCoolerSelfM[c1, c2]) for c1 in m.Component]) for c2 in m.Component]))

        test_model.HpaCoolerLiqMEtlp = Expression(rule=HpaCoolerSelfHL_Calc)

        def HpaCoolerMEtlp(m):
            return (
                    m.HpaCoolerLiqMEtlp * m.HpaCoolerLiqMFlow + m.HpaCoolerVapMEtlp * m.HpaCoolerVapMFlow) / m.FeedSplitterOut3MFlow

        test_model.HpaCoolerMEtlp = Expression(rule=HpaCoolerMEtlp)

        # --------MaCooler Enthalpy---------
        def MaCoolerSelfTr_Calc(m, comp):
            return (m.MaCoolerTemp + 273.15) / m.PR_CrtcTemp[comp]

        test_model.MaCoolerSelfTr = Expression(test_model.Component, rule=MaCoolerSelfTr_Calc)

        def MaCoolerSelfalpha_Calc(m, comp):
            return (1 + m.PR_m[comp] * (1 - m.MaCoolerSelfSqrtTr[comp])) ** 2

        test_model.MaCoolerSelfalpha = Expression(test_model.Component, rule=MaCoolerSelfalpha_Calc)

        def MaCoolerSelfai_Calc(m, comp):
            return 0.45724 * m.MaCoolerSelfalpha[comp] * (m.PR_RGas ** 2) * (m.PR_CrtcTemp[comp] ** 2) / \
                m.PR_CrtcPressure[
                    comp]

        test_model.MaCoolerSelfai = Expression(test_model.Component, rule=MaCoolerSelfai_Calc)

        def MaCoolerSelfaij_Calc(m, comp1, comp2):
            return m.MaCoolerSelfSqrt_ai[comp1, comp2] * (1 - m.PR_Kij[comp1, comp2])

        test_model.MaCoolerSelfaij = Expression(test_model.Component, test_model.Component, rule=MaCoolerSelfaij_Calc)

        def MaCoolerSelfaV_Calc(m):
            return sum(
                [sum([m.MaCoolerSelfaij[c1, c2] * m.FeedMFrac[c1] * m.FeedMFrac[c2] for c1 in m.Component]) for c2 in
                 m.Component])

        test_model.MaCoolerSelfaV = Expression(rule=MaCoolerSelfaV_Calc)

        def MaCoolerSelfbV_Calc(m):
            return sum([m.PR_bi[c] * m.FeedMFrac[c] for c in m.Component])

        test_model.MaCoolerSelfbV = Expression(rule=MaCoolerSelfbV_Calc)

        def MaCoolerSelfAV_Calc(m):
            return m.MaCoolerSelfaV * m.MaCoolerPressure * 1000 / ((m.PR_RGas * (m.MaCoolerTemp + 273.15)) ** 2)

        test_model.MaCoolerSelfAV = Expression(rule=MaCoolerSelfAV_Calc)

        def MaCoolerSelfBV_Calc(m):
            return m.MaCoolerSelfbV * m.MaCoolerPressure * 1000 / (m.PR_RGas * (m.MaCoolerTemp + 273.15))

        test_model.MaCoolerSelfBV = Expression(rule=MaCoolerSelfBV_Calc)

        def MaCoolerSelfb1_Calc(m):
            return m.MaCoolerSelfBV - 1

        test_model.MaCoolerSelfb1 = Expression(rule=MaCoolerSelfb1_Calc)

        def MaCoolerSelfb2_Calc(m):
            return m.MaCoolerSelfAV - 3 * m.MaCoolerSelfBV ** 2 - 2 * m.MaCoolerSelfBV

        test_model.MaCoolerSelfb2 = Expression(rule=MaCoolerSelfb2_Calc)

        def MaCoolerSelfb3_Calc(m):
            return m.MaCoolerSelfBV ** 2 + m.MaCoolerSelfBV ** 3 - m.MaCoolerSelfAV * m.MaCoolerSelfBV

        test_model.MaCoolerSelfb3 = Expression(rule=MaCoolerSelfb3_Calc)

        def MaCoolerSelfH0_Calc(m, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.MaCoolerTemp + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (
                        m.PR_DHFORM[comp] + m.PR_CPIGDP1[comp] * ((m.MaCoolerTemp + 273.15) - m.PR_RefTemp[comp]) + 2 *
                        m.PR_CPIGDP2[comp] * \
                        m.PR_CPIGDP3[comp] * (-exp(2 / (m.MaCoolerTemp + 273.15) * m.PR_CPIGDP3[comp]) + \
                                              exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp])) / (
                                exp(2 / (m.MaCoolerTemp + 273.15) * m.PR_CPIGDP3[comp]) - 1) / (
                                exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp]) - 1) - \
                        2 * m.PR_CPIGDP4[comp] * m.PR_CPIGDP5[comp] * (
                                exp(2 / (m.MaCoolerTemp + 273.15) * m.PR_CPIGDP5[comp]) - exp(
                            2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp])) / \
                        (exp(2 / (m.MaCoolerTemp + 273.15) * m.PR_CPIGDP5[comp]) + 1) / (
                                exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp]) + 1))

        test_model.MaCoolerSelfH0 = Expression(test_model.Component, rule=MaCoolerSelfH0_Calc)

        def MaCoolerSelfM_Calc(m, comp1, comp2):
            return m.PR_m[comp1] * m.MaCoolerSelfSqrtTr[comp1] / (2 * sqrt(m.MaCoolerSelfalpha[comp1]))

        test_model.MaCoolerSelfM = Expression(test_model.Component, test_model.Component,
                                              rule=MaCoolerSelfM_Calc)  # Calculate the enthalpy of liquid phase mixture %unit: J/mol

        def MaCoolerSelfHV_Calc(m):
            return sum([m.FeedMFrac[c] * m.MaCoolerSelfH0[c] for c in m.Component]) + \
                m.PR_RGas * (m.MaCoolerTemp + 273.15) * (
                        (m.MaCoolerSelfZV - 1) - 1 / (2 ** 1.5) * (1 / m.MaCoolerSelfBV) * m.MaCoolerSelfZVdivBV * \
                        sum([sum([m.FeedMFrac[c1] * m.FeedMFrac[c2] * m.MaCoolerSelfAV * (1 + m.MaCoolerSelfM[c2, c1] + \
                                                                                          m.MaCoolerSelfM[c1, c2]) for
                                  c1 in m.Component]) for c2 in m.Component]))

        test_model.MaCoolerMEtlp = Expression(rule=MaCoolerSelfHV_Calc)

        # --------HPC Component Flowrate---------
        def HPCVapLvMCompFlow(m, tray, comp):
            return m.HPCVapLvMFrac[tray, comp] * m.HPCVapLvMFlow[tray]

        test_model.HPCVapLvMCompFlow = Expression(test_model.HPCTrays, test_model.Component, rule=HPCVapLvMCompFlow)

        def HPCLiqLvMCompFlow(m, tray, comp):
            return m.HPCLiqLvMFrac[tray, comp] * m.HPCLiqLvMFlow[tray]

        test_model.HPCLiqLvMCompFlow = Expression(test_model.HPCTrays, test_model.Component, rule=HPCLiqLvMCompFlow)

        # --------HPC Thermo and Enthalpy---------
        def HPCAllTraysTr_Calc(m, tray, comp):
            return (m.HPCTrayTemp[tray] + 273.15) / m.PR_CrtcTemp[comp]

        test_model.HPCAllTraysTr = Expression(test_model.HPCTrays, test_model.Component, rule=HPCAllTraysTr_Calc)

        def HPCAllTraysalpha_Calc(m, tray, comp):
            return (1 + m.PR_m[comp] * (1 - m.HPCAllTraysSqrtTr[tray, comp])) ** 2

        test_model.HPCAllTraysalpha = Expression(test_model.HPCTrays, test_model.Component, rule=HPCAllTraysalpha_Calc)

        def HPCAllTraysai_Calc(m, tray, comp):
            return 0.45724 * m.HPCAllTraysalpha[tray, comp] * (m.PR_RGas ** 2) * (m.PR_CrtcTemp[comp] ** 2) / \
                m.PR_CrtcPressure[comp]

        test_model.HPCAllTraysai = Expression(test_model.HPCTrays, test_model.Component, rule=HPCAllTraysai_Calc)

        def HPCAllTraysaij_Calc(m, tray, comp1, comp2):
            return m.HPCAllTraysSqrt_ai[tray, comp1, comp2] * (1 - m.PR_Kij[comp1, comp2])

        test_model.HPCAllTraysaij = Expression(test_model.HPCTrays, test_model.Component, test_model.Component,
                                               rule=HPCAllTraysaij_Calc)

        def HPCAllTraysaL_Calc(m, tray):
            return sum(
                [sum([m.HPCAllTraysaij[tray, c1, c2] * m.HPCLiqLvMFrac[tray, c1] * m.HPCLiqLvMFrac[tray, c2] for c1
                      in m.Component]) for c2 in m.Component])

        test_model.HPCAllTraysaL = Expression(test_model.HPCTrays, rule=HPCAllTraysaL_Calc)

        def HPCAllTraysbL_Calc(m, tray):
            return sum([m.PR_bi[c] * m.HPCLiqLvMFrac[tray, c] for c in m.Component])

        test_model.HPCAllTraysbL = Expression(test_model.HPCTrays, rule=HPCAllTraysbL_Calc)

        def HPCAllTraysAL_Calc(m, tray):
            return m.HPCAllTraysaL[tray] * m.HPCTrayPressure[tray] * 1000 / (
                    (m.PR_RGas * (m.HPCTrayTemp[tray] + 273.15)) ** 2)

        test_model.HPCAllTraysAL = Expression(test_model.HPCTrays, rule=HPCAllTraysAL_Calc)

        def HPCAllTraysBL_Calc(m, tray):
            return m.HPCAllTraysbL[tray] * m.HPCTrayPressure[tray] * 1000 / (m.PR_RGas * (m.HPCTrayTemp[tray] + 273.15))

        test_model.HPCAllTraysBL = Expression(test_model.HPCTrays, rule=HPCAllTraysBL_Calc)

        def HPCAllTraysa1_Calc(m, tray):
            return m.HPCAllTraysBL[tray] - 1

        test_model.HPCAllTraysa1 = Expression(test_model.HPCTrays, rule=HPCAllTraysa1_Calc)

        def HPCAllTraysa2_Calc(m, tray):
            return m.HPCAllTraysAL[tray] - 3 * m.HPCAllTraysBL[tray] ** 2 - 2 * m.HPCAllTraysBL[tray]

        test_model.HPCAllTraysa2 = Expression(test_model.HPCTrays, rule=HPCAllTraysa2_Calc)

        def HPCAllTraysa3_Calc(m, tray):
            return m.HPCAllTraysBL[tray] ** 2 + m.HPCAllTraysBL[tray] ** 3 - m.HPCAllTraysAL[tray] * m.HPCAllTraysBL[
                tray]

        test_model.HPCAllTraysa3 = Expression(test_model.HPCTrays, rule=HPCAllTraysa3_Calc)

        def HPCAllTraysaV_Calc(m, tray):
            return sum(
                [sum([m.HPCAllTraysaij[tray, c1, c2] * m.HPCVapLvMFrac[tray, c1] * m.HPCVapLvMFrac[tray, c2] for c1
                      in m.Component]) for c2 in m.Component])

        test_model.HPCAllTraysaV = Expression(test_model.HPCTrays, rule=HPCAllTraysaV_Calc)

        def HPCAllTraysbV_Calc(m, tray):
            return sum([m.PR_bi[c] * m.HPCVapLvMFrac[tray, c] for c in m.Component])

        test_model.HPCAllTraysbV = Expression(test_model.HPCTrays, rule=HPCAllTraysbV_Calc)

        def HPCAllTraysAV_Calc(m, tray):
            return m.HPCAllTraysaV[tray] * m.HPCTrayPressure[tray] * 1000 / (
                    (m.PR_RGas * (m.HPCTrayTemp[tray] + 273.15)) ** 2)

        test_model.HPCAllTraysAV = Expression(test_model.HPCTrays, rule=HPCAllTraysAV_Calc)

        def HPCAllTraysBV_Calc(m, tray):
            return m.HPCAllTraysbV[tray] * m.HPCTrayPressure[tray] * 1000 / (m.PR_RGas * (m.HPCTrayTemp[tray] + 273.15))

        test_model.HPCAllTraysBV = Expression(test_model.HPCTrays, rule=HPCAllTraysBV_Calc)

        def HPCAllTraysb1_Calc(m, tray):
            return m.HPCAllTraysBV[tray] - 1

        test_model.HPCAllTraysb1 = Expression(test_model.HPCTrays, rule=HPCAllTraysb1_Calc)

        def HPCAllTraysb2_Calc(m, tray):
            return m.HPCAllTraysAV[tray] - 3 * m.HPCAllTraysBV[tray] ** 2 - 2 * m.HPCAllTraysBV[tray]

        test_model.HPCAllTraysb2 = Expression(test_model.HPCTrays, rule=HPCAllTraysb2_Calc)

        def HPCAllTraysb3_Calc(m, tray):
            return m.HPCAllTraysBV[tray] ** 2 + m.HPCAllTraysBV[tray] ** 3 - m.HPCAllTraysAV[tray] * m.HPCAllTraysBV[
                tray]

        test_model.HPCAllTraysb3 = Expression(test_model.HPCTrays, rule=HPCAllTraysb3_Calc)

        def HPCAllTraysSL_Calc(m, tray, comp):
            return sum([m.HPCAllTraysaij[tray, c, comp] * m.HPCLiqLvMFrac[tray, c] for c in m.Component])

        test_model.HPCAllTraysSL = Expression(test_model.HPCTrays, test_model.Component, rule=HPCAllTraysSL_Calc)

        def HPCAllTraysSV_Calc(m, tray, comp):
            return sum([m.HPCAllTraysaij[tray, c, comp] * m.HPCVapLvMFrac[tray, c] for c in m.Component])

        test_model.HPCAllTraysSV = Expression(test_model.HPCTrays, test_model.Component, rule=HPCAllTraysSV_Calc)

        def HPCAllTraysPhiL_Calc(m, tray, comp):
            return exp(m.PR_bi[comp] / m.HPCAllTraysbL[tray] * (m.HPCAllTraysZL[tray] - 1) - \
                       m.HPCAllTraysZLBL[tray] - m.HPCAllTraysAL[tray] / m.HPCAllTraysBL[tray] / 2 / sqrt(2) * (2 * \
                                                                                                                m.HPCAllTraysSL[
                                                                                                                    tray, comp] /
                                                                                                                m.HPCAllTraysaL[
                                                                                                                    tray] -
                                                                                                                m.PR_bi[
                                                                                                                    comp] /
                                                                                                                m.HPCAllTraysbL[
                                                                                                                    tray]) * \
                       m.HPCAllTraysZLdivBL[tray])

        test_model.HPCAllTraysPhiL = Expression(test_model.HPCTrays, test_model.Component, rule=HPCAllTraysPhiL_Calc)

        def HPCAllTraysPhiV_Calc(m, tray, comp):
            return exp(m.PR_bi[comp] / m.HPCAllTraysbV[tray] * (m.HPCAllTraysZV[tray] - 1) - \
                       m.HPCAllTraysZVBV[tray] - m.HPCAllTraysAV[tray] / m.HPCAllTraysBV[tray] / 2 / sqrt(2) * (2 * \
                                                                                                                m.HPCAllTraysSV[
                                                                                                                    tray, comp] /
                                                                                                                m.HPCAllTraysaV[
                                                                                                                    tray] -
                                                                                                                m.PR_bi[
                                                                                                                    comp] /
                                                                                                                m.HPCAllTraysbV[
                                                                                                                    tray]) * \
                       m.HPCAllTraysZVdivBV[tray])

        test_model.HPCAllTraysPhiV = Expression(test_model.HPCTrays, test_model.Component, rule=HPCAllTraysPhiV_Calc)

        def HPCAllTraysH0_Calc(m, tray, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.HPCTrayTemp[tray] + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (m.PR_DHFORM[comp] + m.PR_CPIGDP1[comp] * (
                        (m.HPCTrayTemp[tray] + 273.15) - m.PR_RefTemp[comp]) + 2 * m.PR_CPIGDP2[comp] * \
                               m.PR_CPIGDP3[comp] * (-exp(2 / (m.HPCTrayTemp[tray] + 273.15) * m.PR_CPIGDP3[comp]) + \
                                                     exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp])) / (
                                       exp(2 / (m.HPCTrayTemp[tray] + 273.15) * m.PR_CPIGDP3[comp]) - 1) / (
                                       exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp]) - 1) - \
                               2 * m.PR_CPIGDP4[comp] * m.PR_CPIGDP5[comp] * (
                                       exp(2 / (m.HPCTrayTemp[tray] + 273.15) * m.PR_CPIGDP5[comp]) - exp(
                                   2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp])) / \
                               (exp(2 / (m.HPCTrayTemp[tray] + 273.15) * m.PR_CPIGDP5[comp]) + 1) / (
                                       exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp]) + 1))

        test_model.HPCAllTraysH0 = Expression(test_model.HPCTrays, test_model.Component, rule=HPCAllTraysH0_Calc)

        def HPCAllTraysM_Calc(m, tray, comp1, comp2):
            return m.PR_m[comp1] * m.HPCAllTraysSqrtTr[tray, comp1] / (2 * sqrt(m.HPCAllTraysalpha[tray, comp1]))

        test_model.HPCAllTraysM = Expression(test_model.HPCTrays, test_model.Component, test_model.Component,
                                             rule=HPCAllTraysM_Calc)  # Calculate the enthalpy of liquid phase mixture %unit: J/mol

        def HPCAllTraysHV_Calc(m, tray):
            return sum([m.HPCVapLvMFrac[tray, c] * m.HPCAllTraysH0[tray, c] for c in m.Component]) + \
                m.PR_RGas * (m.HPCTrayTemp[tray] + 273.15) * (
                        (m.HPCAllTraysZV[tray] - 1) - 1 / (2 ** 1.5) * (1 / m.HPCAllTraysBV[tray]) *
                        m.HPCAllTraysZVdivBV[tray] * \
                        sum([sum([m.HPCVapLvMFrac[tray, c1] * m.HPCVapLvMFrac[tray, c2] * m.HPCAllTraysAV[tray] * (
                                1 + m.HPCAllTraysM[tray, c2, c1] + \
                                m.HPCAllTraysM[tray, c1, c2]) for c1 in m.Component]) for c2 in m.Component]))

        test_model.HPCVapLvMEtlp = Expression(test_model.HPCTrays, rule=HPCAllTraysHV_Calc)

        # Calculate the enthalpy of liquid phase mixture %unit: J/mol
        def HPCAllTraysHL_Calc(m, tray):
            return sum([m.HPCLiqLvMFrac[tray, c] * m.HPCAllTraysH0[tray, c] for c in m.Component]) + \
                m.PR_RGas * (m.HPCTrayTemp[tray] + 273.15) * (
                        (m.HPCAllTraysZL[tray] - 1) - 1 / (2 ** 1.5) * (1 / m.HPCAllTraysBL[tray]) *
                        m.HPCAllTraysZLdivBL[tray] * \
                        sum([sum([m.HPCLiqLvMFrac[tray, c1] * m.HPCLiqLvMFrac[tray, c2] * m.HPCAllTraysAL[tray] * (
                                1 + m.HPCAllTraysM[tray, c2, c1] + \
                                m.HPCAllTraysM[tray, c1, c2]) for c1 in m.Component]) for c2 in m.Component]))

        test_model.HPCLiqLvMEtlp = Expression(test_model.HPCTrays, rule=HPCAllTraysHL_Calc)

        # --------HPCCond Enthalpy---------
        def HPCCondOutletTr_Calc(m, comp):
            return (m.HPCTrayTemp[41] + 273.15) / m.PR_CrtcTemp[comp]

        test_model.HPCCondOutletTr = Expression(test_model.Component, rule=HPCCondOutletTr_Calc)

        def HPCCondOutletalpha_Calc(m, comp):
            return (1 + m.PR_m[comp] * (1 - m.HPCCondOutletSqrtTr[comp])) ** 2

        test_model.HPCCondOutletalpha = Expression(test_model.Component, rule=HPCCondOutletalpha_Calc)

        def HPCCondOutletai_Calc(m, comp):
            return 0.45724 * m.HPCCondOutletalpha[comp] * (m.PR_RGas ** 2) * (m.PR_CrtcTemp[comp] ** 2) / \
                m.PR_CrtcPressure[
                    comp]

        test_model.HPCCondOutletai = Expression(test_model.Component, rule=HPCCondOutletai_Calc)

        def HPCCondOutletaij_Calc(m, comp1, comp2):
            return m.HPCCondOutletSqrt_ai[comp1, comp2] * (1 - m.PR_Kij[comp1, comp2])

        test_model.HPCCondOutletaij = Expression(test_model.Component, test_model.Component, rule=HPCCondOutletaij_Calc)

        def HPCCondOutletaL_Calc(m):
            return sum(
                [sum([m.HPCCondOutletaij[c1, c2] * m.HPCVapLvMFrac[41, c1] * m.HPCVapLvMFrac[41, c2] for c1 in
                      m.Component])
                 for c2 in m.Component])

        test_model.HPCCondOutletaL = Expression(rule=HPCCondOutletaL_Calc)

        def HPCCondOutletbL_Calc(m):
            return sum([m.PR_bi[c] * m.HPCVapLvMFrac[41, c] for c in m.Component])

        test_model.HPCCondOutletbL = Expression(rule=HPCCondOutletbL_Calc)

        def HPCCondOutletAL_Calc(m):
            return m.HPCCondOutletaL * m.HPCCondPressure * 1000 / ((m.PR_RGas * (m.HPCTrayTemp[41] + 273.15)) ** 2)

        test_model.HPCCondOutletAL = Expression(rule=HPCCondOutletAL_Calc)

        def HPCCondOutletBL_Calc(m):
            return m.HPCCondOutletbL * m.HPCCondPressure * 1000 / (m.PR_RGas * (m.HPCTrayTemp[41] + 273.15))

        test_model.HPCCondOutletBL = Expression(rule=HPCCondOutletBL_Calc)

        def HPCCondOutleta1_Calc(m):
            return m.HPCCondOutletBL - 1

        test_model.HPCCondOutleta1 = Expression(rule=HPCCondOutleta1_Calc)

        def HPCCondOutleta2_Calc(m):
            return m.HPCCondOutletAL - 3 * m.HPCCondOutletBL ** 2 - 2 * m.HPCCondOutletBL

        test_model.HPCCondOutleta2 = Expression(rule=HPCCondOutleta2_Calc)

        def HPCCondOutleta3_Calc(m):
            return m.HPCCondOutletBL ** 2 + m.HPCCondOutletBL ** 3 - m.HPCCondOutletAL * m.HPCCondOutletBL

        test_model.HPCCondOutleta3 = Expression(rule=HPCCondOutleta3_Calc)

        def HPCCondOutletH0_Calc(m, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.HPCTrayTemp[41] + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (m.PR_DHFORM[comp] + m.PR_CPIGDP1[comp] * (
                        (m.HPCTrayTemp[41] + 273.15) - m.PR_RefTemp[comp]) + 2 * m.PR_CPIGDP2[comp] * \
                               m.PR_CPIGDP3[comp] * (-exp(2 / (m.HPCTrayTemp[41] + 273.15) * m.PR_CPIGDP3[comp]) + \
                                                     exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp])) / (
                                       exp(2 / (m.HPCTrayTemp[41] + 273.15) * m.PR_CPIGDP3[comp]) - 1) / (
                                       exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp]) - 1) - \
                               2 * m.PR_CPIGDP4[comp] * m.PR_CPIGDP5[comp] * (
                                       exp(2 / (m.HPCTrayTemp[41] + 273.15) * m.PR_CPIGDP5[comp]) - exp(
                                   2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp])) / \
                               (exp(2 / (m.HPCTrayTemp[41] + 273.15) * m.PR_CPIGDP5[comp]) + 1) / (
                                       exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp]) + 1))

        test_model.HPCCondOutletH0 = Expression(test_model.Component, rule=HPCCondOutletH0_Calc)

        def HPCCondOutletM_Calc(m, comp1, comp2):
            return m.PR_m[comp1] * m.HPCCondOutletSqrtTr[comp1] / (2 * sqrt(m.HPCCondOutletalpha[comp1]))

        test_model.HPCCondOutletM = Expression(test_model.Component, test_model.Component,
                                               rule=HPCCondOutletM_Calc)  # Calculate the enthalpy of liquid phase mixture %unit: J/mol

        def HPCCondOutletHL_Calc(m):
            return sum([m.HPCVapLvMFrac[41, c] * m.HPCCondOutletH0[c] for c in m.Component]) + \
                m.PR_RGas * (m.HPCTrayTemp[41] + 273.15) * (
                        (m.HPCCondOutletZL - 1) - 1 / (2 ** 1.5) * (1 / m.HPCCondOutletBL) * m.HPCCondOutletZLdivBL * \
                        sum([sum([m.HPCVapLvMFrac[41, c1] * m.HPCVapLvMFrac[41, c2] * m.HPCCondOutletAL * (
                                1 + m.HPCCondOutletM[c2, c1] + \
                                m.HPCCondOutletM[c1, c2]) for c1 in m.Component]) for c2 in m.Component]))

        test_model.HPCCondOutMEtlp = Expression(rule=HPCCondOutletHL_Calc)

        # --------CoolLin Enthalpy---------
        def CoolLinSelfTr_Calc(m, comp):
            return (m.CoolLinTemp + 273.15) / m.PR_CrtcTemp[comp]

        test_model.CoolLinSelfTr = Expression(test_model.Component, rule=CoolLinSelfTr_Calc)

        def CoolLinSelfalpha_Calc(m, comp):
            return (1 + m.PR_m[comp] * (1 - m.CoolLinSelfSqrtTr[comp])) ** 2

        test_model.CoolLinSelfalpha = Expression(test_model.Component, rule=CoolLinSelfalpha_Calc)

        def CoolLinSelfai_Calc(m, comp):
            return 0.45724 * m.CoolLinSelfalpha[comp] * (m.PR_RGas ** 2) * (m.PR_CrtcTemp[comp] ** 2) / \
                m.PR_CrtcPressure[
                    comp]

        test_model.CoolLinSelfai = Expression(test_model.Component, rule=CoolLinSelfai_Calc)

        def CoolLinSelfaij_Calc(m, comp1, comp2):
            return m.CoolLinSelfSqrt_ai[comp1, comp2] * (1 - m.PR_Kij[comp1, comp2])

        test_model.CoolLinSelfaij = Expression(test_model.Component, test_model.Component, rule=CoolLinSelfaij_Calc)

        def CoolLinSelfaL_Calc(m):
            return sum(
                [sum([m.CoolLinSelfaij[c1, c2] * m.HPCVapLvMFrac[41, c1] * m.HPCVapLvMFrac[41, c2] for c1 in
                      m.Component])
                 for c2 in m.Component])

        test_model.CoolLinSelfaL = Expression(rule=CoolLinSelfaL_Calc)

        def CoolLinSelfbL_Calc(m):
            return sum([m.PR_bi[c] * m.HPCVapLvMFrac[41, c] for c in m.Component])

        test_model.CoolLinSelfbL = Expression(rule=CoolLinSelfbL_Calc)

        def CoolLinSelfAL_Calc(m):
            return m.CoolLinSelfaL * m.CoolLinPressure * 1000 / ((m.PR_RGas * (m.CoolLinTemp + 273.15)) ** 2)

        test_model.CoolLinSelfAL = Expression(rule=CoolLinSelfAL_Calc)

        def CoolLinSelfBL_Calc(m):
            return m.CoolLinSelfbL * m.CoolLinPressure * 1000 / (m.PR_RGas * (m.CoolLinTemp + 273.15))

        test_model.CoolLinSelfBL = Expression(rule=CoolLinSelfBL_Calc)

        def CoolLinSelfa1_Calc(m):
            return m.CoolLinSelfBL - 1

        test_model.CoolLinSelfa1 = Expression(rule=CoolLinSelfa1_Calc)

        def CoolLinSelfa2_Calc(m):
            return m.CoolLinSelfAL - 3 * m.CoolLinSelfBL ** 2 - 2 * m.CoolLinSelfBL

        test_model.CoolLinSelfa2 = Expression(rule=CoolLinSelfa2_Calc)

        def CoolLinSelfa3_Calc(m):
            return m.CoolLinSelfBL ** 2 + m.CoolLinSelfBL ** 3 - m.CoolLinSelfAL * m.CoolLinSelfBL

        test_model.CoolLinSelfa3 = Expression(rule=CoolLinSelfa3_Calc)

        def CoolLinSelfH0_Calc(m, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.CoolLinTemp + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (
                        m.PR_DHFORM[comp] + m.PR_CPIGDP1[comp] * ((m.CoolLinTemp + 273.15) - m.PR_RefTemp[comp]) + 2 *
                        m.PR_CPIGDP2[comp] * \
                        m.PR_CPIGDP3[comp] * (-exp(2 / (m.CoolLinTemp + 273.15) * m.PR_CPIGDP3[comp]) + \
                                              exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp])) / (
                                exp(2 / (m.CoolLinTemp + 273.15) * m.PR_CPIGDP3[comp]) - 1) / (
                                exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp]) - 1) - \
                        2 * m.PR_CPIGDP4[comp] * m.PR_CPIGDP5[comp] * (
                                exp(2 / (m.CoolLinTemp + 273.15) * m.PR_CPIGDP5[comp]) - exp(
                            2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp])) / \
                        (exp(2 / (m.CoolLinTemp + 273.15) * m.PR_CPIGDP5[comp]) + 1) / (
                                exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp]) + 1))

        test_model.CoolLinSelfH0 = Expression(test_model.Component, rule=CoolLinSelfH0_Calc)

        def CoolLinSelfM_Calc(m, comp1, comp2):
            return m.PR_m[comp1] * m.CoolLinSelfSqrtTr[comp1] / (2 * sqrt(m.CoolLinSelfalpha[comp1]))

        test_model.CoolLinSelfM = Expression(test_model.Component, test_model.Component,
                                             rule=CoolLinSelfM_Calc)  # Calculate the enthalpy of liquid phase mixture %unit: J/mol

        def CoolLinSelfHL_Calc(m):
            return sum([m.HPCVapLvMFrac[41, c] * m.CoolLinSelfH0[c] for c in m.Component]) + \
                m.PR_RGas * (m.CoolLinTemp + 273.15) * (
                        (m.CoolLinSelfZL - 1) - 1 / (2 ** 1.5) * (1 / m.CoolLinSelfBL) * m.CoolLinSelfZLdivBL * \
                        sum([sum([m.HPCVapLvMFrac[41, c1] * m.HPCVapLvMFrac[41, c2] * m.CoolLinSelfAL * (
                                1 + m.CoolLinSelfM[c2, c1] + \
                                m.CoolLinSelfM[c1, c2]) for c1 in m.Component]) for c2 in m.Component]))

        test_model.CoolLinMEtlp = Expression(rule=CoolLinSelfHL_Calc)

        # --------Throttle Enthalpy---------
        def ThrottleSelfTr_Calc(m, comp):
            return (m.ThrottleTemp + 273.15) / m.PR_CrtcTemp[comp]

        test_model.ThrottleSelfTr = Expression(test_model.Component, rule=ThrottleSelfTr_Calc)

        def ThrottleSelfalpha_Calc(m, comp):
            return (1 + m.PR_m[comp] * (1 - m.ThrottleSelfSqrtTr[comp])) ** 2

        test_model.ThrottleSelfalpha = Expression(test_model.Component, rule=ThrottleSelfalpha_Calc)

        def ThrottleSelfai_Calc(m, comp):
            return 0.45724 * m.ThrottleSelfalpha[comp] * (m.PR_RGas ** 2) * (m.PR_CrtcTemp[comp] ** 2) / \
                m.PR_CrtcPressure[
                    comp]

        test_model.ThrottleSelfai = Expression(test_model.Component, rule=ThrottleSelfai_Calc)

        def ThrottleSelfaij_Calc(m, comp1, comp2):
            return m.ThrottleSelfSqrt_ai[comp1, comp2] * (1 - m.PR_Kij[comp1, comp2])

        test_model.ThrottleSelfaij = Expression(test_model.Component, test_model.Component, rule=ThrottleSelfaij_Calc)

        def ThrottleSelfaL_Calc(m):
            return sum(
                [sum([m.ThrottleSelfaij[c1, c2] * m.ThrottleLiqMFrac[c1] * m.ThrottleLiqMFrac[c2] for c1 in
                      m.Component])
                 for c2 in m.Component])

        test_model.ThrottleSelfaL = Expression(rule=ThrottleSelfaL_Calc)

        def ThrottleSelfbL_Calc(m):
            return sum([m.PR_bi[c] * m.ThrottleLiqMFrac[c] for c in m.Component])

        test_model.ThrottleSelfbL = Expression(rule=ThrottleSelfbL_Calc)

        def ThrottleSelfAL_Calc(m):
            return m.ThrottleSelfaL * m.ThrottlePressure * 1000 / ((m.PR_RGas * (m.ThrottleTemp + 273.15)) ** 2)

        test_model.ThrottleSelfAL = Expression(rule=ThrottleSelfAL_Calc)

        def ThrottleSelfBL_Calc(m):
            return m.ThrottleSelfbL * m.ThrottlePressure * 1000 / (m.PR_RGas * (m.ThrottleTemp + 273.15))

        test_model.ThrottleSelfBL = Expression(rule=ThrottleSelfBL_Calc)

        def ThrottleSelfa1_Calc(m):
            return m.ThrottleSelfBL - 1

        test_model.ThrottleSelfa1 = Expression(rule=ThrottleSelfa1_Calc)

        def ThrottleSelfa2_Calc(m):
            return m.ThrottleSelfAL - 3 * m.ThrottleSelfBL ** 2 - 2 * m.ThrottleSelfBL

        test_model.ThrottleSelfa2 = Expression(rule=ThrottleSelfa2_Calc)

        def ThrottleSelfa3_Calc(m):
            return m.ThrottleSelfBL ** 2 + m.ThrottleSelfBL ** 3 - m.ThrottleSelfAL * m.ThrottleSelfBL

        test_model.ThrottleSelfa3 = Expression(rule=ThrottleSelfa3_Calc)

        def ThrottleSelfaV_Calc(m):
            return sum(
                [sum([m.ThrottleSelfaij[c1, c2] * m.ThrottleVapMFrac[c1] * m.ThrottleVapMFrac[c2] for c1 in
                      m.Component])
                 for c2 in m.Component])

        test_model.ThrottleSelfaV = Expression(rule=ThrottleSelfaV_Calc)

        def ThrottleSelfbV_Calc(m):
            return sum([m.PR_bi[c] * m.ThrottleVapMFrac[c] for c in m.Component])

        test_model.ThrottleSelfbV = Expression(rule=ThrottleSelfbV_Calc)

        def ThrottleSelfAV_Calc(m):
            return m.ThrottleSelfaV * m.ThrottlePressure * 1000 / ((m.PR_RGas * (m.ThrottleTemp + 273.15)) ** 2)

        test_model.ThrottleSelfAV = Expression(rule=ThrottleSelfAV_Calc)

        def ThrottleSelfBV_Calc(m):
            return m.ThrottleSelfbV * m.ThrottlePressure * 1000 / (m.PR_RGas * (m.ThrottleTemp + 273.15))

        test_model.ThrottleSelfBV = Expression(rule=ThrottleSelfBV_Calc)

        def ThrottleSelfb1_Calc(m):
            return m.ThrottleSelfBV - 1

        test_model.ThrottleSelfb1 = Expression(rule=ThrottleSelfb1_Calc)

        def ThrottleSelfb2_Calc(m):
            return m.ThrottleSelfAV - 3 * m.ThrottleSelfBV ** 2 - 2 * m.ThrottleSelfBV

        test_model.ThrottleSelfb2 = Expression(rule=ThrottleSelfb2_Calc)

        def ThrottleSelfb3_Calc(m):
            return m.ThrottleSelfBV ** 2 + m.ThrottleSelfBV ** 3 - m.ThrottleSelfAV * m.ThrottleSelfBV

        test_model.ThrottleSelfb3 = Expression(rule=ThrottleSelfb3_Calc)

        def ThrottleSelfSL_Calc(m, comp):
            return sum([m.ThrottleSelfaij[c, comp] * m.ThrottleLiqMFrac[c] for c in m.Component])

        test_model.ThrottleSelfSL = Expression(test_model.Component, rule=ThrottleSelfSL_Calc)

        def ThrottleSelfSV_Calc(m, comp):
            return sum([m.ThrottleSelfaij[c, comp] * m.ThrottleVapMFrac[c] for c in m.Component])

        test_model.ThrottleSelfSV = Expression(test_model.Component, rule=ThrottleSelfSV_Calc)

        def ThrottleSelfPhiL_Calc(m, comp):
            return exp(m.PR_bi[comp] / m.ThrottleSelfbL * (m.ThrottleSelfZL - 1) - \
                       m.ThrottleSelfZLBL - m.ThrottleSelfAL / m.ThrottleSelfBL / 2 / sqrt(2) * (2 * \
                                                                                                 m.ThrottleSelfSL[
                                                                                                     comp] / m.ThrottleSelfaL -
                                                                                                 m.PR_bi[
                                                                                                     comp] / m.ThrottleSelfbL) * \
                       m.ThrottleSelfZLdivBL)

        test_model.ThrottleSelfPhiL = Expression(test_model.Component, rule=ThrottleSelfPhiL_Calc)

        def ThrottleSelfPhiV_Calc(m, comp):
            return exp(m.PR_bi[comp] / m.ThrottleSelfbV * (m.ThrottleSelfZV - 1) - \
                       m.ThrottleSelfZVBV - m.ThrottleSelfAV / m.ThrottleSelfBV / 2 / sqrt(2) * (2 * \
                                                                                                 m.ThrottleSelfSV[
                                                                                                     comp] / m.ThrottleSelfaV -
                                                                                                 m.PR_bi[
                                                                                                     comp] / m.ThrottleSelfbV) * \
                       m.ThrottleSelfZVdivBV)

        test_model.ThrottleSelfPhiV = Expression(test_model.Component, rule=ThrottleSelfPhiV_Calc)

        def ThrottleSelfH0_Calc(m, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.ThrottleTemp + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (
                        m.PR_DHFORM[comp] + m.PR_CPIGDP1[comp] * ((m.ThrottleTemp + 273.15) - m.PR_RefTemp[comp]) + 2 *
                        m.PR_CPIGDP2[comp] * \
                        m.PR_CPIGDP3[comp] * (-exp(2 / (m.ThrottleTemp + 273.15) * m.PR_CPIGDP3[comp]) + \
                                              exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp])) / (
                                exp(2 / (m.ThrottleTemp + 273.15) * m.PR_CPIGDP3[comp]) - 1) / (
                                exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp]) - 1) - \
                        2 * m.PR_CPIGDP4[comp] * m.PR_CPIGDP5[comp] * (
                                exp(2 / (m.ThrottleTemp + 273.15) * m.PR_CPIGDP5[comp]) - exp(
                            2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp])) / \
                        (exp(2 / (m.ThrottleTemp + 273.15) * m.PR_CPIGDP5[comp]) + 1) / (
                                exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp]) + 1))

        test_model.ThrottleSelfH0 = Expression(test_model.Component, rule=ThrottleSelfH0_Calc)

        def ThrottleSelfM_Calc(m, comp1, comp2):
            return m.PR_m[comp1] * m.ThrottleSelfSqrtTr[comp1] / (2 * sqrt(m.ThrottleSelfalpha[comp1]))

        test_model.ThrottleSelfM = Expression(test_model.Component, test_model.Component,
                                              rule=ThrottleSelfM_Calc)  # Calculate the enthalpy of liquid phase mixture %unit: J/mol

        def ThrottleSelfHV_Calc(m):
            return sum([m.ThrottleVapMFrac[c] * m.ThrottleSelfH0[c] for c in m.Component]) + \
                m.PR_RGas * (m.ThrottleTemp + 273.15) * (
                        (m.ThrottleSelfZV - 1) - 1 / (2 ** 1.5) * (1 / m.ThrottleSelfBV) * m.ThrottleSelfZVdivBV * \
                        sum([sum([m.ThrottleVapMFrac[c1] * m.ThrottleVapMFrac[c2] * m.ThrottleSelfAV * (
                                1 + m.ThrottleSelfM[c2, c1] + \
                                m.ThrottleSelfM[c1, c2]) for c1 in m.Component]) for c2 in m.Component]))

        test_model.ThrottleVapMEtlp = Expression(rule=ThrottleSelfHV_Calc)

        # Calculate the enthalpy of liquid phase mixture %unit: J/mol
        def ThrottleSelfHL_Calc(m):
            return sum([m.ThrottleLiqMFrac[c] * m.ThrottleSelfH0[c] for c in m.Component]) + \
                m.PR_RGas * (m.ThrottleTemp + 273.15) * (
                        (m.ThrottleSelfZL - 1) - 1 / (2 ** 1.5) * (1 / m.ThrottleSelfBL) * m.ThrottleSelfZLdivBL * \
                        sum([sum([m.ThrottleLiqMFrac[c1] * m.ThrottleLiqMFrac[c2] * m.ThrottleSelfAL * (
                                1 + m.ThrottleSelfM[c2, c1] + \
                                m.ThrottleSelfM[c1, c2]) for c1 in m.Component]) for c2 in m.Component]))

        test_model.ThrottleLiqMEtlp = Expression(rule=ThrottleSelfHL_Calc)

        def ThrottleMEtlp(m):
            return (
                        m.ThrottleLiqMEtlp * m.ThrottleLiqMFlow + m.ThrottleVapMEtlp * m.ThrottleVapMFlow) / m.HPCSumpOutMFlow

        test_model.ThrottleMEtlp = Expression(rule=ThrottleMEtlp)

        # --------HeatLPC Enthalpy---------
        def HeatLPCSelfTr_Calc(m, comp):
            return (m.HeatLPCTemp + 273.15) / m.PR_CrtcTemp[comp]

        test_model.HeatLPCSelfTr = Expression(test_model.Component, rule=HeatLPCSelfTr_Calc)

        def HeatLPCSelfalpha_Calc(m, comp):
            return (1 + m.PR_m[comp] * (1 - m.HeatLPCSelfSqrtTr[comp])) ** 2

        test_model.HeatLPCSelfalpha = Expression(test_model.Component, rule=HeatLPCSelfalpha_Calc)

        def HeatLPCSelfai_Calc(m, comp):
            return 0.45724 * m.HeatLPCSelfalpha[comp] * (m.PR_RGas ** 2) * (m.PR_CrtcTemp[comp] ** 2) / \
                m.PR_CrtcPressure[
                    comp]

        test_model.HeatLPCSelfai = Expression(test_model.Component, rule=HeatLPCSelfai_Calc)

        def HeatLPCSelfaij_Calc(m, comp1, comp2):
            return m.HeatLPCSelfSqrt_ai[comp1, comp2] * (1 - m.PR_Kij[comp1, comp2])

        test_model.HeatLPCSelfaij = Expression(test_model.Component, test_model.Component, rule=HeatLPCSelfaij_Calc)

        def HeatLPCSelfaL_Calc(m):
            return sum(
                [sum([m.HeatLPCSelfaij[c1, c2] * m.HeatLPCLiqMFrac[c1] * m.HeatLPCLiqMFrac[c2] for c1 in m.Component])
                 for
                 c2 in m.Component])

        test_model.HeatLPCSelfaL = Expression(rule=HeatLPCSelfaL_Calc)

        def HeatLPCSelfbL_Calc(m):
            return sum([m.PR_bi[c] * m.HeatLPCLiqMFrac[c] for c in m.Component])

        test_model.HeatLPCSelfbL = Expression(rule=HeatLPCSelfbL_Calc)

        def HeatLPCSelfAL_Calc(m):
            return m.HeatLPCSelfaL * m.HeatLPCPressure * 1000 / ((m.PR_RGas * (m.HeatLPCTemp + 273.15)) ** 2)

        test_model.HeatLPCSelfAL = Expression(rule=HeatLPCSelfAL_Calc)

        def HeatLPCSelfBL_Calc(m):
            return m.HeatLPCSelfbL * m.HeatLPCPressure * 1000 / (m.PR_RGas * (m.HeatLPCTemp + 273.15))

        test_model.HeatLPCSelfBL = Expression(rule=HeatLPCSelfBL_Calc)

        def HeatLPCSelfa1_Calc(m):
            return m.HeatLPCSelfBL - 1

        test_model.HeatLPCSelfa1 = Expression(rule=HeatLPCSelfa1_Calc)

        def HeatLPCSelfa2_Calc(m):
            return m.HeatLPCSelfAL - 3 * m.HeatLPCSelfBL ** 2 - 2 * m.HeatLPCSelfBL

        test_model.HeatLPCSelfa2 = Expression(rule=HeatLPCSelfa2_Calc)

        def HeatLPCSelfa3_Calc(m):
            return m.HeatLPCSelfBL ** 2 + m.HeatLPCSelfBL ** 3 - m.HeatLPCSelfAL * m.HeatLPCSelfBL

        test_model.HeatLPCSelfa3 = Expression(rule=HeatLPCSelfa3_Calc)

        def HeatLPCSelfaV_Calc(m):
            return sum(
                [sum([m.HeatLPCSelfaij[c1, c2] * m.HeatLPCVapMFrac[c1] * m.HeatLPCVapMFrac[c2] for c1 in m.Component])
                 for
                 c2 in m.Component])

        test_model.HeatLPCSelfaV = Expression(rule=HeatLPCSelfaV_Calc)

        def HeatLPCSelfbV_Calc(m):
            return sum([m.PR_bi[c] * m.HeatLPCVapMFrac[c] for c in m.Component])

        test_model.HeatLPCSelfbV = Expression(rule=HeatLPCSelfbV_Calc)

        def HeatLPCSelfAV_Calc(m):
            return m.HeatLPCSelfaV * m.HeatLPCPressure * 1000 / ((m.PR_RGas * (m.HeatLPCTemp + 273.15)) ** 2)

        test_model.HeatLPCSelfAV = Expression(rule=HeatLPCSelfAV_Calc)

        def HeatLPCSelfBV_Calc(m):
            return m.HeatLPCSelfbV * m.HeatLPCPressure * 1000 / (m.PR_RGas * (m.HeatLPCTemp + 273.15))

        test_model.HeatLPCSelfBV = Expression(rule=HeatLPCSelfBV_Calc)

        def HeatLPCSelfb1_Calc(m):
            return m.HeatLPCSelfBV - 1

        test_model.HeatLPCSelfb1 = Expression(rule=HeatLPCSelfb1_Calc)

        def HeatLPCSelfb2_Calc(m):
            return m.HeatLPCSelfAV - 3 * m.HeatLPCSelfBV ** 2 - 2 * m.HeatLPCSelfBV

        test_model.HeatLPCSelfb2 = Expression(rule=HeatLPCSelfb2_Calc)

        def HeatLPCSelfb3_Calc(m):
            return m.HeatLPCSelfBV ** 2 + m.HeatLPCSelfBV ** 3 - m.HeatLPCSelfAV * m.HeatLPCSelfBV

        test_model.HeatLPCSelfb3 = Expression(rule=HeatLPCSelfb3_Calc)

        def HeatLPCSelfSL_Calc(m, comp):
            return sum([m.HeatLPCSelfaij[c, comp] * m.HeatLPCLiqMFrac[c] for c in m.Component])

        test_model.HeatLPCSelfSL = Expression(test_model.Component, rule=HeatLPCSelfSL_Calc)

        def HeatLPCSelfSV_Calc(m, comp):
            return sum([m.HeatLPCSelfaij[c, comp] * m.HeatLPCVapMFrac[c] for c in m.Component])

        test_model.HeatLPCSelfSV = Expression(test_model.Component, rule=HeatLPCSelfSV_Calc)

        def HeatLPCSelfPhiL_Calc(m, comp):
            return exp(m.PR_bi[comp] / m.HeatLPCSelfbL * (m.HeatLPCSelfZL - 1) - \
                       m.HeatLPCSelfZLBL - m.HeatLPCSelfAL / m.HeatLPCSelfBL / 2 / sqrt(2) * (2 * \
                                                                                              m.HeatLPCSelfSL[
                                                                                                  comp] / m.HeatLPCSelfaL -
                                                                                              m.PR_bi[
                                                                                                  comp] / m.HeatLPCSelfbL) * \
                       m.HeatLPCSelfZLdivBL)

        test_model.HeatLPCSelfPhiL = Expression(test_model.Component, rule=HeatLPCSelfPhiL_Calc)

        def HeatLPCSelfPhiV_Calc(m, comp):
            return exp(m.PR_bi[comp] / m.HeatLPCSelfbV * (m.HeatLPCSelfZV - 1) - \
                       m.HeatLPCSelfZVBV - m.HeatLPCSelfAV / m.HeatLPCSelfBV / 2 / sqrt(2) * (2 * \
                                                                                              m.HeatLPCSelfSV[
                                                                                                  comp] / m.HeatLPCSelfaV -
                                                                                              m.PR_bi[
                                                                                                  comp] / m.HeatLPCSelfbV) * \
                       m.HeatLPCSelfZVdivBV)

        test_model.HeatLPCSelfPhiV = Expression(test_model.Component, rule=HeatLPCSelfPhiV_Calc)

        def HeatLPCSelfH0_Calc(m, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.HeatLPCTemp + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (
                        m.PR_DHFORM[comp] + m.PR_CPIGDP1[comp] * ((m.HeatLPCTemp + 273.15) - m.PR_RefTemp[comp]) + 2 *
                        m.PR_CPIGDP2[comp] * \
                        m.PR_CPIGDP3[comp] * (-exp(2 / (m.HeatLPCTemp + 273.15) * m.PR_CPIGDP3[comp]) + \
                                              exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp])) / (
                                exp(2 / (m.HeatLPCTemp + 273.15) * m.PR_CPIGDP3[comp]) - 1) / (
                                exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp]) - 1) - \
                        2 * m.PR_CPIGDP4[comp] * m.PR_CPIGDP5[comp] * (
                                exp(2 / (m.HeatLPCTemp + 273.15) * m.PR_CPIGDP5[comp]) - exp(
                            2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp])) / \
                        (exp(2 / (m.HeatLPCTemp + 273.15) * m.PR_CPIGDP5[comp]) + 1) / (
                                exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp]) + 1))

        test_model.HeatLPCSelfH0 = Expression(test_model.Component, rule=HeatLPCSelfH0_Calc)

        def HeatLPCSelfM_Calc(m, comp1, comp2):
            return m.PR_m[comp1] * m.HeatLPCSelfSqrtTr[comp1] / (2 * sqrt(m.HeatLPCSelfalpha[comp1]))

        test_model.HeatLPCSelfM = Expression(test_model.Component, test_model.Component,
                                             rule=HeatLPCSelfM_Calc)  # Calculate the enthalpy of liquid phase mixture %unit: J/mol

        def HeatLPCSelfHV_Calc(m):
            return sum([m.HeatLPCVapMFrac[c] * m.HeatLPCSelfH0[c] for c in m.Component]) + \
                m.PR_RGas * (m.HeatLPCTemp + 273.15) * (
                        (m.HeatLPCSelfZV - 1) - 1 / (2 ** 1.5) * (1 / m.HeatLPCSelfBV) * m.HeatLPCSelfZVdivBV * \
                        sum([sum([m.HeatLPCVapMFrac[c1] * m.HeatLPCVapMFrac[c2] * m.HeatLPCSelfAV * (
                                1 + m.HeatLPCSelfM[c2, c1] + \
                                m.HeatLPCSelfM[c1, c2]) for c1 in m.Component]) for c2 in m.Component]))

        test_model.HeatLPCVapMEtlp = Expression(rule=HeatLPCSelfHV_Calc)

        # Calculate the enthalpy of liquid phase mixture %unit: J/mol
        def HeatLPCSelfHL_Calc(m):
            return sum([m.HeatLPCLiqMFrac[c] * m.HeatLPCSelfH0[c] for c in m.Component]) + \
                m.PR_RGas * (m.HeatLPCTemp + 273.15) * (
                        (m.HeatLPCSelfZL - 1) - 1 / (2 ** 1.5) * (1 / m.HeatLPCSelfBL) * m.HeatLPCSelfZLdivBL * \
                        sum([sum([m.HeatLPCLiqMFrac[c1] * m.HeatLPCLiqMFrac[c2] * m.HeatLPCSelfAL * (
                                1 + m.HeatLPCSelfM[c2, c1] + \
                                m.HeatLPCSelfM[c1, c2]) for c1 in m.Component]) for c2 in m.Component]))

        test_model.HeatLPCLiqMEtlp = Expression(rule=HeatLPCSelfHL_Calc)

        def HeatLPCMEtlp(m):
            return (m.HeatLPCLiqMEtlp * m.HeatLPCLiqMFlow + m.HeatLPCVapMEtlp * m.HeatLPCVapMFlow) / m.HPCSumpOutMFlow

        test_model.HeatLPCMEtlp = Expression(rule=HeatLPCMEtlp)

        # --------LPC Component Flowrate---------
        def LPCVapLvMCompFlow(m, tray, comp):
            return m.LPCVapLvMFrac[tray, comp] * m.LPCVapLvMFlow[tray]

        test_model.LPCVapLvMCompFlow = Expression(test_model.LPCTrays, test_model.Component, rule=LPCVapLvMCompFlow)

        def LPCLiqLvMCompFlow(m, tray, comp):
            return m.LPCLiqLvMFrac[tray, comp] * m.LPCLiqLvMFlow[tray]

        test_model.LPCLiqLvMCompFlow = Expression(test_model.LPCTrays, test_model.Component, rule=LPCLiqLvMCompFlow)

        # --------LPC Thermo and Enthalpy---------
        def LPCAllTraysTr_Calc(m, tray, comp):
            return (m.LPCTrayTemp[tray] + 273.15) / m.PR_CrtcTemp[comp]

        test_model.LPCAllTraysTr = Expression(test_model.LPCTrays, test_model.Component, rule=LPCAllTraysTr_Calc)

        def LPCAllTraysalpha_Calc(m, tray, comp):
            return (1 + m.PR_m[comp] * (1 - m.LPCAllTraysSqrtTr[tray, comp])) ** 2

        test_model.LPCAllTraysalpha = Expression(test_model.LPCTrays, test_model.Component, rule=LPCAllTraysalpha_Calc)

        def LPCAllTraysai_Calc(m, tray, comp):
            return 0.45724 * m.LPCAllTraysalpha[tray, comp] * (m.PR_RGas ** 2) * (m.PR_CrtcTemp[comp] ** 2) / \
                m.PR_CrtcPressure[comp]

        test_model.LPCAllTraysai = Expression(test_model.LPCTrays, test_model.Component, rule=LPCAllTraysai_Calc)

        def LPCAllTraysaij_Calc(m, tray, comp1, comp2):
            return m.LPCAllTraysSqrt_ai[tray, comp1, comp2] * (1 - m.PR_Kij[comp1, comp2])

        test_model.LPCAllTraysaij = Expression(test_model.LPCTrays, test_model.Component, test_model.Component,
                                               rule=LPCAllTraysaij_Calc)

        def LPCAllTraysaL_Calc(m, tray):
            return sum(
                [sum([m.LPCAllTraysaij[tray, c1, c2] * m.LPCLiqLvMFrac[tray, c1] * m.LPCLiqLvMFrac[tray, c2] for c1
                      in m.Component]) for c2 in m.Component])

        test_model.LPCAllTraysaL = Expression(test_model.LPCTrays, rule=LPCAllTraysaL_Calc)

        def LPCAllTraysbL_Calc(m, tray):
            return sum([m.PR_bi[c] * m.LPCLiqLvMFrac[tray, c] for c in m.Component])

        test_model.LPCAllTraysbL = Expression(test_model.LPCTrays, rule=LPCAllTraysbL_Calc)

        def LPCAllTraysAL_Calc(m, tray):
            return m.LPCAllTraysaL[tray] * m.LPCTrayPressure[tray] * 1000 / (
                    (m.PR_RGas * (m.LPCTrayTemp[tray] + 273.15)) ** 2)

        test_model.LPCAllTraysAL = Expression(test_model.LPCTrays, rule=LPCAllTraysAL_Calc)

        def LPCAllTraysBL_Calc(m, tray):
            return m.LPCAllTraysbL[tray] * m.LPCTrayPressure[tray] * 1000 / (m.PR_RGas * (m.LPCTrayTemp[tray] + 273.15))

        test_model.LPCAllTraysBL = Expression(test_model.LPCTrays, rule=LPCAllTraysBL_Calc)

        def LPCAllTraysa1_Calc(m, tray):
            return m.LPCAllTraysBL[tray] - 1

        test_model.LPCAllTraysa1 = Expression(test_model.LPCTrays, rule=LPCAllTraysa1_Calc)

        def LPCAllTraysa2_Calc(m, tray):
            return m.LPCAllTraysAL[tray] - 3 * m.LPCAllTraysBL[tray] ** 2 - 2 * m.LPCAllTraysBL[tray]

        test_model.LPCAllTraysa2 = Expression(test_model.LPCTrays, rule=LPCAllTraysa2_Calc)

        def LPCAllTraysa3_Calc(m, tray):
            return m.LPCAllTraysBL[tray] ** 2 + m.LPCAllTraysBL[tray] ** 3 - m.LPCAllTraysAL[tray] * m.LPCAllTraysBL[
                tray]

        test_model.LPCAllTraysa3 = Expression(test_model.LPCTrays, rule=LPCAllTraysa3_Calc)

        def LPCAllTraysaV_Calc(m, tray):
            return sum(
                [sum([m.LPCAllTraysaij[tray, c1, c2] * m.LPCVapLvMFrac[tray, c1] * m.LPCVapLvMFrac[tray, c2] for c1
                      in m.Component]) for c2 in m.Component])

        test_model.LPCAllTraysaV = Expression(test_model.LPCTrays, rule=LPCAllTraysaV_Calc)

        def LPCAllTraysbV_Calc(m, tray):
            return sum([m.PR_bi[c] * m.LPCVapLvMFrac[tray, c] for c in m.Component])

        test_model.LPCAllTraysbV = Expression(test_model.LPCTrays, rule=LPCAllTraysbV_Calc)

        def LPCAllTraysAV_Calc(m, tray):
            return m.LPCAllTraysaV[tray] * m.LPCTrayPressure[tray] * 1000 / (
                    (m.PR_RGas * (m.LPCTrayTemp[tray] + 273.15)) ** 2)

        test_model.LPCAllTraysAV = Expression(test_model.LPCTrays, rule=LPCAllTraysAV_Calc)

        def LPCAllTraysBV_Calc(m, tray):
            return m.LPCAllTraysbV[tray] * m.LPCTrayPressure[tray] * 1000 / (m.PR_RGas * (m.LPCTrayTemp[tray] + 273.15))

        test_model.LPCAllTraysBV = Expression(test_model.LPCTrays, rule=LPCAllTraysBV_Calc)

        def LPCAllTraysb1_Calc(m, tray):
            return m.LPCAllTraysBV[tray] - 1

        test_model.LPCAllTraysb1 = Expression(test_model.LPCTrays, rule=LPCAllTraysb1_Calc)

        def LPCAllTraysb2_Calc(m, tray):
            return m.LPCAllTraysAV[tray] - 3 * m.LPCAllTraysBV[tray] ** 2 - 2 * m.LPCAllTraysBV[tray]

        test_model.LPCAllTraysb2 = Expression(test_model.LPCTrays, rule=LPCAllTraysb2_Calc)

        def LPCAllTraysb3_Calc(m, tray):
            return m.LPCAllTraysBV[tray] ** 2 + m.LPCAllTraysBV[tray] ** 3 - m.LPCAllTraysAV[tray] * m.LPCAllTraysBV[
                tray]

        test_model.LPCAllTraysb3 = Expression(test_model.LPCTrays, rule=LPCAllTraysb3_Calc)

        def LPCAllTraysSL_Calc(m, tray, comp):
            return sum([m.LPCAllTraysaij[tray, c, comp] * m.LPCLiqLvMFrac[tray, c] for c in m.Component])

        test_model.LPCAllTraysSL = Expression(test_model.LPCTrays, test_model.Component, rule=LPCAllTraysSL_Calc)

        def LPCAllTraysSV_Calc(m, tray, comp):
            return sum([m.LPCAllTraysaij[tray, c, comp] * m.LPCVapLvMFrac[tray, c] for c in m.Component])

        test_model.LPCAllTraysSV = Expression(test_model.LPCTrays, test_model.Component, rule=LPCAllTraysSV_Calc)

        def LPCAllTraysPhiL_Calc(m, tray, comp):
            return exp(m.PR_bi[comp] / m.LPCAllTraysbL[tray] * (m.LPCAllTraysZL[tray] - 1) - \
                       m.LPCAllTraysZLBL[tray] - m.LPCAllTraysAL[tray] / m.LPCAllTraysBL[tray] / 2 / sqrt(2) * (2 * \
                                                                                                                m.LPCAllTraysSL[
                                                                                                                    tray, comp] /
                                                                                                                m.LPCAllTraysaL[
                                                                                                                    tray] -
                                                                                                                m.PR_bi[
                                                                                                                    comp] /
                                                                                                                m.LPCAllTraysbL[
                                                                                                                    tray]) * \
                       m.LPCAllTraysZLdivBL[tray])

        test_model.LPCAllTraysPhiL = Expression(test_model.LPCTrays, test_model.Component, rule=LPCAllTraysPhiL_Calc)

        def LPCAllTraysPhiV_Calc(m, tray, comp):
            return exp(m.PR_bi[comp] / m.LPCAllTraysbV[tray] * (m.LPCAllTraysZV[tray] - 1) - \
                       m.LPCAllTraysZVBV[tray] - m.LPCAllTraysAV[tray] / m.LPCAllTraysBV[tray] / 2 / sqrt(2) * (2 * \
                                                                                                                m.LPCAllTraysSV[
                                                                                                                    tray, comp] /
                                                                                                                m.LPCAllTraysaV[
                                                                                                                    tray] -
                                                                                                                m.PR_bi[
                                                                                                                    comp] /
                                                                                                                m.LPCAllTraysbV[
                                                                                                                    tray]) * \
                       m.LPCAllTraysZVdivBV[tray])

        test_model.LPCAllTraysPhiV = Expression(test_model.LPCTrays, test_model.Component, rule=LPCAllTraysPhiV_Calc)

        def LPCAllTraysH0_Calc(m, tray, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.LPCTrayTemp[tray] + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (m.PR_DHFORM[comp] + m.PR_CPIGDP1[comp] * (
                        (m.LPCTrayTemp[tray] + 273.15) - m.PR_RefTemp[comp]) + 2 * m.PR_CPIGDP2[comp] * \
                               m.PR_CPIGDP3[comp] * (-exp(2 / (m.LPCTrayTemp[tray] + 273.15) * m.PR_CPIGDP3[comp]) + \
                                                     exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp])) / (
                                       exp(2 / (m.LPCTrayTemp[tray] + 273.15) * m.PR_CPIGDP3[comp]) - 1) / (
                                       exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp]) - 1) - \
                               2 * m.PR_CPIGDP4[comp] * m.PR_CPIGDP5[comp] * (
                                       exp(2 / (m.LPCTrayTemp[tray] + 273.15) * m.PR_CPIGDP5[comp]) - exp(
                                   2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp])) / \
                               (exp(2 / (m.LPCTrayTemp[tray] + 273.15) * m.PR_CPIGDP5[comp]) + 1) / (
                                       exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp]) + 1))

        test_model.LPCAllTraysH0 = Expression(test_model.LPCTrays, test_model.Component, rule=LPCAllTraysH0_Calc)

        def LPCAllTraysM_Calc(m, tray, comp1, comp2):
            return m.PR_m[comp1] * m.LPCAllTraysSqrtTr[tray, comp1] / (2 * sqrt(m.LPCAllTraysalpha[tray, comp1]))

        test_model.LPCAllTraysM = Expression(test_model.LPCTrays, test_model.Component, test_model.Component,
                                             rule=LPCAllTraysM_Calc)  # Calculate the enthalpy of liquid phase mixture %unit: J/mol

        def LPCAllTraysHV_Calc(m, tray):
            return sum([m.LPCVapLvMFrac[tray, c] * m.LPCAllTraysH0[tray, c] for c in m.Component]) + \
                m.PR_RGas * (m.LPCTrayTemp[tray] + 273.15) * (
                        (m.LPCAllTraysZV[tray] - 1) - 1 / (2 ** 1.5) * (1 / m.LPCAllTraysBV[tray]) *
                        m.LPCAllTraysZVdivBV[tray] * \
                        sum([sum([m.LPCVapLvMFrac[tray, c1] * m.LPCVapLvMFrac[tray, c2] * m.LPCAllTraysAV[tray] * (
                                1 + m.LPCAllTraysM[tray, c2, c1] + \
                                m.LPCAllTraysM[tray, c1, c2]) for c1 in m.Component]) for c2 in m.Component]))

        test_model.LPCVapLvMEtlp = Expression(test_model.LPCTrays, rule=LPCAllTraysHV_Calc)

        # Calculate the enthalpy of liquid phase mixture %unit: J/mol
        def LPCAllTraysHL_Calc(m, tray):
            return sum([m.LPCLiqLvMFrac[tray, c] * m.LPCAllTraysH0[tray, c] for c in m.Component]) + \
                m.PR_RGas * (m.LPCTrayTemp[tray] + 273.15) * (
                        (m.LPCAllTraysZL[tray] - 1) - 1 / (2 ** 1.5) * (1 / m.LPCAllTraysBL[tray]) *
                        m.LPCAllTraysZLdivBL[tray] * \
                        sum([sum([m.LPCLiqLvMFrac[tray, c1] * m.LPCLiqLvMFrac[tray, c2] * m.LPCAllTraysAL[tray] * (
                                1 + m.LPCAllTraysM[tray, c2, c1] + \
                                m.LPCAllTraysM[tray, c1, c2]) for c1 in m.Component]) for c2 in m.Component]))

        test_model.LPCLiqLvMEtlp = Expression(test_model.LPCTrays, rule=LPCAllTraysHL_Calc)

        # --------LPC Extractions---------
        def LPCExt46MFracSpec(m, comp):
            return m.LPCVapLvMFrac[46, comp]

        test_model.LPCExt46MFrac = Expression(test_model.Component, rule=LPCExt46MFracSpec)

        def LPCExt46MEtlpSpec(m):
            return m.LPCVapLvMEtlp[46]

        test_model.LPCExt46MEtlp = Expression(rule=LPCExt46MEtlpSpec)

        def LPCExt15MFracSpec(m, comp):
            return m.LPCVapLvMFrac[15, comp]

        test_model.LPCExt15MFrac = Expression(test_model.Component, rule=LPCExt15MFracSpec)

        def LPCExt15MEtlpSpec(m):
            return m.LPCVapLvMEtlp[15]

        test_model.LPCExt15MEtlp = Expression(rule=LPCExt15MEtlpSpec)

        # --------LPCReboiler VLE & Enthalpy---------
        def LPCReboilerHdlpTr_Calc(m, comp):
            return (m.LPCReboilerTemp + 273.15) / m.PR_CrtcTemp[comp]

        test_model.LPCReboilerHdlpTr = Expression(test_model.Component, rule=LPCReboilerHdlpTr_Calc)

        def LPCReboilerHdlpalpha_Calc(m, comp):
            return (1 + m.PR_m[comp] * (1 - m.LPCReboilerHdlpSqrtTr[comp])) ** 2

        test_model.LPCReboilerHdlpalpha = Expression(test_model.Component, rule=LPCReboilerHdlpalpha_Calc)

        def LPCReboilerHdlpai_Calc(m, comp):
            return 0.45724 * m.LPCReboilerHdlpalpha[comp] * (m.PR_RGas ** 2) * (m.PR_CrtcTemp[comp] ** 2) / \
                m.PR_CrtcPressure[comp]

        test_model.LPCReboilerHdlpai = Expression(test_model.Component, rule=LPCReboilerHdlpai_Calc)

        def LPCReboilerHdlpaij_Calc(m, comp1, comp2):
            return m.LPCReboilerHdlpSqrt_ai[comp1, comp2] * (1 - m.PR_Kij[comp1, comp2])

        test_model.LPCReboilerHdlpaij = Expression(test_model.Component, test_model.Component,
                                                   rule=LPCReboilerHdlpaij_Calc)

        def LPCReboilerHdlpaL_Calc(m):
            return sum(
                [sum([m.LPCReboilerHdlpaij[c1, c2] * m.LPCReboilerLiqLvMFrac[c1] * m.LPCReboilerLiqLvMFrac[c2] for c1
                      in m.Component]) for c2 in m.Component])

        test_model.LPCReboilerHdlpaL = Expression(rule=LPCReboilerHdlpaL_Calc)

        def LPCReboilerHdlpbL_Calc(m):
            return sum([m.PR_bi[c] * m.LPCReboilerLiqLvMFrac[c] for c in m.Component])

        test_model.LPCReboilerHdlpbL = Expression(rule=LPCReboilerHdlpbL_Calc)

        def LPCReboilerHdlpAL_Calc(m):
            return m.LPCReboilerHdlpaL * m.LPCReboilerPressure * 1000 / (
                        (m.PR_RGas * (m.LPCReboilerTemp + 273.15)) ** 2)

        test_model.LPCReboilerHdlpAL = Expression(rule=LPCReboilerHdlpAL_Calc)

        def LPCReboilerHdlpBL_Calc(m):
            return m.LPCReboilerHdlpbL * m.LPCReboilerPressure * 1000 / (m.PR_RGas * (m.LPCReboilerTemp + 273.15))

        test_model.LPCReboilerHdlpBL = Expression(rule=LPCReboilerHdlpBL_Calc)

        def LPCReboilerHdlpa1_Calc(m):
            return m.LPCReboilerHdlpBL - 1

        test_model.LPCReboilerHdlpa1 = Expression(rule=LPCReboilerHdlpa1_Calc)

        def LPCReboilerHdlpa2_Calc(m):
            return m.LPCReboilerHdlpAL - 3 * m.LPCReboilerHdlpBL ** 2 - 2 * m.LPCReboilerHdlpBL

        test_model.LPCReboilerHdlpa2 = Expression(rule=LPCReboilerHdlpa2_Calc)

        def LPCReboilerHdlpa3_Calc(m):
            return m.LPCReboilerHdlpBL ** 2 + m.LPCReboilerHdlpBL ** 3 - m.LPCReboilerHdlpAL * m.LPCReboilerHdlpBL

        test_model.LPCReboilerHdlpa3 = Expression(rule=LPCReboilerHdlpa3_Calc)

        def LPCReboilerHdlpaV_Calc(m):
            return sum(
                [sum([m.LPCReboilerHdlpaij[c1, c2] * m.LPCReboilerVapLvMFrac[c1] * m.LPCReboilerVapLvMFrac[c2] for c1
                      in m.Component]) for c2 in m.Component])

        test_model.LPCReboilerHdlpaV = Expression(rule=LPCReboilerHdlpaV_Calc)

        def LPCReboilerHdlpbV_Calc(m):
            return sum([m.PR_bi[c] * m.LPCReboilerVapLvMFrac[c] for c in m.Component])

        test_model.LPCReboilerHdlpbV = Expression(rule=LPCReboilerHdlpbV_Calc)

        def LPCReboilerHdlpAV_Calc(m):
            return m.LPCReboilerHdlpaV * m.LPCReboilerPressure * 1000 / (
                        (m.PR_RGas * (m.LPCReboilerTemp + 273.15)) ** 2)

        test_model.LPCReboilerHdlpAV = Expression(rule=LPCReboilerHdlpAV_Calc)

        def LPCReboilerHdlpBV_Calc(m):
            return m.LPCReboilerHdlpbV * m.LPCReboilerPressure * 1000 / (m.PR_RGas * (m.LPCReboilerTemp + 273.15))

        test_model.LPCReboilerHdlpBV = Expression(rule=LPCReboilerHdlpBV_Calc)

        def LPCReboilerHdlpb1_Calc(m):
            return m.LPCReboilerHdlpBV - 1

        test_model.LPCReboilerHdlpb1 = Expression(rule=LPCReboilerHdlpb1_Calc)

        def LPCReboilerHdlpb2_Calc(m):
            return m.LPCReboilerHdlpAV - 3 * m.LPCReboilerHdlpBV ** 2 - 2 * m.LPCReboilerHdlpBV

        test_model.LPCReboilerHdlpb2 = Expression(rule=LPCReboilerHdlpb2_Calc)

        def LPCReboilerHdlpb3_Calc(m):
            return m.LPCReboilerHdlpBV ** 2 + m.LPCReboilerHdlpBV ** 3 - m.LPCReboilerHdlpAV * m.LPCReboilerHdlpBV

        test_model.LPCReboilerHdlpb3 = Expression(rule=LPCReboilerHdlpb3_Calc)

        def LPCReboilerHdlpSL_Calc(m, comp):
            return sum([m.LPCReboilerHdlpaij[c, comp] * m.LPCReboilerLiqLvMFrac[c] for c in m.Component])

        test_model.LPCReboilerHdlpSL = Expression(test_model.Component, rule=LPCReboilerHdlpSL_Calc)

        def LPCReboilerHdlpSV_Calc(m, comp):
            return sum([m.LPCReboilerHdlpaij[c, comp] * m.LPCReboilerVapLvMFrac[c] for c in m.Component])

        test_model.LPCReboilerHdlpSV = Expression(test_model.Component, rule=LPCReboilerHdlpSV_Calc)

        def LPCReboilerHdlpPhiL_Calc(m, comp):
            return exp(m.PR_bi[comp] / m.LPCReboilerHdlpbL * (m.LPCReboilerHdlpZL - 1) - \
                       m.LPCReboilerHdlpZLBL - m.LPCReboilerHdlpAL / m.LPCReboilerHdlpBL / 2 / sqrt(2) * (2 * \
                                                                                                          m.LPCReboilerHdlpSL[
                                                                                                              comp] / m.LPCReboilerHdlpaL -
                                                                                                          m.PR_bi[
                                                                                                              comp] / m.LPCReboilerHdlpbL) * \
                       m.LPCReboilerHdlpZLdivBL)

        test_model.LPCReboilerHdlpPhiL = Expression(test_model.Component, rule=LPCReboilerHdlpPhiL_Calc)

        def LPCReboilerHdlpPhiV_Calc(m, comp):
            return exp(m.PR_bi[comp] / m.LPCReboilerHdlpbV * (m.LPCReboilerHdlpZV - 1) - \
                       m.LPCReboilerHdlpZVBV - m.LPCReboilerHdlpAV / m.LPCReboilerHdlpBV / 2 / sqrt(2) * (2 * \
                                                                                                          m.LPCReboilerHdlpSV[
                                                                                                              comp] / m.LPCReboilerHdlpaV -
                                                                                                          m.PR_bi[
                                                                                                              comp] / m.LPCReboilerHdlpbV) * \
                       m.LPCReboilerHdlpZVdivBV)

        test_model.LPCReboilerHdlpPhiV = Expression(test_model.Component, rule=LPCReboilerHdlpPhiV_Calc)

        def LPCReboilerHdlpH0_Calc(m, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.LPCReboilerTemp + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (m.PR_DHFORM[comp] + m.PR_CPIGDP1[comp] * (
                        (m.LPCReboilerTemp + 273.15) - m.PR_RefTemp[comp]) + 2 * m.PR_CPIGDP2[comp] * \
                               m.PR_CPIGDP3[comp] * (-exp(2 / (m.LPCReboilerTemp + 273.15) * m.PR_CPIGDP3[comp]) + \
                                                     exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp])) / (
                                       exp(2 / (m.LPCReboilerTemp + 273.15) * m.PR_CPIGDP3[comp]) - 1) / (
                                       exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp]) - 1) - \
                               2 * m.PR_CPIGDP4[comp] * m.PR_CPIGDP5[comp] * (
                                       exp(2 / (m.LPCReboilerTemp + 273.15) * m.PR_CPIGDP5[comp]) - exp(
                                   2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp])) / \
                               (exp(2 / (m.LPCReboilerTemp + 273.15) * m.PR_CPIGDP5[comp]) + 1) / (
                                       exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp]) + 1))

        test_model.LPCReboilerHdlpH0 = Expression(test_model.Component, rule=LPCReboilerHdlpH0_Calc)

        def LPCReboilerHdlpM_Calc(m, comp1, comp2):
            return m.PR_m[comp1] * m.LPCReboilerHdlpSqrtTr[comp1] / (2 * sqrt(m.LPCReboilerHdlpalpha[comp1]))

        test_model.LPCReboilerHdlpM = Expression(test_model.Component, test_model.Component,
                                                 rule=LPCReboilerHdlpM_Calc)  # Calculate the enthalpy of liquid phase mixture %unit: J/mol

        def LPCReboilerHdlpHV_Calc(m):
            return sum([m.LPCReboilerVapLvMFrac[c] * m.LPCReboilerHdlpH0[c] for c in m.Component]) + \
                m.PR_RGas * (m.LPCReboilerTemp + 273.15) * ((m.LPCReboilerHdlpZV - 1) - 1 / (2 ** 1.5) * (
                        1 / m.LPCReboilerHdlpBV) * m.LPCReboilerHdlpZVdivBV * \
                                                            sum([sum([m.LPCReboilerVapLvMFrac[c1] *
                                                                      m.LPCReboilerVapLvMFrac[
                                                                          c2] * m.LPCReboilerHdlpAV * (
                                                                              1 + m.LPCReboilerHdlpM[c2, c1] + \
                                                                              m.LPCReboilerHdlpM[c1, c2]) for c1 in
                                                                      m.Component]) for c2 in m.Component]))

        test_model.LPCReboilerVapLvMEtlp = Expression(rule=LPCReboilerHdlpHV_Calc)

        # Calculate the enthalpy of liquid phase mixture %unit: J/mol
        def LPCReboilerHdlpHL_Calc(m):
            return sum([m.LPCReboilerLiqLvMFrac[c] * m.LPCReboilerHdlpH0[c] for c in m.Component]) + \
                m.PR_RGas * (m.LPCReboilerTemp + 273.15) * ((m.LPCReboilerHdlpZL - 1) - 1 / (2 ** 1.5) * (
                        1 / m.LPCReboilerHdlpBL) * m.LPCReboilerHdlpZLdivBL * \
                                                            sum([sum([m.LPCReboilerLiqLvMFrac[c1] *
                                                                      m.LPCReboilerLiqLvMFrac[
                                                                          c2] * m.LPCReboilerHdlpAL * (
                                                                              1 + m.LPCReboilerHdlpM[c2, c1] + \
                                                                              m.LPCReboilerHdlpM[c1, c2]) for c1 in
                                                                      m.Component]) for c2 in m.Component]))

        test_model.LPCReboilerLiqLvMEtlp = Expression(rule=LPCReboilerHdlpHL_Calc)

        # --------ASC Component Flowrate---------
        def ASCVapLvMCompFlow(m, tray, comp):
            return m.ASCVapLvMFrac[tray, comp] * m.ASCVapLvMFlow[tray]

        test_model.ASCVapLvMCompFlow = Expression(test_model.ASCTrays, test_model.Component, rule=ASCVapLvMCompFlow)

        def ASCLiqLvMCompFlow(m, tray, comp):
            return m.ASCLiqLvMFrac[tray, comp] * m.ASCLiqLvMFlow[tray]

        test_model.ASCLiqLvMCompFlow = Expression(test_model.ASCTrays, test_model.Component, rule=ASCLiqLvMCompFlow)

        # --------ASC Thermo and Enthalpy---------
        def ASCAllTraysTr_Calc(m, tray, comp):
            return (m.ASCTrayTemp[tray] + 273.15) / m.SRK_CrtcTemp[comp]

        test_model.ASCAllTraysTr = Expression(test_model.ASCTrays, test_model.Component, rule=ASCAllTraysTr_Calc)

        def ASCAllTraysalpha_Calc(m, tray, comp):
            return (1 + m.SRK_m[comp] * (1 - m.ASCAllTraysSqrtTr[tray, comp])) ** 2

        test_model.ASCAllTraysalpha = Expression(test_model.ASCTrays, test_model.Component, rule=ASCAllTraysalpha_Calc)

        def ASCAllTraysaij_Calc(m, tray, comp1, comp2):
            return m.ASCAllTraysSqrt_ai[tray, comp1, comp2] * (1 - m.SRK_Kij[comp1, comp2])

        test_model.ASCAllTraysaij = Expression(test_model.ASCTrays, test_model.Component, test_model.Component,
                                               rule=ASCAllTraysaij_Calc)

        def ASCAllTraysam_Calc(m, tray):
            return sum(
                [sum([m.ASCAllTraysaij[tray, c1, c2] * m.ASCLiqLvMFrac[tray, c1] * m.ASCLiqLvMFrac[tray, c2] for c1
                      in m.Component]) for c2 in m.Component])

        test_model.ASCAllTraysam = Expression(test_model.ASCTrays, rule=ASCAllTraysam_Calc)

        def ASCAllTraysbm_Calc(m, tray):
            return sum([m.SRK_bi[c] * m.ASCLiqLvMFrac[tray, c] for c in m.Component])

        test_model.ASCAllTraysbm = Expression(test_model.ASCTrays, rule=ASCAllTraysbm_Calc)

        def ASCAllTraysA_Calc(m, tray):
            return m.ASCAllTraysam[tray] * m.ASCTrayPressure[tray] * 1000 / (
                    (m.SRK_RGas * (m.ASCTrayTemp[tray] + 273.15)) ** 2)

        test_model.ASCAllTraysA = Expression(test_model.ASCTrays, rule=ASCAllTraysA_Calc)

        def ASCAllTraysB_Calc(m, tray):
            return m.ASCAllTraysbm[tray] * m.ASCTrayPressure[tray] * 1000 / (
                        m.SRK_RGas * (m.ASCTrayTemp[tray] + 273.15))

        test_model.ASCAllTraysB = Expression(test_model.ASCTrays, rule=ASCAllTraysB_Calc)

        def ASCAllTraysa1_Calc(m, tray):
            return m.ASCAllTraysA[tray] - m.ASCAllTraysB[tray] - m.ASCAllTraysB[tray] * m.ASCAllTraysB[tray]

        test_model.ASCAllTraysa1 = Expression(test_model.ASCTrays, rule=ASCAllTraysa1_Calc)

        def ASCAllTraysa2_Calc(m, tray):
            return -m.ASCAllTraysA[tray] * m.ASCAllTraysB[tray]

        test_model.ASCAllTraysa2 = Expression(test_model.ASCTrays, rule=ASCAllTraysa2_Calc)

        def ASCAllTraysBi_Calc(m, tray, comp):
            return m.SRK_bi[comp] * m.ASCTrayPressure[tray] * 1000 / (m.SRK_RGas * (m.ASCTrayTemp[tray] + 273.15))

        test_model.ASCAllTraysBi = Expression(test_model.ASCTrays, test_model.Component, rule=ASCAllTraysBi_Calc)

        def ASCAllTraysAijL_Calc(m, tray, comp):
            return sum([m.ASCAllTraysaij[tray, comp, c] * m.ASCLiqLvMFrac[tray, c] for c in m.Component]) * \
                m.ASCTrayPressure[tray] * 1000 / ((m.SRK_RGas * (m.ASCTrayTemp[tray] + 273.15)) ** 2)

        test_model.ASCAllTraysAijL = Expression(test_model.ASCTrays, test_model.Component, rule=ASCAllTraysAijL_Calc)

        def ASCAllTraysAijV_Calc(m, tray, comp):
            return sum([m.ASCAllTraysaij[tray, comp, c] * m.ASCVapLvMFrac[tray, c] for c in m.Component]) * \
                m.ASCTrayPressure[tray] * 1000 / ((m.SRK_RGas * (m.ASCTrayTemp[tray] + 273.15)) ** 2)

        test_model.ASCAllTraysAijV = Expression(test_model.ASCTrays, test_model.Component, rule=ASCAllTraysAijV_Calc)

        def ASCAllTraysPhiL_Calc(m, tray, comp):
            return exp(m.ASCAllTraysBi[tray, comp] / m.ASCAllTraysB[tray] * (m.ASCAllTraysZL[tray] - 1) - \
                       m.ASCAllTraysZLB[tray] - m.ASCAllTraysA[tray] / m.ASCAllTraysB[tray] * (2 * \
                                                                                               m.ASCAllTraysAijL[
                                                                                                   tray, comp] /
                                                                                               m.ASCAllTraysA[tray] -
                                                                                               m.ASCAllTraysBi[
                                                                                                   tray, comp] /
                                                                                               m.ASCAllTraysB[tray]) * \
                       m.ASCAllTraysBZL[tray])

        test_model.ASCAllTraysPhiL = Expression(test_model.ASCTrays, test_model.Component, rule=ASCAllTraysPhiL_Calc)

        def ASCAllTraysPhiV_Calc(m, tray, comp):
            return exp(m.ASCAllTraysBi[tray, comp] / m.ASCAllTraysB[tray] * (m.ASCAllTraysZV[tray] - 1) - \
                       m.ASCAllTraysZVB[tray] - m.ASCAllTraysA[tray] / m.ASCAllTraysB[tray] * (2 * \
                                                                                               m.ASCAllTraysAijV[
                                                                                                   tray, comp] /
                                                                                               m.ASCAllTraysA[tray] -
                                                                                               m.ASCAllTraysBi[
                                                                                                   tray, comp] /
                                                                                               m.ASCAllTraysB[tray]) * \
                       m.ASCAllTraysBZV[tray])

        test_model.ASCAllTraysPhiV = Expression(test_model.ASCTrays, test_model.Component, rule=ASCAllTraysPhiV_Calc)

        def ASCAllTraysH0_Calc(m, tray, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.ASCTrayTemp[tray] + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (m.SRK_DHFORM[comp] + m.SRK_CPIGDP1[comp] * (
                        (m.ASCTrayTemp[tray] + 273.15) - m.SRK_RefTemp[comp]) + 2 * m.SRK_CPIGDP2[comp] * \
                               m.SRK_CPIGDP3[comp] * (-exp(2 / (m.ASCTrayTemp[tray] + 273.15) * m.SRK_CPIGDP3[comp]) + \
                                                      exp(2 / m.SRK_RefTemp[comp] * m.SRK_CPIGDP3[comp])) / (
                                       exp(2 / (m.ASCTrayTemp[tray] + 273.15) * m.SRK_CPIGDP3[comp]) - 1) / (
                                       exp(2 / m.SRK_RefTemp[comp] * m.SRK_CPIGDP3[comp]) - 1) - \
                               2 * m.SRK_CPIGDP4[comp] * m.SRK_CPIGDP5[comp] * (
                                       exp(2 / (m.ASCTrayTemp[tray] + 273.15) * m.SRK_CPIGDP5[comp]) - exp(
                                   2 / m.SRK_RefTemp[comp] * m.SRK_CPIGDP5[comp])) / \
                               (exp(2 / (m.ASCTrayTemp[tray] + 273.15) * m.SRK_CPIGDP5[comp]) + 1) / (
                                       exp(2 / m.SRK_RefTemp[comp] * m.SRK_CPIGDP5[comp]) + 1))

        test_model.ASCAllTraysH0 = Expression(test_model.ASCTrays, test_model.Component, rule=ASCAllTraysH0_Calc)

        # Calculate the enthalpy of liquid phase mixture %unit: J/mol
        def ASCAllTraysHV_Calc(m, tray):
            return sum([m.ASCVapLvMFrac[tray, c] * m.ASCAllTraysH0[tray, c] for c in m.Component]) + \
                m.SRK_RGas * (m.ASCTrayTemp[tray] + 273.15) * (
                        (m.ASCAllTraysZV[tray] - 1) - m.ASCAllTraysBZV[tray] * (m.ASCAllTraysA[tray] / \
                                                                                m.ASCAllTraysB[tray] +
                                                                                m.ASCAllTrayssqrtam[tray] / (
                                                                                        m.SRK_RGas * (m.ASCTrayTemp[
                                                                                                          tray] + 273.15) *
                                                                                        m.ASCAllTraysbm[tray]) *
                                                                                sum([m.ASCVapLvMFrac[tray, c] * m.SRK_m[
                                                                                    c] * m.ASCAllTrayssqrtaiTr[tray, c]
                                                                                     for c in m.Component])))

        test_model.ASCVapLvMEtlp = Expression(test_model.ASCTrays, rule=ASCAllTraysHV_Calc)

        # Calculate the enthalpy of liquid phase mixture %unit: J/mol
        def ASCAllTraysHL_Calc(m, tray):
            return sum([m.ASCLiqLvMFrac[tray, c] * m.ASCAllTraysH0[tray, c] for c in m.Component]) + \
                m.SRK_RGas * (m.ASCTrayTemp[tray] + 273.15) * (
                        (m.ASCAllTraysZL[tray] - 1) - m.ASCAllTraysBZL[tray] * (m.ASCAllTraysA[tray] / \
                                                                                m.ASCAllTraysB[tray] +
                                                                                m.ASCAllTrayssqrtam[tray] / (
                                                                                        m.SRK_RGas * (m.ASCTrayTemp[
                                                                                                          tray] + 273.15) *
                                                                                        m.ASCAllTraysbm[tray]) *
                                                                                sum([m.ASCLiqLvMFrac[tray, c] * m.SRK_m[
                                                                                    c] * m.ASCAllTrayssqrtaiTr[tray, c]
                                                                                     for c in m.Component])))

        test_model.ASCLiqLvMEtlp = Expression(test_model.ASCTrays, rule=ASCAllTraysHL_Calc)

        # --------ASCCond Enthalpy---------
        def ASCCondSelfTr_Calc(m, comp):
            return (m.ASCCondTemp + 273.15) / m.SRK_CrtcTemp[comp]

        test_model.ASCCondSelfTr = Expression(test_model.Component, rule=ASCCondSelfTr_Calc)

        def ASCCondSelfalpha_Calc(m, comp):
            return (1 + m.SRK_m[comp] * (1 - m.ASCCondSelfSqrtTr[comp])) ** 2

        test_model.ASCCondSelfalpha = Expression(test_model.Component, rule=ASCCondSelfalpha_Calc)

        def ASCCondSelfaij_Calc(m, comp1, comp2):
            return m.ASCCondSelfSqrt_ai[comp1, comp2] * (1 - m.SRK_Kij[comp1, comp2])

        test_model.ASCCondSelfaij = Expression(test_model.Component, test_model.Component, rule=ASCCondSelfaij_Calc)

        def ASCCondSelfam_Calc(m):
            return sum(
                [sum([m.ASCCondSelfaij[c1, c2] * m.ASCCondLiqMFrac[c1] * m.ASCCondLiqMFrac[c2] for c1 in m.Component])
                 for
                 c2 in m.Component])

        test_model.ASCCondSelfam = Expression(rule=ASCCondSelfam_Calc)

        def ASCCondSelfbm_Calc(m):
            return sum([m.SRK_bi[c] * m.ASCCondLiqMFrac[c] for c in m.Component])

        test_model.ASCCondSelfbm = Expression(rule=ASCCondSelfbm_Calc)

        def ASCCondSelfA_Calc(m):
            return m.ASCCondSelfam * m.ASCCondPressure * 1000 / ((m.SRK_RGas * (m.ASCCondTemp + 273.15)) ** 2)

        test_model.ASCCondSelfA = Expression(rule=ASCCondSelfA_Calc)

        def ASCCondSelfB_Calc(m):
            return m.ASCCondSelfbm * m.ASCCondPressure * 1000 / (m.SRK_RGas * (m.ASCCondTemp + 273.15))

        test_model.ASCCondSelfB = Expression(rule=ASCCondSelfB_Calc)

        def ASCCondSelfa1_Calc(m):
            return m.ASCCondSelfA - m.ASCCondSelfB - m.ASCCondSelfB * m.ASCCondSelfB

        test_model.ASCCondSelfa1 = Expression(rule=ASCCondSelfa1_Calc)

        def ASCCondSelfa2_Calc(m):
            return -m.ASCCondSelfA * m.ASCCondSelfB

        test_model.ASCCondSelfa2 = Expression(rule=ASCCondSelfa2_Calc)

        def ASCCondSelfBi_Calc(m, comp):
            return m.SRK_bi[comp] * m.ASCCondPressure * 1000 / (m.SRK_RGas * (m.ASCCondTemp + 273.15))

        test_model.ASCCondSelfBi = Expression(test_model.Component, rule=ASCCondSelfBi_Calc)

        def ASCCondSelfAijL_Calc(m, comp):
            return sum(
                [m.ASCCondSelfaij[comp, c] * m.ASCCondLiqMFrac[c] for c in m.Component]) * m.ASCCondPressure * 1000 / (
                    (m.SRK_RGas * (m.ASCCondTemp + 273.15)) ** 2)

        test_model.ASCCondSelfAijL = Expression(test_model.Component, rule=ASCCondSelfAijL_Calc)

        def ASCCondSelfAijV_Calc(m, comp):
            return sum(
                [m.ASCCondSelfaij[comp, c] * m.ASCCondVapMFrac[c] for c in m.Component]) * m.ASCCondPressure * 1000 / (
                    (m.SRK_RGas * (m.ASCCondTemp + 273.15)) ** 2)

        test_model.ASCCondSelfAijV = Expression(test_model.Component, rule=ASCCondSelfAijV_Calc)

        def ASCCondSelfPhiL_Calc(m, comp):
            return exp(m.ASCCondSelfBi[comp] / m.ASCCondSelfB * (m.ASCCondSelfZL - 1) - \
                       m.ASCCondSelfZLB - m.ASCCondSelfA / m.ASCCondSelfB * (2 * \
                                                                             m.ASCCondSelfAijL[comp] / m.ASCCondSelfA -
                                                                             m.ASCCondSelfBi[comp] / m.ASCCondSelfB) * \
                       m.ASCCondSelfBZL)

        test_model.ASCCondSelfPhiL = Expression(test_model.Component, rule=ASCCondSelfPhiL_Calc)

        def ASCCondSelfPhiV_Calc(m, comp):
            return exp(m.ASCCondSelfBi[comp] / m.ASCCondSelfB * (m.ASCCondSelfZV - 1) - \
                       m.ASCCondSelfZVB - m.ASCCondSelfA / m.ASCCondSelfB * (2 * \
                                                                             m.ASCCondSelfAijV[comp] / m.ASCCondSelfA -
                                                                             m.ASCCondSelfBi[comp] / m.ASCCondSelfB) * \
                       m.ASCCondSelfBZV)

        test_model.ASCCondSelfPhiV = Expression(test_model.Component, rule=ASCCondSelfPhiV_Calc)

        def ASCCondSelfH0_Calc(m, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.ASCCondTemp + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (m.SRK_DHFORM[comp] + m.SRK_CPIGDP1[comp] * (
                        (m.ASCCondTemp + 273.15) - m.SRK_RefTemp[comp]) + 2 * m.SRK_CPIGDP2[comp] * \
                               m.SRK_CPIGDP3[comp] * (-exp(2 / (m.ASCCondTemp + 273.15) * m.SRK_CPIGDP3[comp]) + \
                                                      exp(2 / m.SRK_RefTemp[comp] * m.SRK_CPIGDP3[comp])) / (
                                       exp(2 / (m.ASCCondTemp + 273.15) * m.SRK_CPIGDP3[comp]) - 1) / (
                                       exp(2 / m.SRK_RefTemp[comp] * m.SRK_CPIGDP3[comp]) - 1) - \
                               2 * m.SRK_CPIGDP4[comp] * m.SRK_CPIGDP5[comp] * (
                                       exp(2 / (m.ASCCondTemp + 273.15) * m.SRK_CPIGDP5[comp]) - exp(
                                   2 / m.SRK_RefTemp[comp] * m.SRK_CPIGDP5[comp])) / \
                               (exp(2 / (m.ASCCondTemp + 273.15) * m.SRK_CPIGDP5[comp]) + 1) / (
                                       exp(2 / m.SRK_RefTemp[comp] * m.SRK_CPIGDP5[comp]) + 1))

        test_model.ASCCondSelfH0 = Expression(test_model.Component, rule=ASCCondSelfH0_Calc)

        # Calculate the enthalpy of liquid phase mixture %unit: J/mol
        def ASCCondSelfHV_Calc(m):
            return sum([m.ASCCondVapMFrac[c] * m.ASCCondSelfH0[c] for c in m.Component]) + \
                m.SRK_RGas * (m.ASCCondTemp + 273.15) * ((m.ASCCondSelfZV - 1) - m.ASCCondSelfBZV * (m.ASCCondSelfA / \
                                                                                                     m.ASCCondSelfB + m.ASCCondSelfsqrtam / (
                                                                                                             m.SRK_RGas * (
                                                                                                             m.ASCCondTemp + 273.15) * m.ASCCondSelfbm) *
                                                                                                     sum([
                                                                                                             m.ASCCondVapMFrac[
                                                                                                                 c] *
                                                                                                             m.SRK_m[
                                                                                                                 c] *
                                                                                                             m.ASCCondSelfsqrtaiTr[
                                                                                                                 c] for
                                                                                                             c in
                                                                                                             m.Component])))

        test_model.ASCCondVapMEtlp = Expression(rule=ASCCondSelfHV_Calc)

        # Calculate the enthalpy of liquid phase mixture %unit: J/mol
        def ASCCondSelfHL_Calc(m):
            return sum([m.ASCCondLiqMFrac[c] * m.ASCCondSelfH0[c] for c in m.Component]) + \
                m.SRK_RGas * (m.ASCCondTemp + 273.15) * ((m.ASCCondSelfZL - 1) - m.ASCCondSelfBZL * (m.ASCCondSelfA / \
                                                                                                     m.ASCCondSelfB + m.ASCCondSelfsqrtam / (
                                                                                                             m.SRK_RGas * (
                                                                                                             m.ASCCondTemp + 273.15) * m.ASCCondSelfbm) *
                                                                                                     sum([
                                                                                                             m.ASCCondLiqMFrac[
                                                                                                                 c] *
                                                                                                             m.SRK_m[
                                                                                                                 c] *
                                                                                                             m.ASCCondSelfsqrtaiTr[
                                                                                                                 c] for
                                                                                                             c in
                                                                                                             m.Component])))

        test_model.ASCCondLiqMEtlp = Expression(rule=ASCCondSelfHL_Calc)

        # --------HeatGan Enthalpy---------
        def HeatGanSelfTr_Calc(m, comp):
            return (m.HeatGanTemp + 273.15) / m.PR_CrtcTemp[comp]

        test_model.HeatGanSelfTr = Expression(test_model.Component, rule=HeatGanSelfTr_Calc)

        def HeatGanSelfalpha_Calc(m, comp):
            return (1 + m.PR_m[comp] * (1 - m.HeatGanSelfSqrtTr[comp])) ** 2

        test_model.HeatGanSelfalpha = Expression(test_model.Component, rule=HeatGanSelfalpha_Calc)

        def HeatGanSelfai_Calc(m, comp):
            return 0.45724 * m.HeatGanSelfalpha[comp] * (m.PR_RGas ** 2) * (m.PR_CrtcTemp[comp] ** 2) / \
                m.PR_CrtcPressure[
                    comp]

        test_model.HeatGanSelfai = Expression(test_model.Component, rule=HeatGanSelfai_Calc)

        def HeatGanSelfaij_Calc(m, comp1, comp2):
            return m.HeatGanSelfSqrt_ai[comp1, comp2] * (1 - m.PR_Kij[comp1, comp2])

        test_model.HeatGanSelfaij = Expression(test_model.Component, test_model.Component, rule=HeatGanSelfaij_Calc)

        def HeatGanSelfaV_Calc(m):
            return sum(
                [sum([m.HeatGanSelfaij[c1, c2] * m.LPCVapLvMFrac[52, c1] * m.LPCVapLvMFrac[52, c2] for c1 in
                      m.Component])
                 for c2 in m.Component])

        test_model.HeatGanSelfaV = Expression(rule=HeatGanSelfaV_Calc)

        def HeatGanSelfbV_Calc(m):
            return sum([m.PR_bi[c] * m.LPCVapLvMFrac[52, c] for c in m.Component])

        test_model.HeatGanSelfbV = Expression(rule=HeatGanSelfbV_Calc)

        def HeatGanSelfAV_Calc(m):
            return m.HeatGanSelfaV * m.HeatGanPressure * 1000 / ((m.PR_RGas * (m.HeatGanTemp + 273.15)) ** 2)

        test_model.HeatGanSelfAV = Expression(rule=HeatGanSelfAV_Calc)

        def HeatGanSelfBV_Calc(m):
            return m.HeatGanSelfbV * m.HeatGanPressure * 1000 / (m.PR_RGas * (m.HeatGanTemp + 273.15))

        test_model.HeatGanSelfBV = Expression(rule=HeatGanSelfBV_Calc)

        def HeatGanSelfb1_Calc(m):
            return m.HeatGanSelfBV - 1

        test_model.HeatGanSelfb1 = Expression(rule=HeatGanSelfb1_Calc)

        def HeatGanSelfb2_Calc(m):
            return m.HeatGanSelfAV - 3 * m.HeatGanSelfBV ** 2 - 2 * m.HeatGanSelfBV

        test_model.HeatGanSelfb2 = Expression(rule=HeatGanSelfb2_Calc)

        def HeatGanSelfb3_Calc(m):
            return m.HeatGanSelfBV ** 2 + m.HeatGanSelfBV ** 3 - m.HeatGanSelfAV * m.HeatGanSelfBV

        test_model.HeatGanSelfb3 = Expression(rule=HeatGanSelfb3_Calc)

        def HeatGanSelfH0_Calc(m, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.HeatGanTemp + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (
                        m.PR_DHFORM[comp] + m.PR_CPIGDP1[comp] * ((m.HeatGanTemp + 273.15) - m.PR_RefTemp[comp]) + 2 *
                        m.PR_CPIGDP2[comp] * \
                        m.PR_CPIGDP3[comp] * (-exp(2 / (m.HeatGanTemp + 273.15) * m.PR_CPIGDP3[comp]) + \
                                              exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp])) / (
                                exp(2 / (m.HeatGanTemp + 273.15) * m.PR_CPIGDP3[comp]) - 1) / (
                                exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp]) - 1) - \
                        2 * m.PR_CPIGDP4[comp] * m.PR_CPIGDP5[comp] * (
                                exp(2 / (m.HeatGanTemp + 273.15) * m.PR_CPIGDP5[comp]) - exp(
                            2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp])) / \
                        (exp(2 / (m.HeatGanTemp + 273.15) * m.PR_CPIGDP5[comp]) + 1) / (
                                exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp]) + 1))

        test_model.HeatGanSelfH0 = Expression(test_model.Component, rule=HeatGanSelfH0_Calc)

        def HeatGanSelfM_Calc(m, comp1, comp2):
            return m.PR_m[comp1] * m.HeatGanSelfSqrtTr[comp1] / (2 * sqrt(m.HeatGanSelfalpha[comp1]))

        test_model.HeatGanSelfM = Expression(test_model.Component, test_model.Component,
                                             rule=HeatGanSelfM_Calc)  # Calculate the enthalpy of liquid phase mixture %unit: J/mol

        def HeatGanSelfHV_Calc(m):
            return sum([m.LPCVapLvMFrac[52, c] * m.HeatGanSelfH0[c] for c in m.Component]) + \
                m.PR_RGas * (m.HeatGanTemp + 273.15) * (
                        (m.HeatGanSelfZV - 1) - 1 / (2 ** 1.5) * (1 / m.HeatGanSelfBV) * m.HeatGanSelfZVdivBV * \
                        sum([sum([m.LPCVapLvMFrac[52, c1] * m.LPCVapLvMFrac[52, c2] * m.HeatGanSelfAV * (
                                1 + m.HeatGanSelfM[c2, c1] + \
                                m.HeatGanSelfM[c1, c2]) for c1 in m.Component]) for c2 in m.Component]))

        test_model.HeatGanMEtlp = Expression(rule=HeatGanSelfHV_Calc)

        # --------RecoGan Enthalpy---------
        def RecoGanSelfTr_Calc(m, comp):
            return (m.RecoGanTemp + 273.15) / m.PR_CrtcTemp[comp]

        test_model.RecoGanSelfTr = Expression(test_model.Component, rule=RecoGanSelfTr_Calc)

        def RecoGanSelfalpha_Calc(m, comp):
            return (1 + m.PR_m[comp] * (1 - m.RecoGanSelfSqrtTr[comp])) ** 2

        test_model.RecoGanSelfalpha = Expression(test_model.Component, rule=RecoGanSelfalpha_Calc)

        def RecoGanSelfai_Calc(m, comp):
            return 0.45724 * m.RecoGanSelfalpha[comp] * (m.PR_RGas ** 2) * (m.PR_CrtcTemp[comp] ** 2) / \
                m.PR_CrtcPressure[
                    comp]

        test_model.RecoGanSelfai = Expression(test_model.Component, rule=RecoGanSelfai_Calc)

        def RecoGanSelfaij_Calc(m, comp1, comp2):
            return m.RecoGanSelfSqrt_ai[comp1, comp2] * (1 - m.PR_Kij[comp1, comp2])

        test_model.RecoGanSelfaij = Expression(test_model.Component, test_model.Component, rule=RecoGanSelfaij_Calc)

        def RecoGanSelfaV_Calc(m):
            return sum(
                [sum([m.RecoGanSelfaij[c1, c2] * m.LPCVapLvMFrac[52, c1] * m.LPCVapLvMFrac[52, c2] for c1 in
                      m.Component])
                 for c2 in m.Component])

        test_model.RecoGanSelfaV = Expression(rule=RecoGanSelfaV_Calc)

        def RecoGanSelfbV_Calc(m):
            return sum([m.PR_bi[c] * m.LPCVapLvMFrac[52, c] for c in m.Component])

        test_model.RecoGanSelfbV = Expression(rule=RecoGanSelfbV_Calc)

        def RecoGanSelfAV_Calc(m):
            return m.RecoGanSelfaV * m.RecoGanPressure * 1000 / ((m.PR_RGas * (m.RecoGanTemp + 273.15)) ** 2)

        test_model.RecoGanSelfAV = Expression(rule=RecoGanSelfAV_Calc)

        def RecoGanSelfBV_Calc(m):
            return m.RecoGanSelfbV * m.RecoGanPressure * 1000 / (m.PR_RGas * (m.RecoGanTemp + 273.15))

        test_model.RecoGanSelfBV = Expression(rule=RecoGanSelfBV_Calc)

        def RecoGanSelfb1_Calc(m):
            return m.RecoGanSelfBV - 1

        test_model.RecoGanSelfb1 = Expression(rule=RecoGanSelfb1_Calc)

        def RecoGanSelfb2_Calc(m):
            return m.RecoGanSelfAV - 3 * m.RecoGanSelfBV ** 2 - 2 * m.RecoGanSelfBV

        test_model.RecoGanSelfb2 = Expression(rule=RecoGanSelfb2_Calc)

        def RecoGanSelfb3_Calc(m):
            return m.RecoGanSelfBV ** 2 + m.RecoGanSelfBV ** 3 - m.RecoGanSelfAV * m.RecoGanSelfBV

        test_model.RecoGanSelfb3 = Expression(rule=RecoGanSelfb3_Calc)

        def RecoGanSelfH0_Calc(m, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.RecoGanTemp + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (
                        m.PR_DHFORM[comp] + m.PR_CPIGDP1[comp] * ((m.RecoGanTemp + 273.15) - m.PR_RefTemp[comp]) + 2 *
                        m.PR_CPIGDP2[comp] * \
                        m.PR_CPIGDP3[comp] * (-exp(2 / (m.RecoGanTemp + 273.15) * m.PR_CPIGDP3[comp]) + \
                                              exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp])) / (
                                exp(2 / (m.RecoGanTemp + 273.15) * m.PR_CPIGDP3[comp]) - 1) / (
                                exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp]) - 1) - \
                        2 * m.PR_CPIGDP4[comp] * m.PR_CPIGDP5[comp] * (
                                exp(2 / (m.RecoGanTemp + 273.15) * m.PR_CPIGDP5[comp]) - exp(
                            2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp])) / \
                        (exp(2 / (m.RecoGanTemp + 273.15) * m.PR_CPIGDP5[comp]) + 1) / (
                                exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp]) + 1))

        test_model.RecoGanSelfH0 = Expression(test_model.Component, rule=RecoGanSelfH0_Calc)

        def RecoGanSelfM_Calc(m, comp1, comp2):
            return m.PR_m[comp1] * m.RecoGanSelfSqrtTr[comp1] / (2 * sqrt(m.RecoGanSelfalpha[comp1]))

        test_model.RecoGanSelfM = Expression(test_model.Component, test_model.Component,
                                             rule=RecoGanSelfM_Calc)  # Calculate the enthalpy of liquid phase mixture %unit: J/mol

        def RecoGanSelfHV_Calc(m):
            return sum([m.LPCVapLvMFrac[52, c] * m.RecoGanSelfH0[c] for c in m.Component]) + \
                m.PR_RGas * (m.RecoGanTemp + 273.15) * (
                        (m.RecoGanSelfZV - 1) - 1 / (2 ** 1.5) * (1 / m.RecoGanSelfBV) * m.RecoGanSelfZVdivBV * \
                        sum([sum([m.LPCVapLvMFrac[52, c1] * m.LPCVapLvMFrac[52, c2] * m.RecoGanSelfAV * (
                                1 + m.RecoGanSelfM[c2, c1] + \
                                m.RecoGanSelfM[c1, c2]) for c1 in m.Component]) for c2 in m.Component]))

        test_model.RecoGanMEtlp = Expression(rule=RecoGanSelfHV_Calc)

        # --------RecoGar Enthalpy---------
        def RecoGarSelfTr_Calc(m, comp):
            return (m.RecoGarTemp + 273.15) / m.PR_CrtcTemp[comp]

        test_model.RecoGarSelfTr = Expression(test_model.Component, rule=RecoGarSelfTr_Calc)

        def RecoGarSelfalpha_Calc(m, comp):
            return (1 + m.PR_m[comp] * (1 - m.RecoGarSelfSqrtTr[comp])) ** 2

        test_model.RecoGarSelfalpha = Expression(test_model.Component, rule=RecoGarSelfalpha_Calc)

        def RecoGarSelfai_Calc(m, comp):
            return 0.45724 * m.RecoGarSelfalpha[comp] * (m.PR_RGas ** 2) * (m.PR_CrtcTemp[comp] ** 2) / \
                m.PR_CrtcPressure[
                    comp]

        test_model.RecoGarSelfai = Expression(test_model.Component, rule=RecoGarSelfai_Calc)

        def RecoGarSelfaij_Calc(m, comp1, comp2):
            return m.RecoGarSelfSqrt_ai[comp1, comp2] * (1 - m.PR_Kij[comp1, comp2])

        test_model.RecoGarSelfaij = Expression(test_model.Component, test_model.Component, rule=RecoGarSelfaij_Calc)

        def RecoGarSelfaV_Calc(m):
            return sum(
                [sum([m.RecoGarSelfaij[c1, c2] * m.ASCCondVapMFrac[c1] * m.ASCCondVapMFrac[c2] for c1 in m.Component])
                 for
                 c2 in m.Component])

        test_model.RecoGarSelfaV = Expression(rule=RecoGarSelfaV_Calc)

        def RecoGarSelfbV_Calc(m):
            return sum([m.PR_bi[c] * m.ASCCondVapMFrac[c] for c in m.Component])

        test_model.RecoGarSelfbV = Expression(rule=RecoGarSelfbV_Calc)

        def RecoGarSelfAV_Calc(m):
            return m.RecoGarSelfaV * m.RecoGarPressure * 1000 / ((m.PR_RGas * (m.RecoGarTemp + 273.15)) ** 2)

        test_model.RecoGarSelfAV = Expression(rule=RecoGarSelfAV_Calc)

        def RecoGarSelfBV_Calc(m):
            return m.RecoGarSelfbV * m.RecoGarPressure * 1000 / (m.PR_RGas * (m.RecoGarTemp + 273.15))

        test_model.RecoGarSelfBV = Expression(rule=RecoGarSelfBV_Calc)

        def RecoGarSelfb1_Calc(m):
            return m.RecoGarSelfBV - 1

        test_model.RecoGarSelfb1 = Expression(rule=RecoGarSelfb1_Calc)

        def RecoGarSelfb2_Calc(m):
            return m.RecoGarSelfAV - 3 * m.RecoGarSelfBV ** 2 - 2 * m.RecoGarSelfBV

        test_model.RecoGarSelfb2 = Expression(rule=RecoGarSelfb2_Calc)

        def RecoGarSelfb3_Calc(m):
            return m.RecoGarSelfBV ** 2 + m.RecoGarSelfBV ** 3 - m.RecoGarSelfAV * m.RecoGarSelfBV

        test_model.RecoGarSelfb3 = Expression(rule=RecoGarSelfb3_Calc)

        def RecoGarSelfH0_Calc(m, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.RecoGarTemp + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (
                        m.PR_DHFORM[comp] + m.PR_CPIGDP1[comp] * ((m.RecoGarTemp + 273.15) - m.PR_RefTemp[comp]) + 2 *
                        m.PR_CPIGDP2[comp] * \
                        m.PR_CPIGDP3[comp] * (-exp(2 / (m.RecoGarTemp + 273.15) * m.PR_CPIGDP3[comp]) + \
                                              exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp])) / (
                                exp(2 / (m.RecoGarTemp + 273.15) * m.PR_CPIGDP3[comp]) - 1) / (
                                exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp]) - 1) - \
                        2 * m.PR_CPIGDP4[comp] * m.PR_CPIGDP5[comp] * (
                                exp(2 / (m.RecoGarTemp + 273.15) * m.PR_CPIGDP5[comp]) - exp(
                            2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp])) / \
                        (exp(2 / (m.RecoGarTemp + 273.15) * m.PR_CPIGDP5[comp]) + 1) / (
                                exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp]) + 1))

        test_model.RecoGarSelfH0 = Expression(test_model.Component, rule=RecoGarSelfH0_Calc)

        def RecoGarSelfM_Calc(m, comp1, comp2):
            return m.PR_m[comp1] * m.RecoGarSelfSqrtTr[comp1] / (2 * sqrt(m.RecoGarSelfalpha[comp1]))

        test_model.RecoGarSelfM = Expression(test_model.Component, test_model.Component,
                                             rule=RecoGarSelfM_Calc)  # Calculate the enthalpy of liquid phase mixture %unit: J/mol

        def RecoGarSelfHV_Calc(m):
            return sum([m.ASCCondVapMFrac[c] * m.RecoGarSelfH0[c] for c in m.Component]) + \
                m.PR_RGas * (m.RecoGarTemp + 273.15) * (
                        (m.RecoGarSelfZV - 1) - 1 / (2 ** 1.5) * (1 / m.RecoGarSelfBV) * m.RecoGarSelfZVdivBV * \
                        sum([sum([m.ASCCondVapMFrac[c1] * m.ASCCondVapMFrac[c2] * m.RecoGarSelfAV * (
                                1 + m.RecoGarSelfM[c2, c1] + \
                                m.RecoGarSelfM[c1, c2]) for c1 in m.Component]) for c2 in m.Component]))

        test_model.RecoGarMEtlp = Expression(rule=RecoGarSelfHV_Calc)

        # --------RecoWn Enthalpy---------
        def RecoWnSelfTr_Calc(m, comp):
            return (m.RecoWnTemp + 273.15) / m.PR_CrtcTemp[comp]

        test_model.RecoWnSelfTr = Expression(test_model.Component, rule=RecoWnSelfTr_Calc)

        def RecoWnSelfalpha_Calc(m, comp):
            return (1 + m.PR_m[comp] * (1 - m.RecoWnSelfSqrtTr[comp])) ** 2

        test_model.RecoWnSelfalpha = Expression(test_model.Component, rule=RecoWnSelfalpha_Calc)

        def RecoWnSelfai_Calc(m, comp):
            return 0.45724 * m.RecoWnSelfalpha[comp] * (m.PR_RGas ** 2) * (m.PR_CrtcTemp[comp] ** 2) / \
                m.PR_CrtcPressure[
                    comp]

        test_model.RecoWnSelfai = Expression(test_model.Component, rule=RecoWnSelfai_Calc)

        def RecoWnSelfaij_Calc(m, comp1, comp2):
            return m.RecoWnSelfSqrt_ai[comp1, comp2] * (1 - m.PR_Kij[comp1, comp2])

        test_model.RecoWnSelfaij = Expression(test_model.Component, test_model.Component, rule=RecoWnSelfaij_Calc)

        def RecoWnSelfaV_Calc(m):
            return sum(
                [sum([m.RecoWnSelfaij[c1, c2] * m.LPCExt46MFrac[c1] * m.LPCExt46MFrac[c2] for c1 in m.Component]) for c2
                 in
                 m.Component])

        test_model.RecoWnSelfaV = Expression(rule=RecoWnSelfaV_Calc)

        def RecoWnSelfbV_Calc(m):
            return sum([m.PR_bi[c] * m.LPCExt46MFrac[c] for c in m.Component])

        test_model.RecoWnSelfbV = Expression(rule=RecoWnSelfbV_Calc)

        def RecoWnSelfAV_Calc(m):
            return m.RecoWnSelfaV * m.RecoWnPressure * 1000 / ((m.PR_RGas * (m.RecoWnTemp + 273.15)) ** 2)

        test_model.RecoWnSelfAV = Expression(rule=RecoWnSelfAV_Calc)

        def RecoWnSelfBV_Calc(m):
            return m.RecoWnSelfbV * m.RecoWnPressure * 1000 / (m.PR_RGas * (m.RecoWnTemp + 273.15))

        test_model.RecoWnSelfBV = Expression(rule=RecoWnSelfBV_Calc)

        def RecoWnSelfb1_Calc(m):
            return m.RecoWnSelfBV - 1

        test_model.RecoWnSelfb1 = Expression(rule=RecoWnSelfb1_Calc)

        def RecoWnSelfb2_Calc(m):
            return m.RecoWnSelfAV - 3 * m.RecoWnSelfBV ** 2 - 2 * m.RecoWnSelfBV

        test_model.RecoWnSelfb2 = Expression(rule=RecoWnSelfb2_Calc)

        def RecoWnSelfb3_Calc(m):
            return m.RecoWnSelfBV ** 2 + m.RecoWnSelfBV ** 3 - m.RecoWnSelfAV * m.RecoWnSelfBV

        test_model.RecoWnSelfb3 = Expression(rule=RecoWnSelfb3_Calc)

        def RecoWnSelfH0_Calc(m, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.RecoWnTemp + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (m.PR_DHFORM[comp] + m.PR_CPIGDP1[comp] * (
                            (m.RecoWnTemp + 273.15) - m.PR_RefTemp[comp]) + 2 *
                               m.PR_CPIGDP2[comp] * \
                               m.PR_CPIGDP3[comp] * (-exp(2 / (m.RecoWnTemp + 273.15) * m.PR_CPIGDP3[comp]) + \
                                                     exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp])) / (
                                       exp(2 / (m.RecoWnTemp + 273.15) * m.PR_CPIGDP3[comp]) - 1) / (
                                       exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp]) - 1) - \
                               2 * m.PR_CPIGDP4[comp] * m.PR_CPIGDP5[comp] * (
                                       exp(2 / (m.RecoWnTemp + 273.15) * m.PR_CPIGDP5[comp]) - exp(
                                   2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp])) / \
                               (exp(2 / (m.RecoWnTemp + 273.15) * m.PR_CPIGDP5[comp]) + 1) / (
                                       exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp]) + 1))

        test_model.RecoWnSelfH0 = Expression(test_model.Component, rule=RecoWnSelfH0_Calc)

        def RecoWnSelfM_Calc(m, comp1, comp2):
            return m.PR_m[comp1] * m.RecoWnSelfSqrtTr[comp1] / (2 * sqrt(m.RecoWnSelfalpha[comp1]))

        test_model.RecoWnSelfM = Expression(test_model.Component, test_model.Component,
                                            rule=RecoWnSelfM_Calc)  # Calculate the enthalpy of liquid phase mixture %unit: J/mol

        def RecoWnSelfHV_Calc(m):
            return sum([m.LPCExt46MFrac[c] * m.RecoWnSelfH0[c] for c in m.Component]) + \
                m.PR_RGas * (m.RecoWnTemp + 273.15) * (
                        (m.RecoWnSelfZV - 1) - 1 / (2 ** 1.5) * (1 / m.RecoWnSelfBV) * m.RecoWnSelfZVdivBV * \
                        sum([sum([m.LPCExt46MFrac[c1] * m.LPCExt46MFrac[c2] * m.RecoWnSelfAV * (
                                1 + m.RecoWnSelfM[c2, c1] + \
                                m.RecoWnSelfM[c1, c2]) for c1 in m.Component]) for c2 in m.Component]))

        test_model.RecoWnMEtlp = Expression(rule=RecoWnSelfHV_Calc)

        # --------RecoGox Enthalpy---------
        def RecoGoxSelfTr_Calc(m, comp):
            return (m.RecoGoxTemp + 273.15) / m.PR_CrtcTemp[comp]

        test_model.RecoGoxSelfTr = Expression(test_model.Component, rule=RecoGoxSelfTr_Calc)

        def RecoGoxSelfalpha_Calc(m, comp):
            return (1 + m.PR_m[comp] * (1 - m.RecoGoxSelfSqrtTr[comp])) ** 2

        test_model.RecoGoxSelfalpha = Expression(test_model.Component, rule=RecoGoxSelfalpha_Calc)

        def RecoGoxSelfai_Calc(m, comp):
            return 0.45724 * m.RecoGoxSelfalpha[comp] * (m.PR_RGas ** 2) * (m.PR_CrtcTemp[comp] ** 2) / \
                m.PR_CrtcPressure[
                    comp]

        test_model.RecoGoxSelfai = Expression(test_model.Component, rule=RecoGoxSelfai_Calc)

        def RecoGoxSelfaij_Calc(m, comp1, comp2):
            return m.RecoGoxSelfSqrt_ai[comp1, comp2] * (1 - m.PR_Kij[comp1, comp2])

        test_model.RecoGoxSelfaij = Expression(test_model.Component, test_model.Component, rule=RecoGoxSelfaij_Calc)

        def RecoGoxSelfaV_Calc(m):
            return sum(
                [sum([m.RecoGoxSelfaij[c1, c2] * m.LPCReboilerLiqLvMFrac[c1] * m.LPCReboilerLiqLvMFrac[c2] for c1 in
                      m.Component]) for c2 in m.Component])

        test_model.RecoGoxSelfaV = Expression(rule=RecoGoxSelfaV_Calc)

        def RecoGoxSelfbV_Calc(m):
            return sum([m.PR_bi[c] * m.LPCReboilerLiqLvMFrac[c] for c in m.Component])

        test_model.RecoGoxSelfbV = Expression(rule=RecoGoxSelfbV_Calc)

        def RecoGoxSelfAV_Calc(m):
            return m.RecoGoxSelfaV * m.RecoGoxPressure * 1000 / ((m.PR_RGas * (m.RecoGoxTemp + 273.15)) ** 2)

        test_model.RecoGoxSelfAV = Expression(rule=RecoGoxSelfAV_Calc)

        def RecoGoxSelfBV_Calc(m):
            return m.RecoGoxSelfbV * m.RecoGoxPressure * 1000 / (m.PR_RGas * (m.RecoGoxTemp + 273.15))

        test_model.RecoGoxSelfBV = Expression(rule=RecoGoxSelfBV_Calc)

        def RecoGoxSelfb1_Calc(m):
            return m.RecoGoxSelfBV - 1

        test_model.RecoGoxSelfb1 = Expression(rule=RecoGoxSelfb1_Calc)

        def RecoGoxSelfb2_Calc(m):
            return m.RecoGoxSelfAV - 3 * m.RecoGoxSelfBV ** 2 - 2 * m.RecoGoxSelfBV

        test_model.RecoGoxSelfb2 = Expression(rule=RecoGoxSelfb2_Calc)

        def RecoGoxSelfb3_Calc(m):
            return m.RecoGoxSelfBV ** 2 + m.RecoGoxSelfBV ** 3 - m.RecoGoxSelfAV * m.RecoGoxSelfBV

        test_model.RecoGoxSelfb3 = Expression(rule=RecoGoxSelfb3_Calc)

        def RecoGoxSelfH0_Calc(m, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.RecoGoxTemp + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (
                        m.PR_DHFORM[comp] + m.PR_CPIGDP1[comp] * ((m.RecoGoxTemp + 273.15) - m.PR_RefTemp[comp]) + 2 *
                        m.PR_CPIGDP2[comp] * \
                        m.PR_CPIGDP3[comp] * (-exp(2 / (m.RecoGoxTemp + 273.15) * m.PR_CPIGDP3[comp]) + \
                                              exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp])) / (
                                exp(2 / (m.RecoGoxTemp + 273.15) * m.PR_CPIGDP3[comp]) - 1) / (
                                exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp]) - 1) - \
                        2 * m.PR_CPIGDP4[comp] * m.PR_CPIGDP5[comp] * (
                                exp(2 / (m.RecoGoxTemp + 273.15) * m.PR_CPIGDP5[comp]) - exp(
                            2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp])) / \
                        (exp(2 / (m.RecoGoxTemp + 273.15) * m.PR_CPIGDP5[comp]) + 1) / (
                                exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp]) + 1))

        test_model.RecoGoxSelfH0 = Expression(test_model.Component, rule=RecoGoxSelfH0_Calc)

        def RecoGoxSelfM_Calc(m, comp1, comp2):
            return m.PR_m[comp1] * m.RecoGoxSelfSqrtTr[comp1] / (2 * sqrt(m.RecoGoxSelfalpha[comp1]))

        test_model.RecoGoxSelfM = Expression(test_model.Component, test_model.Component,
                                             rule=RecoGoxSelfM_Calc)  # Calculate the enthalpy of liquid phase mixture %unit: J/mol

        def RecoGoxSelfHV_Calc(m):
            return sum([m.LPCReboilerLiqLvMFrac[c] * m.RecoGoxSelfH0[c] for c in m.Component]) + \
                m.PR_RGas * (m.RecoGoxTemp + 273.15) * (
                        (m.RecoGoxSelfZV - 1) - 1 / (2 ** 1.5) * (1 / m.RecoGoxSelfBV) * m.RecoGoxSelfZVdivBV * \
                        sum([sum([m.LPCReboilerLiqLvMFrac[c1] * m.LPCReboilerLiqLvMFrac[c2] * m.RecoGoxSelfAV * (
                                1 + m.RecoGoxSelfM[c2, c1] + \
                                m.RecoGoxSelfM[c1, c2]) for c1 in m.Component]) for c2 in m.Component]))

        test_model.RecoGoxMEtlp = Expression(rule=RecoGoxSelfHV_Calc)

        # ===================================
        #
        #         Constraints
        #
        # ===================================

        # ----------------------------------
        #           Feed
        # ----------------------------------

        # ----------------------------------
        #           FeedSplitter
        # ----------------------------------
        # --------FeedSplitter Mass Balance---------
        def FeedSplitterMassBlnc(m):
            return (
                    m.FeedMFlow - m.FeedSplitterOut1MFlow - m.FeedSplitterOut2MFlow - m.FeedSplitterOut3MFlow) * 0.001000 == 0

        test_model.FeedSplitterMassBlnc = Constraint(rule=FeedSplitterMassBlnc)

        def FeedSplitterOut1Spec(m):
            return (m.FeedSplitterOut1MFlow - m.FeedMFlow * m.FeedSplitterOut1Ratio) * 0.001000 == 0

        test_model.FeedSplitterOut1Spec = Constraint(rule=FeedSplitterOut1Spec)

        def FeedSplitterOut2Spec(m):
            return (m.FeedSplitterOut2MFlow - m.FeedMFlow * m.FeedSplitterOut2Ratio) * 0.001000 == 0

        test_model.FeedSplitterOut2Spec = Constraint(rule=FeedSplitterOut2Spec)

        # ----------------------------------
        #           TaCooler
        # ----------------------------------
        # --------TaCooler Mass Balance---------
        def TaCoolerMassBlnc(m, comp):
            return m.TaCoolerVapMFlow * m.TaCoolerVapMFrac[comp] + \
                m.TaCoolerLiqMFlow * m.TaCoolerLiqMFrac[comp] - m.FeedSplitterOut1MFlow * m.FeedMFrac[comp] == 0

        test_model.TaCoolerMassBlnc = Constraint(test_model.Component, rule=TaCoolerMassBlnc)

        # --------TaCooler Energy Balance---------
        def TaCoolerEngBlnc(m):
            return ((m.FeedMEtlp - m.TaCoolerMEtlp) * m.FeedSplitterOut1MFlow - m.TaCoolerMHeat) * 1.000000e-05 == 0

        test_model.TaCoolerEngBlnc = Constraint(rule=TaCoolerEngBlnc)

        # --------TaCooler Phase Equilibrium---------
        def TaCoolerSelfSqrtTr_Calc(m, comp):
            return m.TaCoolerSelfSqrtTr[comp] ** 2 - m.TaCoolerSelfTr[comp] == 0

        test_model.TaCoolerSelfSqrtTr_Calc = Constraint(test_model.Component, rule=TaCoolerSelfSqrtTr_Calc)

        def TaCoolerSelfSqrt_ai_Calc(m, comp1, comp2):
            return m.TaCoolerSelfSqrt_ai[comp1, comp2] ** 2 - m.TaCoolerSelfai[comp1] * m.TaCoolerSelfai[comp2] == 0

        test_model.TaCoolerSelfSqrt_ai_Calc = Constraint(test_model.Component, test_model.Component,
                                                         rule=TaCoolerSelfSqrt_ai_Calc)

        def TaCoolerSelfZL_Calc(m):
            return m.TaCoolerSelfZL ** 3 + m.TaCoolerSelfa1 * m.TaCoolerSelfZL ** 2 + \
                m.TaCoolerSelfa2 * m.TaCoolerSelfZL + m.TaCoolerSelfa3 == 0

        test_model.TaCoolerSelfZL_Calc = Constraint(rule=TaCoolerSelfZL_Calc)

        def TaCoolerSelfZL_Con1(m):
            return 3 * m.TaCoolerSelfZL ** 2 + 2 * m.TaCoolerSelfa1 * m.TaCoolerSelfZL + m.TaCoolerSelfa2 >= 0

        test_model.TaCoolerSelfZL_Con1 = Constraint(rule=TaCoolerSelfZL_Con1)

        def TaCoolerSelfZL_Con2(m):
            return 6 * m.TaCoolerSelfZL + 2 * m.TaCoolerSelfa1 <= 0

        test_model.TaCoolerSelfZL_Con2 = Constraint(rule=TaCoolerSelfZL_Con2)

        def TaCoolerSelfZV_Calc(m):
            return m.TaCoolerSelfZV ** 3 + m.TaCoolerSelfb1 * m.TaCoolerSelfZV ** 2 + \
                m.TaCoolerSelfb2 * m.TaCoolerSelfZV + m.TaCoolerSelfb3 == 0

        test_model.TaCoolerSelfZV_Calc = Constraint(rule=TaCoolerSelfZV_Calc)

        def TaCoolerSelfZV_Con1(m):
            return 3 * m.TaCoolerSelfZV ** 2 + 2 * m.TaCoolerSelfb1 * m.TaCoolerSelfZV + m.TaCoolerSelfb2 >= 0

        test_model.TaCoolerSelfZV_Con1 = Constraint(rule=TaCoolerSelfZV_Con1)

        def TaCoolerSelfZV_Con2(m):
            return 6 * m.TaCoolerSelfZV + 2 * m.TaCoolerSelfb1 >= 0

        test_model.TaCoolerSelfZV_Con2 = Constraint(rule=TaCoolerSelfZV_Con2)

        def TaCoolerSelfZLBL_Calc(m):
            return exp(m.TaCoolerSelfZLBL) - (m.TaCoolerSelfZL - m.TaCoolerSelfBL) == 0

        test_model.TaCoolerSelfZLBL_Calc = Constraint(rule=TaCoolerSelfZLBL_Calc)

        def TaCoolerSelfZLdivBL_Calc(m):
            return exp(m.TaCoolerSelfZLdivBL) - (m.TaCoolerSelfZL + 2.414 * m.TaCoolerSelfBL) / (
                    m.TaCoolerSelfZL - 0.414 * m.TaCoolerSelfBL) == 0

        test_model.TaCoolerSelfZLdivBL_Calc = Constraint(rule=TaCoolerSelfZLdivBL_Calc)

        def TaCoolerSelfZVBV_Calc(m):
            return exp(m.TaCoolerSelfZVBV) - (m.TaCoolerSelfZV - m.TaCoolerSelfBV) == 0

        test_model.TaCoolerSelfZVBV_Calc = Constraint(rule=TaCoolerSelfZVBV_Calc)

        def TaCoolerSelfZVdivBV_Calc(m):
            return exp(m.TaCoolerSelfZVdivBV) - (m.TaCoolerSelfZV + 2.414 * m.TaCoolerSelfBV) / (
                    m.TaCoolerSelfZV - 0.414 * m.TaCoolerSelfBV) == 0

        test_model.TaCoolerSelfZVdivBV_Calc = Constraint(rule=TaCoolerSelfZVdivBV_Calc)

        def TaCoolerSelfK_Calc(m, comp):
            return m.TaCoolerSelfPhiL[comp] / m.TaCoolerSelfPhiV[comp] - m.TaCoolerSelfVLE_K[comp] == 0

        test_model.TaCoolerSelfK_Calc = Constraint(test_model.Component, rule=TaCoolerSelfK_Calc)

        def TaCoolerSelfVLEEqu(m, comp):
            return m.TaCoolerLiqMFrac[comp] * m.TaCoolerSelfVLE_K[comp] - m.TaCoolerVapMFrac[comp] == 0

        test_model.TaCoolerSelfVLEEqu = Constraint(test_model.Component, rule=TaCoolerSelfVLEEqu)

        # --------TaCooler Summation---------
        def TaCoolerLiqSum(m):
            return sum([m.TaCoolerLiqMFrac[c] for c in m.Component]) == 1

        test_model.TaCoolerLiqSum = Constraint(rule=TaCoolerLiqSum)

        def TaCoolerVapSum(m):
            return sum([m.TaCoolerVapMFrac[c] for c in m.Component]) == 1

        test_model.TaCoolerVapSum = Constraint(rule=TaCoolerVapSum)

        # ----------------------------------
        #           HpaCooler
        # ----------------------------------
        # --------HpaCooler Mass Balance---------
        def HpaCoolerMassBlnc(m, comp):
            return m.HpaCoolerVapMFlow * m.HpaCoolerVapMFrac[comp] + \
                m.HpaCoolerLiqMFlow * m.HpaCoolerLiqMFrac[comp] - m.FeedSplitterOut3MFlow * m.FeedMFrac[comp] == 0

        test_model.HpaCoolerMassBlnc = Constraint(test_model.Component, rule=HpaCoolerMassBlnc)

        def HpaCoolerVF_Spec(m):
            return m.HpaCoolerVF * m.FeedSplitterOut3MFlow - m.HpaCoolerVapMFlow == 0

        test_model.HpaCoolerVF_Spec = Constraint(rule=HpaCoolerVF_Spec)

        # --------HpaCooler Energy Balance---------
        def HpaCoolerEngBlnc(m):
            return ((m.FeedMEtlp - m.HpaCoolerMEtlp) * m.FeedSplitterOut3MFlow - m.HpaCoolerMHeat) * 1.000000e-05 == 0

        test_model.HpaCoolerEngBlnc = Constraint(rule=HpaCoolerEngBlnc)

        # --------HpaCooler Phase Equilibrium---------
        def HpaCoolerSelfSqrtTr_Calc(m, comp):
            return m.HpaCoolerSelfSqrtTr[comp] ** 2 - m.HpaCoolerSelfTr[comp] == 0

        test_model.HpaCoolerSelfSqrtTr_Calc = Constraint(test_model.Component, rule=HpaCoolerSelfSqrtTr_Calc)

        def HpaCoolerSelfSqrt_ai_Calc(m, comp1, comp2):
            return m.HpaCoolerSelfSqrt_ai[comp1, comp2] ** 2 - m.HpaCoolerSelfai[comp1] * m.HpaCoolerSelfai[comp2] == 0

        test_model.HpaCoolerSelfSqrt_ai_Calc = Constraint(test_model.Component, test_model.Component,
                                                          rule=HpaCoolerSelfSqrt_ai_Calc)

        def HpaCoolerSelfZL_Calc(m):
            return m.HpaCoolerSelfZL ** 3 + m.HpaCoolerSelfa1 * m.HpaCoolerSelfZL ** 2 + \
                m.HpaCoolerSelfa2 * m.HpaCoolerSelfZL + m.HpaCoolerSelfa3 == 0

        test_model.HpaCoolerSelfZL_Calc = Constraint(rule=HpaCoolerSelfZL_Calc)

        def HpaCoolerSelfZL_Con1(m):
            return 3 * m.HpaCoolerSelfZL ** 2 + 2 * m.HpaCoolerSelfa1 * m.HpaCoolerSelfZL + m.HpaCoolerSelfa2 >= 0

        test_model.HpaCoolerSelfZL_Con1 = Constraint(rule=HpaCoolerSelfZL_Con1)

        def HpaCoolerSelfZL_Con2(m):
            return 6 * m.HpaCoolerSelfZL + 2 * m.HpaCoolerSelfa1 <= 0

        test_model.HpaCoolerSelfZL_Con2 = Constraint(rule=HpaCoolerSelfZL_Con2)

        def HpaCoolerSelfZV_Calc(m):
            return m.HpaCoolerSelfZV ** 3 + m.HpaCoolerSelfb1 * m.HpaCoolerSelfZV ** 2 + \
                m.HpaCoolerSelfb2 * m.HpaCoolerSelfZV + m.HpaCoolerSelfb3 == 0

        test_model.HpaCoolerSelfZV_Calc = Constraint(rule=HpaCoolerSelfZV_Calc)

        def HpaCoolerSelfZV_Con1(m):
            return 3 * m.HpaCoolerSelfZV ** 2 + 2 * m.HpaCoolerSelfb1 * m.HpaCoolerSelfZV + m.HpaCoolerSelfb2 >= 0

        test_model.HpaCoolerSelfZV_Con1 = Constraint(rule=HpaCoolerSelfZV_Con1)

        def HpaCoolerSelfZV_Con2(m):
            return 6 * m.HpaCoolerSelfZV + 2 * m.HpaCoolerSelfb1 >= 0

        test_model.HpaCoolerSelfZV_Con2 = Constraint(rule=HpaCoolerSelfZV_Con2)

        def HpaCoolerSelfZLBL_Calc(m):
            return exp(m.HpaCoolerSelfZLBL) - (m.HpaCoolerSelfZL - m.HpaCoolerSelfBL) == 0

        test_model.HpaCoolerSelfZLBL_Calc = Constraint(rule=HpaCoolerSelfZLBL_Calc)

        def HpaCoolerSelfZLdivBL_Calc(m):
            return exp(m.HpaCoolerSelfZLdivBL) - (m.HpaCoolerSelfZL + 2.414 * m.HpaCoolerSelfBL) / (
                    m.HpaCoolerSelfZL - 0.414 * m.HpaCoolerSelfBL) == 0

        test_model.HpaCoolerSelfZLdivBL_Calc = Constraint(rule=HpaCoolerSelfZLdivBL_Calc)

        def HpaCoolerSelfZVBV_Calc(m):
            return exp(m.HpaCoolerSelfZVBV) - (m.HpaCoolerSelfZV - m.HpaCoolerSelfBV) == 0

        test_model.HpaCoolerSelfZVBV_Calc = Constraint(rule=HpaCoolerSelfZVBV_Calc)

        def HpaCoolerSelfZVdivBV_Calc(m):
            return exp(m.HpaCoolerSelfZVdivBV) - (m.HpaCoolerSelfZV + 2.414 * m.HpaCoolerSelfBV) / (
                    m.HpaCoolerSelfZV - 0.414 * m.HpaCoolerSelfBV) == 0

        test_model.HpaCoolerSelfZVdivBV_Calc = Constraint(rule=HpaCoolerSelfZVdivBV_Calc)

        def HpaCoolerSelfK_Calc(m, comp):
            return m.HpaCoolerSelfPhiL[comp] / m.HpaCoolerSelfPhiV[comp] - m.HpaCoolerSelfVLE_K[comp] == 0

        test_model.HpaCoolerSelfK_Calc = Constraint(test_model.Component, rule=HpaCoolerSelfK_Calc)

        def HpaCoolerSelfVLEEqu(m, comp):
            return m.HpaCoolerLiqMFrac[comp] * m.HpaCoolerSelfVLE_K[comp] - m.HpaCoolerVapMFrac[comp] == 0

        test_model.HpaCoolerSelfVLEEqu = Constraint(test_model.Component, rule=HpaCoolerSelfVLEEqu)

        # --------HpaCooler Summation---------
        def HpaCoolerLiqSum(m):
            return sum([m.HpaCoolerLiqMFrac[c] for c in m.Component]) == 1

        test_model.HpaCoolerLiqSum = Constraint(rule=HpaCoolerLiqSum)

        def HpaCoolerVapSum(m):
            return sum([m.HpaCoolerVapMFrac[c] for c in m.Component]) == 1

        test_model.HpaCoolerVapSum = Constraint(rule=HpaCoolerVapSum)

        # ----------------------------------
        #           MaCooler
        # ----------------------------------
        # --------MaCooler Enthalpy---------
        def MaCoolerSelfSqrtTr_Calc(m, comp):
            return m.MaCoolerSelfSqrtTr[comp] ** 2 - m.MaCoolerSelfTr[comp] == 0

        test_model.MaCoolerSelfSqrtTr_Calc = Constraint(test_model.Component, rule=MaCoolerSelfSqrtTr_Calc)

        def MaCoolerSelfSqrt_ai_Calc(m, comp1, comp2):
            return m.MaCoolerSelfSqrt_ai[comp1, comp2] ** 2 - m.MaCoolerSelfai[comp1] * m.MaCoolerSelfai[comp2] == 0

        test_model.MaCoolerSelfSqrt_ai_Calc = Constraint(test_model.Component, test_model.Component,
                                                         rule=MaCoolerSelfSqrt_ai_Calc)

        def MaCoolerSelfZV_Calc(m):
            return m.MaCoolerSelfZV ** 3 + m.MaCoolerSelfb1 * m.MaCoolerSelfZV ** 2 + \
                m.MaCoolerSelfb2 * m.MaCoolerSelfZV + m.MaCoolerSelfb3 == 0

        test_model.MaCoolerSelfZV_Calc = Constraint(rule=MaCoolerSelfZV_Calc)

        def MaCoolerSelfZV_Con1(m):
            return 3 * m.MaCoolerSelfZV ** 2 + 2 * m.MaCoolerSelfb1 * m.MaCoolerSelfZV + m.MaCoolerSelfb2 >= 0

        test_model.MaCoolerSelfZV_Con1 = Constraint(rule=MaCoolerSelfZV_Con1)

        def MaCoolerSelfZV_Con2(m):
            return 6 * m.MaCoolerSelfZV + 2 * m.MaCoolerSelfb1 >= 0

        test_model.MaCoolerSelfZV_Con2 = Constraint(rule=MaCoolerSelfZV_Con2)

        def MaCoolerSelfZVdivBV_Calc(m):
            return exp(m.MaCoolerSelfZVdivBV) - (m.MaCoolerSelfZV + 2.414 * m.MaCoolerSelfBV) / (
                    m.MaCoolerSelfZV - 0.414 * m.MaCoolerSelfBV) == 0

        test_model.MaCoolerSelfZVdivBV_Calc = Constraint(rule=MaCoolerSelfZVdivBV_Calc)

        # --------MaCooler Energy Balance---------
        def MaCoolerEngBlnc(m):
            return ((m.FeedMEtlp - m.MaCoolerMEtlp) * m.FeedSplitterOut2MFlow - m.MaCoolerMHeat) * 1.000000e-04 == 0

        test_model.MaCoolerEngBlnc = Constraint(rule=MaCoolerEngBlnc)

        # ----------------------------------
        #           AirMixer
        # ----------------------------------
        # --------AirMixer Mass Balance---------
        def AirMixerOutMassBlnc(m):
            return (m.AirMixerOutMFlow - m.FeedSplitterOut2MFlow - m.FeedSplitterOut1MFlow) * 0.001000 == 0

        test_model.AirMixerOutMassBlnc = Constraint(rule=AirMixerOutMassBlnc)

        def AirMixerOutCompMassBlnc(m, comp):
            return (m.AirMixerOutMFrac[comp] * m.AirMixerOutMFlow - m.FeedMFrac[comp] * m.FeedSplitterOut2MFlow -
                    m.FeedMFrac[comp] * m.FeedSplitterOut1MFlow) * 0.001000 == 0

        test_model.AirMixerOutCompMassBlnc = Constraint(test_model.Component, rule=AirMixerOutCompMassBlnc)

        # --------AirMixer Energy Balance---------
        def AirMixerEnergyBlnc(m):
            return (
                    m.AirMixerOutMEtlp * m.AirMixerOutMFlow - m.MaCoolerMEtlp * m.FeedSplitterOut2MFlow - m.TaCoolerMEtlp * m.FeedSplitterOut1MFlow) * 1.000000e-08 == 0

        test_model.AirMixerEnergyBlnc = Constraint(rule=AirMixerEnergyBlnc)

        # ----------------------------------
        #           HPCZeroReboiled
        # ----------------------------------

        # ----------------------------------
        #           HPC
        # ----------------------------------
        # --------HPC Mass Balance---------
        def HPCMassBlnc(m, tray, comp):
            if tray == 1:
                return (m.HPCZeroReboiledNull * m.HPCZeroReboiledNull + m.HPCLiqLvMCompFlow[tray + 1, comp] -
                        m.HPCLiqLvMCompFlow[tray, comp] - m.HPCVapLvMCompFlow[tray, comp] + m.AirMixerOutMFlow *
                        m.AirMixerOutMFrac[comp]) * 0.001000 == 0
            elif tray == 4:
                return (m.HPCLiqLvMCompFlow[tray + 1, comp] + m.HPCVapLvMCompFlow[tray - 1, comp] - m.HPCLiqLvMCompFlow[
                    tray, comp] - m.HPCVapLvMCompFlow[tray, comp] + m.FeedSplitterOut3MFlow * m.FeedMFrac[
                            comp]) * 0.001000 == 0
            elif tray == 41:
                return (m.HPCCondRefMFlow * m.HPCVapLvMFrac[41, comp] + m.HPCVapLvMCompFlow[tray - 1, comp] -
                        m.HPCLiqLvMCompFlow[tray, comp] - m.HPCVapLvMCompFlow[tray, comp]) * 0.001000 == 0
            else:
                return (m.HPCLiqLvMCompFlow[tray + 1, comp] + m.HPCVapLvMCompFlow[tray - 1, comp] - m.HPCLiqLvMCompFlow[
                    tray, comp] - m.HPCVapLvMCompFlow[tray, comp]) * 0.001000 == 0

        test_model.HPCMassBlnc = Constraint(test_model.HPCTrays, test_model.Component, rule=HPCMassBlnc)

        # --------HPC Energy Balance---------
        def HPCEngBlnc(m, tray):
            if tray == 1:
                return (m.HPCZeroReboiledNull * m.HPCZeroReboiledNull \
                        + m.HPCLiqLvMFlow[tray + 1] * m.HPCLiqLvMEtlp[tray + 1] \
                        - m.HPCLiqLvMFlow[tray] * m.HPCLiqLvMEtlp[tray] \
                        - m.HPCVapLvMFlow[tray] * m.HPCVapLvMEtlp[tray] \
                        + m.AirMixerOutMFlow * m.AirMixerOutMEtlp) * 0.000010 == 0
            elif tray == 4:
                return (m.HPCLiqLvMFlow[tray + 1] * m.HPCLiqLvMEtlp[tray + 1] \
                        + m.HPCVapLvMFlow[tray - 1] * m.HPCVapLvMEtlp[tray - 1] \
                        - m.HPCLiqLvMFlow[tray] * m.HPCLiqLvMEtlp[tray] \
                        - m.HPCVapLvMFlow[tray] * m.HPCVapLvMEtlp[tray] \
                        + m.FeedSplitterOut3MFlow * m.HpaCoolerMEtlp) * 0.000010 == 0
            elif tray == 41:
                return (m.HPCCondRefMFlow * m.HPCCondOutMEtlp \
                        + m.HPCVapLvMFlow[tray - 1] * m.HPCVapLvMEtlp[tray - 1] \
                        - m.HPCLiqLvMFlow[tray] * m.HPCLiqLvMEtlp[tray] \
                        - m.HPCVapLvMFlow[tray] * m.HPCVapLvMEtlp[tray]) * 0.000010 == 0
            else:
                return (m.HPCLiqLvMFlow[tray + 1] * m.HPCLiqLvMEtlp[tray + 1] \
                        + m.HPCVapLvMFlow[tray - 1] * m.HPCVapLvMEtlp[tray - 1] \
                        - m.HPCLiqLvMFlow[tray] * m.HPCLiqLvMEtlp[tray] \
                        - m.HPCVapLvMFlow[tray] * m.HPCVapLvMEtlp[tray]) * 0.000010 == 0

        test_model.HPCEngBlnc = Constraint(test_model.HPCTrays, rule=HPCEngBlnc)

        # --------HPC Phase Equilibrium & System Parts---------
        def HPCAllTraysSqrtTr_Calc(m, tray, comp):
            return m.HPCAllTraysSqrtTr[tray, comp] ** 2 - m.HPCAllTraysTr[tray, comp] == 0

        test_model.HPCAllTraysSqrtTr_Calc = Constraint(test_model.HPCTrays, test_model.Component,
                                                       rule=HPCAllTraysSqrtTr_Calc)

        def HPCAllTraysSqrt_ai_Calc(m, tray, comp1, comp2):
            return m.HPCAllTraysSqrt_ai[tray, comp1, comp2] ** 2 - m.HPCAllTraysai[tray, comp1] * m.HPCAllTraysai[
                tray, comp2] == 0

        test_model.HPCAllTraysSqrt_ai_Calc = Constraint(test_model.HPCTrays, test_model.Component, test_model.Component,
                                                        rule=HPCAllTraysSqrt_ai_Calc)

        def HPCAllTraysZL_Calc(m, tray):
            return m.HPCAllTraysZL[tray] ** 3 + m.HPCAllTraysa1[tray] * m.HPCAllTraysZL[tray] ** 2 + \
                m.HPCAllTraysa2[tray] * m.HPCAllTraysZL[tray] + m.HPCAllTraysa3[tray] == 0

        test_model.HPCAllTraysZL_Calc = Constraint(test_model.HPCTrays, rule=HPCAllTraysZL_Calc)

        def HPCAllTraysZL_Con1(m, tray):
            return 3 * m.HPCAllTraysZL[tray] ** 2 + 2 * m.HPCAllTraysa1[tray] * m.HPCAllTraysZL[tray] + m.HPCAllTraysa2[
                tray] >= 0

        test_model.HPCAllTraysZL_Con1 = Constraint(test_model.HPCTrays, rule=HPCAllTraysZL_Con1)

        def HPCAllTraysZL_Con2(m, tray):
            return 6 * m.HPCAllTraysZL[tray] + 2 * m.HPCAllTraysa1[tray] <= 0

        test_model.HPCAllTraysZL_Con2 = Constraint(test_model.HPCTrays, rule=HPCAllTraysZL_Con2)

        def HPCAllTraysZV_Calc(m, tray):
            return m.HPCAllTraysZV[tray] ** 3 + m.HPCAllTraysb1[tray] * m.HPCAllTraysZV[tray] ** 2 + \
                m.HPCAllTraysb2[tray] * m.HPCAllTraysZV[tray] + m.HPCAllTraysb3[tray] == 0

        test_model.HPCAllTraysZV_Calc = Constraint(test_model.HPCTrays, rule=HPCAllTraysZV_Calc)

        def HPCAllTraysZV_Con1(m, tray):
            return 3 * m.HPCAllTraysZV[tray] ** 2 + 2 * m.HPCAllTraysb1[tray] * m.HPCAllTraysZV[tray] + m.HPCAllTraysb2[
                tray] >= 0

        test_model.HPCAllTraysZV_Con1 = Constraint(test_model.HPCTrays, rule=HPCAllTraysZV_Con1)

        def HPCAllTraysZV_Con2(m, tray):
            return 6 * m.HPCAllTraysZV[tray] + 2 * m.HPCAllTraysb1[tray] >= 0

        test_model.HPCAllTraysZV_Con2 = Constraint(test_model.HPCTrays, rule=HPCAllTraysZV_Con2)

        def HPCAllTraysZLBL_Calc(m, tray):
            return exp(m.HPCAllTraysZLBL[tray]) - (m.HPCAllTraysZL[tray] - m.HPCAllTraysBL[tray]) == 0

        test_model.HPCAllTraysZLBL_Calc = Constraint(test_model.HPCTrays, rule=HPCAllTraysZLBL_Calc)

        def HPCAllTraysZLdivBL_Calc(m, tray):
            return exp(m.HPCAllTraysZLdivBL[tray]) - (m.HPCAllTraysZL[tray] + 2.414 * m.HPCAllTraysBL[tray]) / (
                    m.HPCAllTraysZL[tray] - 0.414 * m.HPCAllTraysBL[tray]) == 0

        test_model.HPCAllTraysZLdivBL_Calc = Constraint(test_model.HPCTrays, rule=HPCAllTraysZLdivBL_Calc)

        def HPCAllTraysZVBV_Calc(m, tray):
            return exp(m.HPCAllTraysZVBV[tray]) - (m.HPCAllTraysZV[tray] - m.HPCAllTraysBV[tray]) == 0

        test_model.HPCAllTraysZVBV_Calc = Constraint(test_model.HPCTrays, rule=HPCAllTraysZVBV_Calc)

        def HPCAllTraysZVdivBV_Calc(m, tray):
            return exp(m.HPCAllTraysZVdivBV[tray]) - (m.HPCAllTraysZV[tray] + 2.414 * m.HPCAllTraysBV[tray]) / (
                    m.HPCAllTraysZV[tray] - 0.414 * m.HPCAllTraysBV[tray]) == 0

        test_model.HPCAllTraysZVdivBV_Calc = Constraint(test_model.HPCTrays, rule=HPCAllTraysZVdivBV_Calc)

        def HPCAllTraysK_Calc(m, tray, comp):
            return m.HPCAllTraysPhiL[tray, comp] / m.HPCAllTraysPhiV[tray, comp] - m.HPCAllTraysVLE_K[tray, comp] == 0

        test_model.HPCAllTraysK_Calc = Constraint(test_model.HPCTrays, test_model.Component, rule=HPCAllTraysK_Calc)

        def HPCAllTraysVLEEqu(m, tray, comp):
            return m.HPCLiqLvMFrac[tray, comp] * m.HPCAllTraysVLE_K[tray, comp] - m.HPCVapLvMFrac[tray, comp] == 0

        test_model.HPCAllTraysVLEEqu = Constraint(test_model.HPCTrays, test_model.Component, rule=HPCAllTraysVLEEqu)

        # --------HPC Summation---------
        def HPCLiqSum(m, tray):
            return sum([m.HPCLiqLvMFrac[tray, c] for c in m.Component]) == 1

        test_model.HPCLiqSum = Constraint(test_model.HPCTrays, rule=HPCLiqSum)

        def HPCVapSum(m, tray):
            return sum([m.HPCVapLvMFrac[tray, c] for c in m.Component]) == 1

        test_model.HPCVapSum = Constraint(test_model.HPCTrays, rule=HPCVapSum)

        # --------HPC Pressure Profile---------
        def HPCPProf(m, tray):
            return ((m.HPCTopPressure - m.HPCBtmPressure) / 40 * (tray - 1) + m.HPCBtmPressure - m.HPCTrayPressure[
                tray]) * 0.010000 == 0

        test_model.HPCPProf = Constraint(test_model.HPCTrays, rule=HPCPProf)

        # ----------------------------------
        #           HPCSump
        # ----------------------------------
        # --------HPCSump Mass Balance---------
        def HPCSumpMassBlnc(m):
            return (m.HPCLiqLvMFlow[1] - m.HPCSumpOutMFlow) * 0.001000 == 0

        test_model.HPCSumpMassBlnc = Constraint(rule=HPCSumpMassBlnc)

        def HPCSumpHldpMFracSpec(m, comp):
            return m.HPCSumpHldpMFrac[comp] - m.HPCLiqLvMFrac[1, comp] == 0

        test_model.HPCSumpHldpMFracSpec = Constraint(test_model.Component, rule=HPCSumpHldpMFracSpec)

        # ----------------------------------
        #           HPCCond
        # ----------------------------------
        # --------HPCCond Mass Balance---------
        def HPCCondMassBlnc(m):
            return (m.HPCCondRefMFlow + m.HPCCondPrdtMFlow - m.HPCVapLvMFlow[41]) * 0.001000 == 0

        test_model.HPCCondMassBlnc = Constraint(rule=HPCCondMassBlnc)

        # --------HPCCond Energy Balance---------
        def HPCCondEngBlnc(m):
            return (m.HPCVapLvMFlow[41] * (m.HPCVapLvMEtlp[41] - m.HPCCondOutMEtlp) - m.HPCCondMHeatOut) * 0.000010 == 0

        test_model.HPCCondEngBlnc = Constraint(rule=HPCCondEngBlnc)

        # --------HPCCond Enthalpy---------
        def HPCCondOutletSqrtTr_Calc(m, comp):
            return m.HPCCondOutletSqrtTr[comp] ** 2 - m.HPCCondOutletTr[comp] == 0

        test_model.HPCCondOutletSqrtTr_Calc = Constraint(test_model.Component, rule=HPCCondOutletSqrtTr_Calc)

        def HPCCondOutletSqrt_ai_Calc(m, comp1, comp2):
            return m.HPCCondOutletSqrt_ai[comp1, comp2] ** 2 - m.HPCCondOutletai[comp1] * m.HPCCondOutletai[comp2] == 0

        test_model.HPCCondOutletSqrt_ai_Calc = Constraint(test_model.Component, test_model.Component,
                                                          rule=HPCCondOutletSqrt_ai_Calc)

        def HPCCondOutletZL_Calc(m):
            return m.HPCCondOutletZL ** 3 + m.HPCCondOutleta1 * m.HPCCondOutletZL ** 2 + \
                m.HPCCondOutleta2 * m.HPCCondOutletZL + m.HPCCondOutleta3 == 0

        test_model.HPCCondOutletZL_Calc = Constraint(rule=HPCCondOutletZL_Calc)

        def HPCCondOutletZL_Con1(m):
            return 3 * m.HPCCondOutletZL ** 2 + 2 * m.HPCCondOutleta1 * m.HPCCondOutletZL + m.HPCCondOutleta2 >= 0

        test_model.HPCCondOutletZL_Con1 = Constraint(rule=HPCCondOutletZL_Con1)

        def HPCCondOutletZL_Con2(m):
            return 6 * m.HPCCondOutletZL + 2 * m.HPCCondOutleta1 <= 0

        test_model.HPCCondOutletZL_Con2 = Constraint(rule=HPCCondOutletZL_Con2)

        def HPCCondOutletZLdivBL_Calc(m):
            return exp(m.HPCCondOutletZLdivBL) - (m.HPCCondOutletZL + 2.414 * m.HPCCondOutletBL) / (
                    m.HPCCondOutletZL - 0.414 * m.HPCCondOutletBL) == 0

        test_model.HPCCondOutletZLdivBL_Calc = Constraint(rule=HPCCondOutletZLdivBL_Calc)

        # ----------------------------------
        #           CoolLin
        # ----------------------------------
        # --------CoolLin DeltaT or DeltaP---------
        def CoolLinDeltaPSpec(m):
            return (m.HPCCondPressure - m.CoolLinDeltaP - m.CoolLinPressure) * 1e-2 == 0

        test_model.CoolLinDeltaPSpec = Constraint(rule=CoolLinDeltaPSpec)

        # --------CoolLin Enthalpy---------
        def CoolLinSelfSqrtTr_Calc(m, comp):
            return m.CoolLinSelfSqrtTr[comp] ** 2 - m.CoolLinSelfTr[comp] == 0

        test_model.CoolLinSelfSqrtTr_Calc = Constraint(test_model.Component, rule=CoolLinSelfSqrtTr_Calc)

        def CoolLinSelfSqrt_ai_Calc(m, comp1, comp2):
            return m.CoolLinSelfSqrt_ai[comp1, comp2] ** 2 - m.CoolLinSelfai[comp1] * m.CoolLinSelfai[comp2] == 0

        test_model.CoolLinSelfSqrt_ai_Calc = Constraint(test_model.Component, test_model.Component,
                                                        rule=CoolLinSelfSqrt_ai_Calc)

        def CoolLinSelfZL_Calc(m):
            return m.CoolLinSelfZL ** 3 + m.CoolLinSelfa1 * m.CoolLinSelfZL ** 2 + \
                m.CoolLinSelfa2 * m.CoolLinSelfZL + m.CoolLinSelfa3 == 0

        test_model.CoolLinSelfZL_Calc = Constraint(rule=CoolLinSelfZL_Calc)

        def CoolLinSelfZL_Con1(m):
            return 3 * m.CoolLinSelfZL ** 2 + 2 * m.CoolLinSelfa1 * m.CoolLinSelfZL + m.CoolLinSelfa2 >= 0

        test_model.CoolLinSelfZL_Con1 = Constraint(rule=CoolLinSelfZL_Con1)

        def CoolLinSelfZL_Con2(m):
            return 6 * m.CoolLinSelfZL + 2 * m.CoolLinSelfa1 <= 0

        test_model.CoolLinSelfZL_Con2 = Constraint(rule=CoolLinSelfZL_Con2)

        def CoolLinSelfZLdivBL_Calc(m):
            return exp(m.CoolLinSelfZLdivBL) - (m.CoolLinSelfZL + 2.414 * m.CoolLinSelfBL) / (
                    m.CoolLinSelfZL - 0.414 * m.CoolLinSelfBL) == 0

        test_model.CoolLinSelfZLdivBL_Calc = Constraint(rule=CoolLinSelfZLdivBL_Calc)

        # --------CoolLin Energy Balance---------
        def CoolLinEngBlnc(m):
            return ((m.HPCCondOutMEtlp - m.CoolLinMEtlp) * m.HPCCondPrdtMFlow - m.CoolLinMHeat) * 1.000000e-04 == 0

        test_model.CoolLinEngBlnc = Constraint(rule=CoolLinEngBlnc)

        # ----------------------------------
        #           QTrtl
        # ----------------------------------

        # ----------------------------------
        #           Throttle
        # ----------------------------------
        # --------Throttle Mass Balance---------
        def ThrottleMassBlnc(m, comp):
            return m.ThrottleVapMFlow * m.ThrottleVapMFrac[comp] + \
                m.ThrottleLiqMFlow * m.ThrottleLiqMFrac[comp] - m.HPCSumpOutMFlow * m.HPCLiqLvMFrac[1, comp] == 0

        test_model.ThrottleMassBlnc = Constraint(test_model.Component, rule=ThrottleMassBlnc)

        # --------Throttle Energy Balance---------
        def ThrottleEngBlnc(m):
            return ((m.HPCLiqLvMEtlp[1] - m.ThrottleMEtlp) * m.HPCSumpOutMFlow + m.QTrtlHeat) * 1.000000e-04 == 0

        test_model.ThrottleEngBlnc = Constraint(rule=ThrottleEngBlnc)

        # --------Throttle Phase Equilibrium---------
        def ThrottleSelfSqrtTr_Calc(m, comp):
            return m.ThrottleSelfSqrtTr[comp] ** 2 - m.ThrottleSelfTr[comp] == 0

        test_model.ThrottleSelfSqrtTr_Calc = Constraint(test_model.Component, rule=ThrottleSelfSqrtTr_Calc)

        def ThrottleSelfSqrt_ai_Calc(m, comp1, comp2):
            return m.ThrottleSelfSqrt_ai[comp1, comp2] ** 2 - m.ThrottleSelfai[comp1] * m.ThrottleSelfai[comp2] == 0

        test_model.ThrottleSelfSqrt_ai_Calc = Constraint(test_model.Component, test_model.Component,
                                                         rule=ThrottleSelfSqrt_ai_Calc)

        def ThrottleSelfZL_Calc(m):
            return m.ThrottleSelfZL ** 3 + m.ThrottleSelfa1 * m.ThrottleSelfZL ** 2 + \
                m.ThrottleSelfa2 * m.ThrottleSelfZL + m.ThrottleSelfa3 == 0

        test_model.ThrottleSelfZL_Calc = Constraint(rule=ThrottleSelfZL_Calc)

        def ThrottleSelfZL_Con1(m):
            return 3 * m.ThrottleSelfZL ** 2 + 2 * m.ThrottleSelfa1 * m.ThrottleSelfZL + m.ThrottleSelfa2 >= 0

        test_model.ThrottleSelfZL_Con1 = Constraint(rule=ThrottleSelfZL_Con1)

        def ThrottleSelfZL_Con2(m):
            return 6 * m.ThrottleSelfZL + 2 * m.ThrottleSelfa1 <= 0

        test_model.ThrottleSelfZL_Con2 = Constraint(rule=ThrottleSelfZL_Con2)

        def ThrottleSelfZV_Calc(m):
            return m.ThrottleSelfZV ** 3 + m.ThrottleSelfb1 * m.ThrottleSelfZV ** 2 + \
                m.ThrottleSelfb2 * m.ThrottleSelfZV + m.ThrottleSelfb3 == 0

        test_model.ThrottleSelfZV_Calc = Constraint(rule=ThrottleSelfZV_Calc)

        def ThrottleSelfZV_Con1(m):
            return 3 * m.ThrottleSelfZV ** 2 + 2 * m.ThrottleSelfb1 * m.ThrottleSelfZV + m.ThrottleSelfb2 >= 0

        test_model.ThrottleSelfZV_Con1 = Constraint(rule=ThrottleSelfZV_Con1)

        def ThrottleSelfZV_Con2(m):
            return 6 * m.ThrottleSelfZV + 2 * m.ThrottleSelfb1 >= 0

        test_model.ThrottleSelfZV_Con2 = Constraint(rule=ThrottleSelfZV_Con2)

        def ThrottleSelfZLBL_Calc(m):
            return exp(m.ThrottleSelfZLBL) - (m.ThrottleSelfZL - m.ThrottleSelfBL) == 0

        test_model.ThrottleSelfZLBL_Calc = Constraint(rule=ThrottleSelfZLBL_Calc)

        def ThrottleSelfZLdivBL_Calc(m):
            return exp(m.ThrottleSelfZLdivBL) - (m.ThrottleSelfZL + 2.414 * m.ThrottleSelfBL) / (
                    m.ThrottleSelfZL - 0.414 * m.ThrottleSelfBL) == 0

        test_model.ThrottleSelfZLdivBL_Calc = Constraint(rule=ThrottleSelfZLdivBL_Calc)

        def ThrottleSelfZVBV_Calc(m):
            return exp(m.ThrottleSelfZVBV) - (m.ThrottleSelfZV - m.ThrottleSelfBV) == 0

        test_model.ThrottleSelfZVBV_Calc = Constraint(rule=ThrottleSelfZVBV_Calc)

        def ThrottleSelfZVdivBV_Calc(m):
            return exp(m.ThrottleSelfZVdivBV) - (m.ThrottleSelfZV + 2.414 * m.ThrottleSelfBV) / (
                    m.ThrottleSelfZV - 0.414 * m.ThrottleSelfBV) == 0

        test_model.ThrottleSelfZVdivBV_Calc = Constraint(rule=ThrottleSelfZVdivBV_Calc)

        def ThrottleSelfK_Calc(m, comp):
            return m.ThrottleSelfPhiL[comp] / m.ThrottleSelfPhiV[comp] - m.ThrottleSelfVLE_K[comp] == 0

        test_model.ThrottleSelfK_Calc = Constraint(test_model.Component, rule=ThrottleSelfK_Calc)

        def ThrottleSelfVLEEqu(m, comp):
            return m.ThrottleLiqMFrac[comp] * m.ThrottleSelfVLE_K[comp] - m.ThrottleVapMFrac[comp] == 0

        test_model.ThrottleSelfVLEEqu = Constraint(test_model.Component, rule=ThrottleSelfVLEEqu)

        # --------Throttle Summation---------
        def ThrottleLiqSum(m):
            return sum([m.ThrottleLiqMFrac[c] for c in m.Component]) == 1

        test_model.ThrottleLiqSum = Constraint(rule=ThrottleLiqSum)

        def ThrottleVapSum(m):
            return sum([m.ThrottleVapMFrac[c] for c in m.Component]) == 1

        test_model.ThrottleVapSum = Constraint(rule=ThrottleVapSum)

        # ----------------------------------
        #           HeatLPC
        # ----------------------------------
        # --------HeatLPC DeltaT or DeltaP---------
        def HeatLPCDeltaPSpec(m):
            return (m.ThrottlePressure - m.HeatLPCDeltaP - m.HeatLPCPressure) * 1e-2 == 0

        test_model.HeatLPCDeltaPSpec = Constraint(rule=HeatLPCDeltaPSpec)

        # --------HeatLPC Mass Balance---------
        def HeatLPCMassBlnc(m, comp):
            return m.HeatLPCVapMFlow * m.HeatLPCVapMFrac[comp] + \
                m.HeatLPCLiqMFlow * m.HeatLPCLiqMFrac[comp] - m.HPCSumpOutMFlow * m.HPCLiqLvMFrac[1, comp] == 0

        test_model.HeatLPCMassBlnc = Constraint(test_model.Component, rule=HeatLPCMassBlnc)

        # --------HeatLPC Energy Balance---------
        def HeatLPCEngBlnc(m):
            return ((m.ThrottleMEtlp - m.HeatLPCMEtlp) * m.HPCSumpOutMFlow + m.ASCCondMHeat) * 1.000000e-04 == 0

        test_model.HeatLPCEngBlnc = Constraint(rule=HeatLPCEngBlnc)

        # --------HeatLPC Phase Equilibrium---------
        def HeatLPCSelfSqrtTr_Calc(m, comp):
            return m.HeatLPCSelfSqrtTr[comp] ** 2 - m.HeatLPCSelfTr[comp] == 0

        test_model.HeatLPCSelfSqrtTr_Calc = Constraint(test_model.Component, rule=HeatLPCSelfSqrtTr_Calc)

        def HeatLPCSelfSqrt_ai_Calc(m, comp1, comp2):
            return m.HeatLPCSelfSqrt_ai[comp1, comp2] ** 2 - m.HeatLPCSelfai[comp1] * m.HeatLPCSelfai[comp2] == 0

        test_model.HeatLPCSelfSqrt_ai_Calc = Constraint(test_model.Component, test_model.Component,
                                                        rule=HeatLPCSelfSqrt_ai_Calc)

        def HeatLPCSelfZL_Calc(m):
            return m.HeatLPCSelfZL ** 3 + m.HeatLPCSelfa1 * m.HeatLPCSelfZL ** 2 + \
                m.HeatLPCSelfa2 * m.HeatLPCSelfZL + m.HeatLPCSelfa3 == 0

        test_model.HeatLPCSelfZL_Calc = Constraint(rule=HeatLPCSelfZL_Calc)

        def HeatLPCSelfZL_Con1(m):
            return 3 * m.HeatLPCSelfZL ** 2 + 2 * m.HeatLPCSelfa1 * m.HeatLPCSelfZL + m.HeatLPCSelfa2 >= 0

        test_model.HeatLPCSelfZL_Con1 = Constraint(rule=HeatLPCSelfZL_Con1)

        def HeatLPCSelfZL_Con2(m):
            return 6 * m.HeatLPCSelfZL + 2 * m.HeatLPCSelfa1 <= 0

        test_model.HeatLPCSelfZL_Con2 = Constraint(rule=HeatLPCSelfZL_Con2)

        def HeatLPCSelfZV_Calc(m):
            return m.HeatLPCSelfZV ** 3 + m.HeatLPCSelfb1 * m.HeatLPCSelfZV ** 2 + \
                m.HeatLPCSelfb2 * m.HeatLPCSelfZV + m.HeatLPCSelfb3 == 0

        test_model.HeatLPCSelfZV_Calc = Constraint(rule=HeatLPCSelfZV_Calc)

        def HeatLPCSelfZV_Con1(m):
            return 3 * m.HeatLPCSelfZV ** 2 + 2 * m.HeatLPCSelfb1 * m.HeatLPCSelfZV + m.HeatLPCSelfb2 >= 0

        test_model.HeatLPCSelfZV_Con1 = Constraint(rule=HeatLPCSelfZV_Con1)

        def HeatLPCSelfZV_Con2(m):
            return 6 * m.HeatLPCSelfZV + 2 * m.HeatLPCSelfb1 >= 0

        test_model.HeatLPCSelfZV_Con2 = Constraint(rule=HeatLPCSelfZV_Con2)

        def HeatLPCSelfZLBL_Calc(m):
            return exp(m.HeatLPCSelfZLBL) - (m.HeatLPCSelfZL - m.HeatLPCSelfBL) == 0

        test_model.HeatLPCSelfZLBL_Calc = Constraint(rule=HeatLPCSelfZLBL_Calc)

        def HeatLPCSelfZLdivBL_Calc(m):
            return exp(m.HeatLPCSelfZLdivBL) - (m.HeatLPCSelfZL + 2.414 * m.HeatLPCSelfBL) / (
                    m.HeatLPCSelfZL - 0.414 * m.HeatLPCSelfBL) == 0

        test_model.HeatLPCSelfZLdivBL_Calc = Constraint(rule=HeatLPCSelfZLdivBL_Calc)

        def HeatLPCSelfZVBV_Calc(m):
            return exp(m.HeatLPCSelfZVBV) - (m.HeatLPCSelfZV - m.HeatLPCSelfBV) == 0

        test_model.HeatLPCSelfZVBV_Calc = Constraint(rule=HeatLPCSelfZVBV_Calc)

        def HeatLPCSelfZVdivBV_Calc(m):
            return exp(m.HeatLPCSelfZVdivBV) - (m.HeatLPCSelfZV + 2.414 * m.HeatLPCSelfBV) / (
                    m.HeatLPCSelfZV - 0.414 * m.HeatLPCSelfBV) == 0

        test_model.HeatLPCSelfZVdivBV_Calc = Constraint(rule=HeatLPCSelfZVdivBV_Calc)

        def HeatLPCSelfK_Calc(m, comp):
            return m.HeatLPCSelfPhiL[comp] / m.HeatLPCSelfPhiV[comp] - m.HeatLPCSelfVLE_K[comp] == 0

        test_model.HeatLPCSelfK_Calc = Constraint(test_model.Component, rule=HeatLPCSelfK_Calc)

        def HeatLPCSelfVLEEqu(m, comp):
            return m.HeatLPCLiqMFrac[comp] * m.HeatLPCSelfVLE_K[comp] - m.HeatLPCVapMFrac[comp] == 0

        test_model.HeatLPCSelfVLEEqu = Constraint(test_model.Component, rule=HeatLPCSelfVLEEqu)

        # --------HeatLPC Summation---------
        def HeatLPCLiqSum(m):
            return sum([m.HeatLPCLiqMFrac[c] for c in m.Component]) == 1

        test_model.HeatLPCLiqSum = Constraint(rule=HeatLPCLiqSum)

        def HeatLPCVapSum(m):
            return sum([m.HeatLPCVapMFrac[c] for c in m.Component]) == 1

        test_model.HeatLPCVapSum = Constraint(rule=HeatLPCVapSum)

        # ----------------------------------
        #           LPCZeroReflux
        # ----------------------------------

        # ----------------------------------
        #           LPC
        # ----------------------------------
        # --------LPC Mass Balance---------
        def LPCMassBlnc(m, tray, comp):
            if tray == 1:
                return (m.LPCReboilerVapLvMFlow * m.LPCReboilerVapLvMFrac[comp] + m.LPCLiqLvMCompFlow[tray + 1, comp] -
                        m.LPCLiqLvMCompFlow[tray, comp] - m.LPCVapLvMCompFlow[tray, comp]) * 0.001000 == 0
            elif tray == 15:
                return (m.LPCLiqLvMCompFlow[tray + 1, comp] + m.LPCVapLvMCompFlow[tray - 1, comp] - m.LPCLiqLvMCompFlow[
                    tray, comp] - m.LPCVapLvMCompFlow[tray, comp] + m.ASCSumpOutMFlow * m.ASCLiqLvMFrac[
                            1, comp] - m.LPCExt15MFlow * m.LPCExt15MFrac[comp]) * 0.001000 == 0
            elif tray == 29:
                return (m.LPCLiqLvMCompFlow[tray + 1, comp] + m.LPCVapLvMCompFlow[tray - 1, comp] - m.LPCLiqLvMCompFlow[
                    tray, comp] - m.LPCVapLvMCompFlow[tray, comp] + m.HPCSumpOutMFlow * m.HPCLiqLvMFrac[
                            1, comp]) * 0.001000 == 0
            elif tray == 46:
                return (m.LPCLiqLvMCompFlow[tray + 1, comp] + m.LPCVapLvMCompFlow[tray - 1, comp] - m.LPCLiqLvMCompFlow[
                    tray, comp] - m.LPCVapLvMCompFlow[tray, comp] - m.LPCExt46MFlow * m.LPCExt46MFrac[
                            comp]) * 0.001000 == 0
            elif tray == 52:
                return (m.LPCZeroRefluxNull * m.LPCZeroRefluxNull + m.LPCVapLvMCompFlow[tray - 1, comp] -
                        m.LPCLiqLvMCompFlow[tray, comp] - m.LPCVapLvMCompFlow[tray, comp] + m.HPCCondPrdtMFlow *
                        m.HPCVapLvMFrac[41, comp]) * 0.001000 == 0
            else:
                return (m.LPCLiqLvMCompFlow[tray + 1, comp] + m.LPCVapLvMCompFlow[tray - 1, comp] - m.LPCLiqLvMCompFlow[
                    tray, comp] - m.LPCVapLvMCompFlow[tray, comp]) * 0.001000 == 0

        test_model.LPCMassBlnc = Constraint(test_model.LPCTrays, test_model.Component, rule=LPCMassBlnc)

        # --------LPC Energy Balance---------
        def LPCEngBlnc(m, tray):
            if tray == 1:
                return (m.LPCReboilerVapLvMFlow * m.LPCReboilerVapLvMEtlp \
                        + m.LPCLiqLvMFlow[tray + 1] * m.LPCLiqLvMEtlp[tray + 1] \
                        - m.LPCLiqLvMFlow[tray] * m.LPCLiqLvMEtlp[tray] \
                        - m.LPCVapLvMFlow[tray] * m.LPCVapLvMEtlp[tray]) * 0.000010 == 0
            elif tray == 15:
                return (m.LPCLiqLvMFlow[tray + 1] * m.LPCLiqLvMEtlp[tray + 1] \
                        + m.LPCVapLvMFlow[tray - 1] * m.LPCVapLvMEtlp[tray - 1] \
                        - m.LPCLiqLvMFlow[tray] * m.LPCLiqLvMEtlp[tray] \
                        - m.LPCVapLvMFlow[tray] * m.LPCVapLvMEtlp[tray] \
                        + m.ASCSumpOutMFlow * m.ASCLiqLvMEtlp[1] \
                        - m.LPCExt15MFlow * m.LPCExt15MEtlp) * 0.000010 == 0
            elif tray == 29:
                return (m.LPCLiqLvMFlow[tray + 1] * m.LPCLiqLvMEtlp[tray + 1] \
                        + m.LPCVapLvMFlow[tray - 1] * m.LPCVapLvMEtlp[tray - 1] \
                        - m.LPCLiqLvMFlow[tray] * m.LPCLiqLvMEtlp[tray] \
                        - m.LPCVapLvMFlow[tray] * m.LPCVapLvMEtlp[tray] \
                        + m.HPCSumpOutMFlow * m.HeatLPCMEtlp) * 0.000010 == 0
            elif tray == 46:
                return (m.LPCLiqLvMFlow[tray + 1] * m.LPCLiqLvMEtlp[tray + 1] \
                        + m.LPCVapLvMFlow[tray - 1] * m.LPCVapLvMEtlp[tray - 1] \
                        - m.LPCLiqLvMFlow[tray] * m.LPCLiqLvMEtlp[tray] \
                        - m.LPCVapLvMFlow[tray] * m.LPCVapLvMEtlp[tray] \
                        - m.LPCExt46MFlow * m.LPCExt46MEtlp) * 0.000010 == 0
            elif tray == 52:
                return (m.LPCZeroRefluxNull * m.LPCZeroRefluxNull \
                        + m.LPCVapLvMFlow[tray - 1] * m.LPCVapLvMEtlp[tray - 1] \
                        - m.LPCLiqLvMFlow[tray] * m.LPCLiqLvMEtlp[tray] \
                        - m.LPCVapLvMFlow[tray] * m.LPCVapLvMEtlp[tray] \
                        + m.HPCCondPrdtMFlow * m.CoolLinMEtlp) * 0.000010 == 0
            else:
                return (m.LPCLiqLvMFlow[tray + 1] * m.LPCLiqLvMEtlp[tray + 1] \
                        + m.LPCVapLvMFlow[tray - 1] * m.LPCVapLvMEtlp[tray - 1] \
                        - m.LPCLiqLvMFlow[tray] * m.LPCLiqLvMEtlp[tray] \
                        - m.LPCVapLvMFlow[tray] * m.LPCVapLvMEtlp[tray]) * 0.000010 == 0

        test_model.LPCEngBlnc = Constraint(test_model.LPCTrays, rule=LPCEngBlnc)

        # --------LPC Phase Equilibrium & System Parts---------
        def LPCAllTraysSqrtTr_Calc(m, tray, comp):
            return m.LPCAllTraysSqrtTr[tray, comp] ** 2 - m.LPCAllTraysTr[tray, comp] == 0

        test_model.LPCAllTraysSqrtTr_Calc = Constraint(test_model.LPCTrays, test_model.Component,
                                                       rule=LPCAllTraysSqrtTr_Calc)

        def LPCAllTraysSqrt_ai_Calc(m, tray, comp1, comp2):
            return m.LPCAllTraysSqrt_ai[tray, comp1, comp2] ** 2 - m.LPCAllTraysai[tray, comp1] * m.LPCAllTraysai[
                tray, comp2] == 0

        test_model.LPCAllTraysSqrt_ai_Calc = Constraint(test_model.LPCTrays, test_model.Component, test_model.Component,
                                                        rule=LPCAllTraysSqrt_ai_Calc)

        def LPCAllTraysZL_Calc(m, tray):
            return m.LPCAllTraysZL[tray] ** 3 + m.LPCAllTraysa1[tray] * m.LPCAllTraysZL[tray] ** 2 + \
                m.LPCAllTraysa2[tray] * m.LPCAllTraysZL[tray] + m.LPCAllTraysa3[tray] == 0

        test_model.LPCAllTraysZL_Calc = Constraint(test_model.LPCTrays, rule=LPCAllTraysZL_Calc)

        def LPCAllTraysZL_Con1(m, tray):
            return 3 * m.LPCAllTraysZL[tray] ** 2 + 2 * m.LPCAllTraysa1[tray] * m.LPCAllTraysZL[tray] + m.LPCAllTraysa2[
                tray] >= 0

        test_model.LPCAllTraysZL_Con1 = Constraint(test_model.LPCTrays, rule=LPCAllTraysZL_Con1)

        def LPCAllTraysZL_Con2(m, tray):
            return 6 * m.LPCAllTraysZL[tray] + 2 * m.LPCAllTraysa1[tray] <= 0

        test_model.LPCAllTraysZL_Con2 = Constraint(test_model.LPCTrays, rule=LPCAllTraysZL_Con2)

        def LPCAllTraysZV_Calc(m, tray):
            return m.LPCAllTraysZV[tray] ** 3 + m.LPCAllTraysb1[tray] * m.LPCAllTraysZV[tray] ** 2 + \
                m.LPCAllTraysb2[tray] * m.LPCAllTraysZV[tray] + m.LPCAllTraysb3[tray] == 0

        test_model.LPCAllTraysZV_Calc = Constraint(test_model.LPCTrays, rule=LPCAllTraysZV_Calc)

        def LPCAllTraysZV_Con1(m, tray):
            return 3 * m.LPCAllTraysZV[tray] ** 2 + 2 * m.LPCAllTraysb1[tray] * m.LPCAllTraysZV[tray] + m.LPCAllTraysb2[
                tray] >= 0

        test_model.LPCAllTraysZV_Con1 = Constraint(test_model.LPCTrays, rule=LPCAllTraysZV_Con1)

        def LPCAllTraysZV_Con2(m, tray):
            return 6 * m.LPCAllTraysZV[tray] + 2 * m.LPCAllTraysb1[tray] >= 0

        test_model.LPCAllTraysZV_Con2 = Constraint(test_model.LPCTrays, rule=LPCAllTraysZV_Con2)

        def LPCAllTraysZLBL_Calc(m, tray):
            return exp(m.LPCAllTraysZLBL[tray]) - (m.LPCAllTraysZL[tray] - m.LPCAllTraysBL[tray]) == 0

        test_model.LPCAllTraysZLBL_Calc = Constraint(test_model.LPCTrays, rule=LPCAllTraysZLBL_Calc)

        def LPCAllTraysZLdivBL_Calc(m, tray):
            return exp(m.LPCAllTraysZLdivBL[tray]) - (m.LPCAllTraysZL[tray] + 2.414 * m.LPCAllTraysBL[tray]) / (
                    m.LPCAllTraysZL[tray] - 0.414 * m.LPCAllTraysBL[tray]) == 0

        test_model.LPCAllTraysZLdivBL_Calc = Constraint(test_model.LPCTrays, rule=LPCAllTraysZLdivBL_Calc)

        def LPCAllTraysZVBV_Calc(m, tray):
            return exp(m.LPCAllTraysZVBV[tray]) - (m.LPCAllTraysZV[tray] - m.LPCAllTraysBV[tray]) == 0

        test_model.LPCAllTraysZVBV_Calc = Constraint(test_model.LPCTrays, rule=LPCAllTraysZVBV_Calc)

        def LPCAllTraysZVdivBV_Calc(m, tray):
            return exp(m.LPCAllTraysZVdivBV[tray]) - (m.LPCAllTraysZV[tray] + 2.414 * m.LPCAllTraysBV[tray]) / (
                    m.LPCAllTraysZV[tray] - 0.414 * m.LPCAllTraysBV[tray]) == 0

        test_model.LPCAllTraysZVdivBV_Calc = Constraint(test_model.LPCTrays, rule=LPCAllTraysZVdivBV_Calc)

        def LPCAllTraysK_Calc(m, tray, comp):
            return m.LPCAllTraysPhiL[tray, comp] / m.LPCAllTraysPhiV[tray, comp] - m.LPCAllTraysVLE_K[tray, comp] == 0

        test_model.LPCAllTraysK_Calc = Constraint(test_model.LPCTrays, test_model.Component, rule=LPCAllTraysK_Calc)

        def LPCAllTraysVLEEqu(m, tray, comp):
            return m.LPCLiqLvMFrac[tray, comp] * m.LPCAllTraysVLE_K[tray, comp] - m.LPCVapLvMFrac[tray, comp] == 0

        test_model.LPCAllTraysVLEEqu = Constraint(test_model.LPCTrays, test_model.Component, rule=LPCAllTraysVLEEqu)

        # --------LPC Summation---------
        def LPCLiqSum(m, tray):
            return sum([m.LPCLiqLvMFrac[tray, c] for c in m.Component]) == 1

        test_model.LPCLiqSum = Constraint(test_model.LPCTrays, rule=LPCLiqSum)

        def LPCVapSum(m, tray):
            return sum([m.LPCVapLvMFrac[tray, c] for c in m.Component]) == 1

        test_model.LPCVapSum = Constraint(test_model.LPCTrays, rule=LPCVapSum)

        # --------LPC Pressure Profile---------
        def LPCPProf(m, tray):
            return ((m.LPCTopPressure - m.LPCBtmPressure) / 51 * (tray - 1) + m.LPCBtmPressure - m.LPCTrayPressure[
                tray]) * 0.010000 == 0

        test_model.LPCPProf = Constraint(test_model.LPCTrays, rule=LPCPProf)

        # ----------------------------------
        #           LPCReboiler
        # ----------------------------------
        # --------LPCReboiler Mass Balance---------
        def LPCReboilerMassBlnc(m, comp):
            return (m.LPCLiqLvMFlow[1] * m.LPCLiqLvMFrac[1, comp] - \
                    m.LPCReboilerLiqLvMFlow * m.LPCReboilerLiqLvMFrac[comp] - \
                    m.LPCReboilerVapLvMFlow * m.LPCReboilerVapLvMFrac[comp]) * 0.001000 == 0

        test_model.LPCReboilerMassBlnc = Constraint(test_model.Component, rule=LPCReboilerMassBlnc)

        # --------LPCReboiler Energy Balance---------
        def LPCReboilerEngBlnc(m):
            return (m.LPCLiqLvMFlow[1] * m.LPCLiqLvMEtlp[1] + m.HPCCondMHeatOut - \
                    m.LPCReboilerLiqLvMFlow * m.LPCReboilerLiqLvMEtlp - \
                    m.LPCReboilerVapLvMFlow * m.LPCReboilerVapLvMEtlp) * 0.000010 == 0

        test_model.LPCReboilerEngBlnc = Constraint(rule=LPCReboilerEngBlnc)

        # --------LPCReboiler Phase Equilibrium---------
        def LPCReboilerHdlpSqrtTr_Calc(m, comp):
            return m.LPCReboilerHdlpSqrtTr[comp] ** 2 - m.LPCReboilerHdlpTr[comp] == 0

        test_model.LPCReboilerHdlpSqrtTr_Calc = Constraint(test_model.Component, rule=LPCReboilerHdlpSqrtTr_Calc)

        def LPCReboilerHdlpSqrt_ai_Calc(m, comp1, comp2):
            return m.LPCReboilerHdlpSqrt_ai[comp1, comp2] ** 2 - m.LPCReboilerHdlpai[comp1] * m.LPCReboilerHdlpai[
                comp2] == 0

        test_model.LPCReboilerHdlpSqrt_ai_Calc = Constraint(test_model.Component, test_model.Component,
                                                            rule=LPCReboilerHdlpSqrt_ai_Calc)

        def LPCReboilerHdlpZL_Calc(m):
            return m.LPCReboilerHdlpZL ** 3 + m.LPCReboilerHdlpa1 * m.LPCReboilerHdlpZL ** 2 + \
                m.LPCReboilerHdlpa2 * m.LPCReboilerHdlpZL + m.LPCReboilerHdlpa3 == 0

        test_model.LPCReboilerHdlpZL_Calc = Constraint(rule=LPCReboilerHdlpZL_Calc)

        def LPCReboilerHdlpZL_Con1(m):
            return 3 * m.LPCReboilerHdlpZL ** 2 + 2 * m.LPCReboilerHdlpa1 * m.LPCReboilerHdlpZL + m.LPCReboilerHdlpa2 >= 0

        test_model.LPCReboilerHdlpZL_Con1 = Constraint(rule=LPCReboilerHdlpZL_Con1)

        def LPCReboilerHdlpZL_Con2(m):
            return 6 * m.LPCReboilerHdlpZL + 2 * m.LPCReboilerHdlpa1 <= 0

        test_model.LPCReboilerHdlpZL_Con2 = Constraint(rule=LPCReboilerHdlpZL_Con2)

        def LPCReboilerHdlpZV_Calc(m):
            return m.LPCReboilerHdlpZV ** 3 + m.LPCReboilerHdlpb1 * m.LPCReboilerHdlpZV ** 2 + \
                m.LPCReboilerHdlpb2 * m.LPCReboilerHdlpZV + m.LPCReboilerHdlpb3 == 0

        test_model.LPCReboilerHdlpZV_Calc = Constraint(rule=LPCReboilerHdlpZV_Calc)

        def LPCReboilerHdlpZV_Con1(m):
            return 3 * m.LPCReboilerHdlpZV ** 2 + 2 * m.LPCReboilerHdlpb1 * m.LPCReboilerHdlpZV + m.LPCReboilerHdlpb2 >= 0

        test_model.LPCReboilerHdlpZV_Con1 = Constraint(rule=LPCReboilerHdlpZV_Con1)

        def LPCReboilerHdlpZV_Con2(m):
            return 6 * m.LPCReboilerHdlpZV + 2 * m.LPCReboilerHdlpb1 >= 0

        test_model.LPCReboilerHdlpZV_Con2 = Constraint(rule=LPCReboilerHdlpZV_Con2)

        def LPCReboilerHdlpZLBL_Calc(m):
            return exp(m.LPCReboilerHdlpZLBL) - (m.LPCReboilerHdlpZL - m.LPCReboilerHdlpBL) == 0

        test_model.LPCReboilerHdlpZLBL_Calc = Constraint(rule=LPCReboilerHdlpZLBL_Calc)

        def LPCReboilerHdlpZLdivBL_Calc(m):
            return exp(m.LPCReboilerHdlpZLdivBL) - (m.LPCReboilerHdlpZL + 2.414 * m.LPCReboilerHdlpBL) / (
                    m.LPCReboilerHdlpZL - 0.414 * m.LPCReboilerHdlpBL) == 0

        test_model.LPCReboilerHdlpZLdivBL_Calc = Constraint(rule=LPCReboilerHdlpZLdivBL_Calc)

        def LPCReboilerHdlpZVBV_Calc(m):
            return exp(m.LPCReboilerHdlpZVBV) - (m.LPCReboilerHdlpZV - m.LPCReboilerHdlpBV) == 0

        test_model.LPCReboilerHdlpZVBV_Calc = Constraint(rule=LPCReboilerHdlpZVBV_Calc)

        def LPCReboilerHdlpZVdivBV_Calc(m):
            return exp(m.LPCReboilerHdlpZVdivBV) - (m.LPCReboilerHdlpZV + 2.414 * m.LPCReboilerHdlpBV) / (
                    m.LPCReboilerHdlpZV - 0.414 * m.LPCReboilerHdlpBV) == 0

        test_model.LPCReboilerHdlpZVdivBV_Calc = Constraint(rule=LPCReboilerHdlpZVdivBV_Calc)

        def LPCReboilerHdlpK_Calc(m, comp):
            return m.LPCReboilerHdlpPhiL[comp] / m.LPCReboilerHdlpPhiV[comp] - m.LPCReboilerHdlpVLE_K[comp] == 0

        test_model.LPCReboilerHdlpK_Calc = Constraint(test_model.Component, rule=LPCReboilerHdlpK_Calc)

        def LPCReboilerHdlpVLEEqu(m, comp):
            return m.LPCReboilerLiqLvMFrac[comp] * m.LPCReboilerHdlpVLE_K[comp] - m.LPCReboilerVapLvMFrac[comp] == 0

        test_model.LPCReboilerHdlpVLEEqu = Constraint(test_model.Component, rule=LPCReboilerHdlpVLEEqu)

        # --------LPCReboiler Summation---------
        def LPCReboilerLiqSum(m):
            return sum([m.LPCReboilerLiqLvMFrac[c] for c in m.Component]) == 1

        test_model.LPCReboilerLiqSum = Constraint(rule=LPCReboilerLiqSum)

        def LPCReboilerVapSum(m):
            return sum([m.LPCReboilerVapLvMFrac[c] for c in m.Component]) == 1

        test_model.LPCReboilerVapSum = Constraint(rule=LPCReboilerVapSum)

        # ----------------------------------
        #           ZeroReboiledASC
        # ----------------------------------

        # ----------------------------------
        #           ASC
        # ----------------------------------
        # --------ASC Mass Balance---------
        def ASCMassBlnc(m, tray, comp):
            if tray == 1:
                return (m.ZeroReboiledASCNull * m.ZeroReboiledASCNull + m.ASCLiqLvMCompFlow[tray + 1, comp] -
                        m.ASCLiqLvMCompFlow[tray, comp] - m.ASCVapLvMCompFlow[tray, comp] + m.LPCExt15MFlow *
                        m.LPCExt15MFrac[comp]) * 0.010000 == 0
            elif tray == 190:
                return (m.ASCCondLiqMFlow * m.ASCCondLiqMFrac[comp] + m.ASCVapLvMCompFlow[tray - 1, comp] -
                        m.ASCLiqLvMCompFlow[tray, comp] - m.ASCVapLvMCompFlow[tray, comp]) * 0.010000 == 0
            else:
                return (m.ASCLiqLvMCompFlow[tray + 1, comp] + m.ASCVapLvMCompFlow[tray - 1, comp] - m.ASCLiqLvMCompFlow[
                    tray, comp] - m.ASCVapLvMCompFlow[tray, comp]) * 0.010000 == 0

        test_model.ASCMassBlnc = Constraint(test_model.ASCTrays, test_model.Component, rule=ASCMassBlnc)

        # --------ASC Energy Balance---------
        def ASCEngBlnc(m, tray):
            if tray == 1:
                return (m.ZeroReboiledASCNull * m.ZeroReboiledASCNull \
                        + m.ASCLiqLvMFlow[tray + 1] * m.ASCLiqLvMEtlp[tray + 1] \
                        - m.ASCLiqLvMFlow[tray] * m.ASCLiqLvMEtlp[tray] \
                        - m.ASCVapLvMFlow[tray] * m.ASCVapLvMEtlp[tray] \
                        + m.LPCExt15MFlow * m.LPCExt15MEtlp) * 0.000100 == 0
            elif tray == 190:
                return (m.ASCCondLiqMFlow * m.ASCCondLiqMEtlp \
                        + m.ASCVapLvMFlow[tray - 1] * m.ASCVapLvMEtlp[tray - 1] \
                        - m.ASCLiqLvMFlow[tray] * m.ASCLiqLvMEtlp[tray] \
                        - m.ASCVapLvMFlow[tray] * m.ASCVapLvMEtlp[tray]) * 0.000100 == 0
            else:
                return (m.ASCLiqLvMFlow[tray + 1] * m.ASCLiqLvMEtlp[tray + 1] \
                        + m.ASCVapLvMFlow[tray - 1] * m.ASCVapLvMEtlp[tray - 1] \
                        - m.ASCLiqLvMFlow[tray] * m.ASCLiqLvMEtlp[tray] \
                        - m.ASCVapLvMFlow[tray] * m.ASCVapLvMEtlp[tray]) * 0.000100 == 0

        test_model.ASCEngBlnc = Constraint(test_model.ASCTrays, rule=ASCEngBlnc)

        # --------ASC Phase Equilibrium & System Parts---------
        def ASCAllTraysSqrtTr_Calc(m, tray, comp):
            return m.ASCAllTraysSqrtTr[tray, comp] ** 2 - m.ASCAllTraysTr[tray, comp] == 0

        test_model.ASCAllTraysSqrtTr_Calc = Constraint(test_model.ASCTrays, test_model.Component,
                                                       rule=ASCAllTraysSqrtTr_Calc)

        def ASCAllTraysSqrt_ai_Calc(m, tray, comp1, comp2):
            return m.ASCAllTraysSqrt_ai[tray, comp1, comp2] ** 2 - m.SRK_ai[comp1] * m.ASCAllTraysalpha[tray, comp1] * \
                m.SRK_ai[comp2] * m.ASCAllTraysalpha[tray, comp2] == 0

        test_model.ASCAllTraysSqrt_ai_Calc = Constraint(test_model.ASCTrays, test_model.Component, test_model.Component,
                                                        rule=ASCAllTraysSqrt_ai_Calc)

        def ASCAllTraysZL_Calc(m, tray):
            return m.ASCAllTraysZL[tray] ** 3 - m.ASCAllTraysZL[tray] ** 2 + \
                m.ASCAllTraysa1[tray] * m.ASCAllTraysZL[tray] + m.ASCAllTraysa2[tray] == 0

        test_model.ASCAllTraysZL_Calc = Constraint(test_model.ASCTrays, rule=ASCAllTraysZL_Calc)

        def ASCAllTraysZL_Con1(m, tray):
            return 3 * m.ASCAllTraysZL[tray] ** 2 - 2 * m.ASCAllTraysZL[tray] + m.ASCAllTraysa1[tray] >= 0

        test_model.ASCAllTraysZL_Con1 = Constraint(test_model.ASCTrays, rule=ASCAllTraysZL_Con1)

        def ASCAllTraysZL_Con2(m, tray):
            return 6 * m.ASCAllTraysZL[tray] - 2 <= 0

        test_model.ASCAllTraysZL_Con2 = Constraint(test_model.ASCTrays, rule=ASCAllTraysZL_Con2)

        def ASCAllTraysZV_Calc(m, tray):
            return m.ASCAllTraysZV[tray] ** 3 - m.ASCAllTraysZV[tray] ** 2 + \
                m.ASCAllTraysa1[tray] * m.ASCAllTraysZV[tray] + m.ASCAllTraysa2[tray] == 0

        test_model.ASCAllTraysZV_Calc = Constraint(test_model.ASCTrays, rule=ASCAllTraysZV_Calc)

        def ASCAllTraysZV_Con1(m, tray):
            return 3 * m.ASCAllTraysZV[tray] ** 2 - 2 * m.ASCAllTraysZV[tray] + m.ASCAllTraysa1[tray] >= 0

        test_model.ASCAllTraysZV_Con1 = Constraint(test_model.ASCTrays, rule=ASCAllTraysZV_Con1)

        def ASCAllTraysZV_Con2(m, tray):
            return 6 * m.ASCAllTraysZV[tray] - 2 >= 0

        test_model.ASCAllTraysZV_Con2 = Constraint(test_model.ASCTrays, rule=ASCAllTraysZV_Con2)

        def ASCAllTraysZLB_Calc(m, tray):
            return exp(m.ASCAllTraysZLB[tray]) - (m.ASCAllTraysZL[tray] - m.ASCAllTraysB[tray]) == 0

        test_model.ASCAllTraysZLB_Calc = Constraint(test_model.ASCTrays, rule=ASCAllTraysZLB_Calc)

        def ASCAllTraysBZL_Calc(m, tray):
            return exp(m.ASCAllTraysBZL[tray]) - (1 + m.ASCAllTraysB[tray] / m.ASCAllTraysZL[tray]) == 0

        test_model.ASCAllTraysBZL_Calc = Constraint(test_model.ASCTrays, rule=ASCAllTraysBZL_Calc)

        def ASCAllTraysZVB_Calc(m, tray):
            return exp(m.ASCAllTraysZVB[tray]) - (m.ASCAllTraysZV[tray] - m.ASCAllTraysB[tray]) == 0

        test_model.ASCAllTraysZVB_Calc = Constraint(test_model.ASCTrays, rule=ASCAllTraysZVB_Calc)

        def ASCAllTraysBZV_Calc(m, tray):
            return exp(m.ASCAllTraysBZV[tray]) - (1 + m.ASCAllTraysB[tray] / m.ASCAllTraysZV[tray]) == 0

        test_model.ASCAllTraysBZV_Calc = Constraint(test_model.ASCTrays, rule=ASCAllTraysBZV_Calc)

        def ASCAllTraysK_Calc(m, tray, comp):
            return m.ASCAllTraysPhiL[tray, comp] / m.ASCAllTraysPhiV[tray, comp] - m.ASCAllTraysVLE_K[tray, comp] == 0

        test_model.ASCAllTraysK_Calc = Constraint(test_model.ASCTrays, test_model.Component, rule=ASCAllTraysK_Calc)

        def ASCAllTraysVLEEqu(m, tray, comp):
            return m.ASCLiqLvMFrac[tray, comp] * m.ASCAllTraysVLE_K[tray, comp] - m.ASCVapLvMFrac[tray, comp] == 0

        test_model.ASCAllTraysVLEEqu = Constraint(test_model.ASCTrays, test_model.Component, rule=ASCAllTraysVLEEqu)

        def ASCAllTrayssqrtaiTr_Calc(m, tray, comp):
            return m.ASCAllTrayssqrtaiTr[tray, comp] ** 2 - m.SRK_ai[comp] * m.ASCAllTraysTr[tray, comp] == 0

        test_model.ASCAllTrayssqrtaiTr_Calc = Constraint(test_model.ASCTrays, test_model.Component,
                                                         rule=ASCAllTrayssqrtaiTr_Calc)

        def ASCAllTrayssqrtam_Calc(m, tray):
            return m.ASCAllTrayssqrtam[tray] ** 2 - m.ASCAllTraysam[tray] == 0

        test_model.ASCAllTrayssqrtam_Calc = Constraint(test_model.ASCTrays, rule=ASCAllTrayssqrtam_Calc)

        # --------ASC Summation---------
        def ASCLiqSum(m, tray):
            return sum([m.ASCLiqLvMFrac[tray, c] for c in m.Component]) == 1

        test_model.ASCLiqSum = Constraint(test_model.ASCTrays, rule=ASCLiqSum)

        def ASCVapSum(m, tray):
            return sum([m.ASCVapLvMFrac[tray, c] for c in m.Component]) == 1

        test_model.ASCVapSum = Constraint(test_model.ASCTrays, rule=ASCVapSum)

        # --------ASC Pressure Profile---------
        def ASCPProf(m, tray):
            return ((m.ASCTopPressure - m.ASCBtmPressure) / 189 * (tray - 1) + m.ASCBtmPressure - m.ASCTrayPressure[
                tray]) * 0.010000 == 0

        test_model.ASCPProf = Constraint(test_model.ASCTrays, rule=ASCPProf)

        # ----------------------------------
        #           ASCSump
        # ----------------------------------
        # --------ASCSump Mass Balance---------
        def ASCSumpMassBlnc(m):
            return (m.ASCLiqLvMFlow[1] - m.ASCSumpOutMFlow) * 0.010000 == 0

        test_model.ASCSumpMassBlnc = Constraint(rule=ASCSumpMassBlnc)

        def ASCSumpHldpMFracSpec(m, comp):
            return m.ASCSumpHldpMFrac[comp] - m.ASCLiqLvMFrac[1, comp] == 0

        test_model.ASCSumpHldpMFracSpec = Constraint(test_model.Component, rule=ASCSumpHldpMFracSpec)

        # ----------------------------------
        #           ASCCond
        # ----------------------------------
        # --------ASCCond Mass Balance---------
        def ASCCondMassBlnc(m, comp):
            return (m.ASCVapLvMFlow[190] * m.ASCVapLvMFrac[190, comp] - \
                    m.ASCCondLiqMFlow * m.ASCCondLiqMFrac[comp] - \
                    m.ASCCondVapMFlow * m.ASCCondVapMFrac[comp]) * 1.000000e-02 == 0

        test_model.ASCCondMassBlnc = Constraint(test_model.Component, rule=ASCCondMassBlnc)

        def ASCCondRefSpec(m):
            return (m.ASCCondLiqMFlow - m.ASCCondVapMFlow * m.ASCCondRefRatio) * 1.000000e-02 == 0

        test_model.ASCCondRefSpec = Constraint(rule=ASCCondRefSpec)

        # --------ASCCond Energy Balance---------
        def ASCCondEngBlnc(m):
            return (m.ASCVapLvMEtlp[190] * m.ASCVapLvMFlow[190] - \
                    m.ASCCondLiqMEtlp * m.ASCCondLiqMFlow - \
                    m.ASCCondVapMEtlp * m.ASCCondVapMFlow - m.ASCCondMHeat) * 1.000000e-04 == 0

        test_model.ASCCondEngBlnc = Constraint(rule=ASCCondEngBlnc)

        # --------ASCCond Phase Equilibrium---------
        def ASCCondSelfSqrtTr_Calc(m, comp):
            return m.ASCCondSelfSqrtTr[comp] ** 2 - m.ASCCondSelfTr[comp] == 0

        test_model.ASCCondSelfSqrtTr_Calc = Constraint(test_model.Component, rule=ASCCondSelfSqrtTr_Calc)

        def ASCCondSelfSqrt_ai_Calc(m, comp1, comp2):
            return m.ASCCondSelfSqrt_ai[comp1, comp2] ** 2 - m.SRK_ai[comp1] * m.ASCCondSelfalpha[comp1] * m.SRK_ai[
                comp2] * \
                m.ASCCondSelfalpha[comp2] == 0

        test_model.ASCCondSelfSqrt_ai_Calc = Constraint(test_model.Component, test_model.Component,
                                                        rule=ASCCondSelfSqrt_ai_Calc)

        def ASCCondSelfZL_Calc(m):
            return m.ASCCondSelfZL ** 3 - m.ASCCondSelfZL ** 2 + \
                m.ASCCondSelfa1 * m.ASCCondSelfZL + m.ASCCondSelfa2 == 0

        test_model.ASCCondSelfZL_Calc = Constraint(rule=ASCCondSelfZL_Calc)

        def ASCCondSelfZL_Con1(m):
            return 3 * m.ASCCondSelfZL ** 2 - 2 * m.ASCCondSelfZL + m.ASCCondSelfa1 >= 0

        test_model.ASCCondSelfZL_Con1 = Constraint(rule=ASCCondSelfZL_Con1)

        def ASCCondSelfZL_Con2(m):
            return 6 * m.ASCCondSelfZL - 2 <= 0

        test_model.ASCCondSelfZL_Con2 = Constraint(rule=ASCCondSelfZL_Con2)

        def ASCCondSelfZV_Calc(m):
            return m.ASCCondSelfZV ** 3 - m.ASCCondSelfZV ** 2 + \
                m.ASCCondSelfa1 * m.ASCCondSelfZV + m.ASCCondSelfa2 == 0

        test_model.ASCCondSelfZV_Calc = Constraint(rule=ASCCondSelfZV_Calc)

        def ASCCondSelfZV_Con1(m):
            return 3 * m.ASCCondSelfZV ** 2 - 2 * m.ASCCondSelfZV + m.ASCCondSelfa1 >= 0

        test_model.ASCCondSelfZV_Con1 = Constraint(rule=ASCCondSelfZV_Con1)

        def ASCCondSelfZV_Con2(m):
            return 6 * m.ASCCondSelfZV - 2 >= 0

        test_model.ASCCondSelfZV_Con2 = Constraint(rule=ASCCondSelfZV_Con2)

        def ASCCondSelfZLB_Calc(m):
            return exp(m.ASCCondSelfZLB) - (m.ASCCondSelfZL - m.ASCCondSelfB) == 0

        test_model.ASCCondSelfZLB_Calc = Constraint(rule=ASCCondSelfZLB_Calc)

        def ASCCondSelfBZL_Calc(m):
            return exp(m.ASCCondSelfBZL) - (1 + m.ASCCondSelfB / m.ASCCondSelfZL) == 0

        test_model.ASCCondSelfBZL_Calc = Constraint(rule=ASCCondSelfBZL_Calc)

        def ASCCondSelfZVB_Calc(m):
            return exp(m.ASCCondSelfZVB) - (m.ASCCondSelfZV - m.ASCCondSelfB) == 0

        test_model.ASCCondSelfZVB_Calc = Constraint(rule=ASCCondSelfZVB_Calc)

        def ASCCondSelfBZV_Calc(m):
            return exp(m.ASCCondSelfBZV) - (1 + m.ASCCondSelfB / m.ASCCondSelfZV) == 0

        test_model.ASCCondSelfBZV_Calc = Constraint(rule=ASCCondSelfBZV_Calc)

        def ASCCondSelfK_Calc(m, comp):
            return m.ASCCondSelfPhiL[comp] / m.ASCCondSelfPhiV[comp] - m.ASCCondSelfVLE_K[comp] == 0

        test_model.ASCCondSelfK_Calc = Constraint(test_model.Component, rule=ASCCondSelfK_Calc)

        def ASCCondSelfVLEEqu(m, comp):
            return m.ASCCondLiqMFrac[comp] * m.ASCCondSelfVLE_K[comp] - m.ASCCondVapMFrac[comp] == 0

        test_model.ASCCondSelfVLEEqu = Constraint(test_model.Component, rule=ASCCondSelfVLEEqu)

        def ASCCondSelfsqrtaiTr_Calc(m, comp):
            return m.ASCCondSelfsqrtaiTr[comp] ** 2 - m.SRK_ai[comp] * m.ASCCondSelfTr[comp] == 0

        test_model.ASCCondSelfsqrtaiTr_Calc = Constraint(test_model.Component, rule=ASCCondSelfsqrtaiTr_Calc)

        def ASCCondSelfsqrtam_Calc(m):
            return m.ASCCondSelfsqrtam ** 2 - m.ASCCondSelfam == 0

        test_model.ASCCondSelfsqrtam_Calc = Constraint(rule=ASCCondSelfsqrtam_Calc)

        # --------ASCCond Summation---------
        def ASCCondLiqSum(m):
            return sum([m.ASCCondLiqMFrac[c] for c in m.Component]) == 1

        test_model.ASCCondLiqSum = Constraint(rule=ASCCondLiqSum)

        def ASCCondVapSum(m):
            return sum([m.ASCCondVapMFrac[c] for c in m.Component]) == 1

        test_model.ASCCondVapSum = Constraint(rule=ASCCondVapSum)

        # ----------------------------------
        #           HeatGan
        # ----------------------------------
        # --------HeatGan DeltaT or DeltaP---------
        def HeatGanDeltaPSpec(m):
            return (m.LPCTrayPressure[52] - m.HeatGanDeltaP - m.HeatGanPressure) * 1e-2 == 0

        test_model.HeatGanDeltaPSpec = Constraint(rule=HeatGanDeltaPSpec)

        # --------HeatGan Enthalpy---------
        def HeatGanSelfSqrtTr_Calc(m, comp):
            return m.HeatGanSelfSqrtTr[comp] ** 2 - m.HeatGanSelfTr[comp] == 0

        test_model.HeatGanSelfSqrtTr_Calc = Constraint(test_model.Component, rule=HeatGanSelfSqrtTr_Calc)

        def HeatGanSelfSqrt_ai_Calc(m, comp1, comp2):
            return m.HeatGanSelfSqrt_ai[comp1, comp2] ** 2 - m.HeatGanSelfai[comp1] * m.HeatGanSelfai[comp2] == 0

        test_model.HeatGanSelfSqrt_ai_Calc = Constraint(test_model.Component, test_model.Component,
                                                        rule=HeatGanSelfSqrt_ai_Calc)

        def HeatGanSelfZV_Calc(m):
            return m.HeatGanSelfZV ** 3 + m.HeatGanSelfb1 * m.HeatGanSelfZV ** 2 + \
                m.HeatGanSelfb2 * m.HeatGanSelfZV + m.HeatGanSelfb3 == 0

        test_model.HeatGanSelfZV_Calc = Constraint(rule=HeatGanSelfZV_Calc)

        def HeatGanSelfZV_Con1(m):
            return 3 * m.HeatGanSelfZV ** 2 + 2 * m.HeatGanSelfb1 * m.HeatGanSelfZV + m.HeatGanSelfb2 >= 0

        test_model.HeatGanSelfZV_Con1 = Constraint(rule=HeatGanSelfZV_Con1)

        def HeatGanSelfZV_Con2(m):
            return 6 * m.HeatGanSelfZV + 2 * m.HeatGanSelfb1 >= 0

        test_model.HeatGanSelfZV_Con2 = Constraint(rule=HeatGanSelfZV_Con2)

        def HeatGanSelfZVdivBV_Calc(m):
            return exp(m.HeatGanSelfZVdivBV) - (m.HeatGanSelfZV + 2.414 * m.HeatGanSelfBV) / (
                    m.HeatGanSelfZV - 0.414 * m.HeatGanSelfBV) == 0

        test_model.HeatGanSelfZVdivBV_Calc = Constraint(rule=HeatGanSelfZVdivBV_Calc)

        # --------HeatGan Energy Balance---------
        def HeatGanEngBlnc(m):
            return ((m.LPCVapLvMEtlp[52] - m.HeatGanMEtlp) * m.LPCVapLvMFlow[52] + m.CoolLinMHeat) * 1.000000e-04 == 0

        test_model.HeatGanEngBlnc = Constraint(rule=HeatGanEngBlnc)

        # ----------------------------------
        #           RecoGan
        # ----------------------------------
        # --------RecoGan Enthalpy---------
        def RecoGanSelfSqrtTr_Calc(m, comp):
            return m.RecoGanSelfSqrtTr[comp] ** 2 - m.RecoGanSelfTr[comp] == 0

        test_model.RecoGanSelfSqrtTr_Calc = Constraint(test_model.Component, rule=RecoGanSelfSqrtTr_Calc)

        def RecoGanSelfSqrt_ai_Calc(m, comp1, comp2):
            return m.RecoGanSelfSqrt_ai[comp1, comp2] ** 2 - m.RecoGanSelfai[comp1] * m.RecoGanSelfai[comp2] == 0

        test_model.RecoGanSelfSqrt_ai_Calc = Constraint(test_model.Component, test_model.Component,
                                                        rule=RecoGanSelfSqrt_ai_Calc)

        def RecoGanSelfZV_Calc(m):
            return m.RecoGanSelfZV ** 3 + m.RecoGanSelfb1 * m.RecoGanSelfZV ** 2 + \
                m.RecoGanSelfb2 * m.RecoGanSelfZV + m.RecoGanSelfb3 == 0

        test_model.RecoGanSelfZV_Calc = Constraint(rule=RecoGanSelfZV_Calc)

        def RecoGanSelfZV_Con1(m):
            return 3 * m.RecoGanSelfZV ** 2 + 2 * m.RecoGanSelfb1 * m.RecoGanSelfZV + m.RecoGanSelfb2 >= 0

        test_model.RecoGanSelfZV_Con1 = Constraint(rule=RecoGanSelfZV_Con1)

        def RecoGanSelfZV_Con2(m):
            return 6 * m.RecoGanSelfZV + 2 * m.RecoGanSelfb1 >= 0

        test_model.RecoGanSelfZV_Con2 = Constraint(rule=RecoGanSelfZV_Con2)

        def RecoGanSelfZVdivBV_Calc(m):
            return exp(m.RecoGanSelfZVdivBV) - (m.RecoGanSelfZV + 2.414 * m.RecoGanSelfBV) / (
                    m.RecoGanSelfZV - 0.414 * m.RecoGanSelfBV) == 0

        test_model.RecoGanSelfZVdivBV_Calc = Constraint(rule=RecoGanSelfZVdivBV_Calc)

        # --------RecoGan Energy Balance---------
        def RecoGanEngBlnc(m):
            return ((m.HeatGanMEtlp - m.RecoGanMEtlp) * m.LPCVapLvMFlow[52] + m.RecoGanMHeat) * 1.000000e-04 == 0

        test_model.RecoGanEngBlnc = Constraint(rule=RecoGanEngBlnc)

        # ----------------------------------
        #           RecoGar
        # ----------------------------------
        # --------RecoGar Enthalpy---------
        def RecoGarSelfSqrtTr_Calc(m, comp):
            return m.RecoGarSelfSqrtTr[comp] ** 2 - m.RecoGarSelfTr[comp] == 0

        test_model.RecoGarSelfSqrtTr_Calc = Constraint(test_model.Component, rule=RecoGarSelfSqrtTr_Calc)

        def RecoGarSelfSqrt_ai_Calc(m, comp1, comp2):
            return m.RecoGarSelfSqrt_ai[comp1, comp2] ** 2 - m.RecoGarSelfai[comp1] * m.RecoGarSelfai[comp2] == 0

        test_model.RecoGarSelfSqrt_ai_Calc = Constraint(test_model.Component, test_model.Component,
                                                        rule=RecoGarSelfSqrt_ai_Calc)

        def RecoGarSelfZV_Calc(m):
            return m.RecoGarSelfZV ** 3 + m.RecoGarSelfb1 * m.RecoGarSelfZV ** 2 + \
                m.RecoGarSelfb2 * m.RecoGarSelfZV + m.RecoGarSelfb3 == 0

        test_model.RecoGarSelfZV_Calc = Constraint(rule=RecoGarSelfZV_Calc)

        def RecoGarSelfZV_Con1(m):
            return 3 * m.RecoGarSelfZV ** 2 + 2 * m.RecoGarSelfb1 * m.RecoGarSelfZV + m.RecoGarSelfb2 >= 0

        test_model.RecoGarSelfZV_Con1 = Constraint(rule=RecoGarSelfZV_Con1)

        def RecoGarSelfZV_Con2(m):
            return 6 * m.RecoGarSelfZV + 2 * m.RecoGarSelfb1 >= 0

        test_model.RecoGarSelfZV_Con2 = Constraint(rule=RecoGarSelfZV_Con2)

        def RecoGarSelfZVdivBV_Calc(m):
            return exp(m.RecoGarSelfZVdivBV) - (m.RecoGarSelfZV + 2.414 * m.RecoGarSelfBV) / (
                    m.RecoGarSelfZV - 0.414 * m.RecoGarSelfBV) == 0

        test_model.RecoGarSelfZVdivBV_Calc = Constraint(rule=RecoGarSelfZVdivBV_Calc)

        # --------RecoGar Energy Balance---------
        def RecoGarEngBlnc(m):
            return ((m.ASCCondVapMEtlp - m.RecoGarMEtlp) * m.ASCCondVapMFlow + m.RecoGarMHeat) * 1.000000e-04 == 0

        test_model.RecoGarEngBlnc = Constraint(rule=RecoGarEngBlnc)

        # ----------------------------------
        #           RecoWn
        # ----------------------------------
        # --------RecoWn Enthalpy---------
        def RecoWnSelfSqrtTr_Calc(m, comp):
            return m.RecoWnSelfSqrtTr[comp] ** 2 - m.RecoWnSelfTr[comp] == 0

        test_model.RecoWnSelfSqrtTr_Calc = Constraint(test_model.Component, rule=RecoWnSelfSqrtTr_Calc)

        def RecoWnSelfSqrt_ai_Calc(m, comp1, comp2):
            return m.RecoWnSelfSqrt_ai[comp1, comp2] ** 2 - m.RecoWnSelfai[comp1] * m.RecoWnSelfai[comp2] == 0

        test_model.RecoWnSelfSqrt_ai_Calc = Constraint(test_model.Component, test_model.Component,
                                                       rule=RecoWnSelfSqrt_ai_Calc)

        def RecoWnSelfZV_Calc(m):
            return m.RecoWnSelfZV ** 3 + m.RecoWnSelfb1 * m.RecoWnSelfZV ** 2 + \
                m.RecoWnSelfb2 * m.RecoWnSelfZV + m.RecoWnSelfb3 == 0

        test_model.RecoWnSelfZV_Calc = Constraint(rule=RecoWnSelfZV_Calc)

        def RecoWnSelfZV_Con1(m):
            return 3 * m.RecoWnSelfZV ** 2 + 2 * m.RecoWnSelfb1 * m.RecoWnSelfZV + m.RecoWnSelfb2 >= 0

        test_model.RecoWnSelfZV_Con1 = Constraint(rule=RecoWnSelfZV_Con1)

        def RecoWnSelfZV_Con2(m):
            return 6 * m.RecoWnSelfZV + 2 * m.RecoWnSelfb1 >= 0

        test_model.RecoWnSelfZV_Con2 = Constraint(rule=RecoWnSelfZV_Con2)

        def RecoWnSelfZVdivBV_Calc(m):
            return exp(m.RecoWnSelfZVdivBV) - (m.RecoWnSelfZV + 2.414 * m.RecoWnSelfBV) / (
                    m.RecoWnSelfZV - 0.414 * m.RecoWnSelfBV) == 0

        test_model.RecoWnSelfZVdivBV_Calc = Constraint(rule=RecoWnSelfZVdivBV_Calc)

        # --------RecoWn Energy Balance---------
        def RecoWnEngBlnc(m):
            return ((m.LPCExt46MEtlp - m.RecoWnMEtlp) * m.LPCExt46MFlow + m.RecoWnMHeat) * 1.000000e-04 == 0

        test_model.RecoWnEngBlnc = Constraint(rule=RecoWnEngBlnc)

        # ----------------------------------
        #           OxSplitter
        # ----------------------------------
        # --------OxSplitter Mass Balance---------
        def OxSplitterOutASpec(m):
            return (m.OxSplitterOutAMFlow - m.LPCReboilerLiqLvMFlow * m.OxSplitterOutARatio) * 0.001000 == 0

        test_model.OxSplitterOutASpec = Constraint(rule=OxSplitterOutASpec)

        def OxSplitterOutBSpec(m):
            return (m.OxSplitterOutBMFlow - m.LPCReboilerLiqLvMFlow * (1 - m.OxSplitterOutARatio)) * 0.001000 == 0

        test_model.OxSplitterOutBSpec = Constraint(rule=OxSplitterOutBSpec)

        # ----------------------------------
        #           RecoGox
        # ----------------------------------
        # --------RecoGox Enthalpy---------
        def RecoGoxSelfSqrtTr_Calc(m, comp):
            return m.RecoGoxSelfSqrtTr[comp] ** 2 - m.RecoGoxSelfTr[comp] == 0

        test_model.RecoGoxSelfSqrtTr_Calc = Constraint(test_model.Component, rule=RecoGoxSelfSqrtTr_Calc)

        def RecoGoxSelfSqrt_ai_Calc(m, comp1, comp2):
            return m.RecoGoxSelfSqrt_ai[comp1, comp2] ** 2 - m.RecoGoxSelfai[comp1] * m.RecoGoxSelfai[comp2] == 0

        test_model.RecoGoxSelfSqrt_ai_Calc = Constraint(test_model.Component, test_model.Component,
                                                        rule=RecoGoxSelfSqrt_ai_Calc)

        def RecoGoxSelfZV_Calc(m):
            return m.RecoGoxSelfZV ** 3 + m.RecoGoxSelfb1 * m.RecoGoxSelfZV ** 2 + \
                m.RecoGoxSelfb2 * m.RecoGoxSelfZV + m.RecoGoxSelfb3 == 0

        test_model.RecoGoxSelfZV_Calc = Constraint(rule=RecoGoxSelfZV_Calc)

        def RecoGoxSelfZV_Con1(m):
            return 3 * m.RecoGoxSelfZV ** 2 + 2 * m.RecoGoxSelfb1 * m.RecoGoxSelfZV + m.RecoGoxSelfb2 >= 0

        test_model.RecoGoxSelfZV_Con1 = Constraint(rule=RecoGoxSelfZV_Con1)

        def RecoGoxSelfZV_Con2(m):
            return 6 * m.RecoGoxSelfZV + 2 * m.RecoGoxSelfb1 >= 0

        test_model.RecoGoxSelfZV_Con2 = Constraint(rule=RecoGoxSelfZV_Con2)

        def RecoGoxSelfZVdivBV_Calc(m):
            return exp(m.RecoGoxSelfZVdivBV) - (m.RecoGoxSelfZV + 2.414 * m.RecoGoxSelfBV) / (
                    m.RecoGoxSelfZV - 0.414 * m.RecoGoxSelfBV) == 0

        test_model.RecoGoxSelfZVdivBV_Calc = Constraint(rule=RecoGoxSelfZVdivBV_Calc)

        # --------RecoGox Energy Balance---------
        def RecoGoxEngBlnc(m):
            return ((
                                m.LPCReboilerLiqLvMEtlp - m.RecoGoxMEtlp) * m.OxSplitterOutBMFlow + m.RecoGoxMHeat) * 1.000000e-04 == 0

        test_model.RecoGoxEngBlnc = Constraint(rule=RecoGoxEngBlnc)

    def build_rto(self, test_model, cv_func):
        LOX = 17000

        def PurityCon1(m):
            return 0.999 - m.LPCReboilerLiqLvMFrac['Oxygen']

        test_model.PurityCon1 = Expression(rule=PurityCon1)

        def PurityCon2(m):
            return m.LPCReboilerLiqLvMFrac['Nitrogen'] - 8e-6

        test_model.PurityCon2 = Expression(rule=PurityCon2)

        def PurityCon3(m):
            return 0.9999 - m.LPCVapLvMFrac[52, 'Nitrogen']

        test_model.PurityCon3 = Expression(rule=PurityCon3)

        def PurityCon4(m):
            return m.LPCVapLvMFrac[52, 'Oxygen'] - 8e-6

        test_model.PurityCon4 = Expression(rule=PurityCon4)

        def PurityCon5(m):
            return m.ASCCondVapMFrac['Oxygen'] - 8e-6

        test_model.PurityCon5 = Expression(rule=PurityCon5)

        def PurityCon6(m):
            return m.LPCVapLvMFrac[15, 'Nitrogen'] - 8e-4

        test_model.PurityCon6 = Expression(rule=PurityCon6)

        # Product Flowrate - adjusted according to LOX load
        def ProductCon1(m):
            return LOX / 22.4 - m.OxSplitterOutBMFlow

        test_model.ProductCon1 = Expression(rule=ProductCon1)

        def ProductCon2(m):
            return 16000 * 2 / 22.4 - m.LPCVapLvMFlow[52]

        test_model.ProductCon2 = Expression(rule=ProductCon2)

        def ProductCon3(m):
            return 16000 / 32 / 22.4 - m.ASCCondVapMFlow

        test_model.ProductCon3 = Expression(rule=ProductCon3)

        # Temperature Difference
        def TempDiffCon1(m):
            return 1.5 - (m.HPCTrayTemp[41] - m.LPCReboilerTemp)

        test_model.TempDiffCon1 = Expression(rule=TempDiffCon1)

        def TempDiffCon2(m):
            return 3 - (m.ASCTrayTemp[190] - m.HeatLPCTemp)

        test_model.TempDiffCon2 = Expression(rule=TempDiffCon2)

        # Turbine
        def TurbineCon1(m):
            return 1.32315 * 0.8 * m.OxSplitterOutBMFlow - m.FeedSplitterOut3MFlow

        test_model.TurbineCon1 = Expression(rule=TurbineCon1)

        def TurbineCon2(m):
            return m.FeedSplitterOut3MFlow - 1.32315 * 1.2 * m.OxSplitterOutBMFlow

        test_model.TurbineCon2 = Expression(rule=TurbineCon2)

        def TurbineCon3(m):
            return 0.6 * 4461.574 - (m.FeedSplitterOut3MFlow + m.FeedSplitterOut1MFlow + m.FeedSplitterOut2MFlow)

        test_model.TurbineCon3 = Expression(rule=TurbineCon3)

        def TurbineCon4(m):
            return m.FeedSplitterOut3MFlow + m.FeedSplitterOut1MFlow + m.FeedSplitterOut2MFlow - 1.5 * 4461.574

        test_model.TurbineCon4 = Expression(rule=TurbineCon4)

        def TurbineCon7(m):
            return 0.6 * 1236.75 - m.FeedSplitterOut3MFlow

        test_model.TurbineCon7 = Expression(rule=TurbineCon7)

        def TurbineCon8(m):
            return m.FeedSplitterOut3MFlow - 1.5 * 1236.75

        test_model.TurbineCon8 = Expression(rule=TurbineCon8)

        def TurbineCon9(m):
            return 0.6 * 1236.3 - m.FeedSplitterOut1MFlow

        test_model.TurbineCon9 = Expression(rule=TurbineCon9)

        def TurbineCon10(m):
            return m.FeedSplitterOut1MFlow - 1.5 * 1236.3

        test_model.TurbineCon10 = Expression(rule=TurbineCon10)

        # Overall Heat Balance
        def OverallHeatBlnc(m):
            return m.RecoGanMHeat + m.RecoWnMHeat + m.RecoGarMHeat + m.RecoGoxMHeat - m.TaCoolerMHeat - m.MaCoolerMHeat - m.HpaCoolerMHeat

        test_model.OverallHeatBlnc = Expression(rule=OverallHeatBlnc)

        # Objective
        def ObjExpr(m):
            return -(16.32 * m.OxSplitterOutAMFlow) + 0.8 * (
                    1.70919 * m.FeedSplitterOut2MFlow + 4.06931 * m.FeedSplitterOut3MFlow + 2.92443 * m.FeedSplitterOut1MFlow)

        test_model.OBJ = Expression(rule=ObjExpr)