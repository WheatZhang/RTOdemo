from rtolib.core.pyomo_model import PyomoModel
from pyomo.environ import *
import os

class RTO_Plant_HPC(PyomoModel):
    def __init__(self):
        self.output_variables = {
            'Purity_GAN': (lambda m: m.ColumnVapLvMFrac[60, 'Nitrogen']),
            'T_GAN': (lambda m: m.ColumnTrayTemp[60]),
            'T_30thTray': (lambda m: m.ColumnTrayTemp[30]),
            'WN_Flow':(lambda m: m.ReboilerVapLvMFlow),
            'Drain_Flow':(lambda m: m.ReboilerLiqLvMFlow),
            'obj': (lambda m: m.obj),
            "feed_con1": (lambda m: m.feed_con1),
            "feed_con2": (lambda m: m.feed_con2),
            "purity_con": (lambda m: m.purity_con_noise_free),
            "drain_con": (lambda m: m.drain_con_noise_free),
        }
        self.noised_outputs = {
            'Purity_GAN': (lambda m: m.ColumnVapLvMFrac[60, 'Nitrogen']+m.Purity_GAN_noise),
            'T_GAN': (lambda m: m.ColumnTrayTemp[60]+m.T_GAN_noise),
            'T_30thTray': (lambda m: m.ColumnTrayTemp[30]+m.T_30thTray_noise),
            'WN_Flow':(lambda m: m.ReboilerVapLvMFlow+m.WN_Flow_noise),
            'Drain_Flow':(lambda m: m.ReboilerLiqLvMFlow+m.Drain_Flow_noise),
            'obj': (lambda m: m.obj),
            "feed_con1": (lambda m: m.feed_con1),
            "feed_con2": (lambda m: m.feed_con2),
            "purity_con": (lambda m: m.purity_con),
            "drain_con": (lambda m: m.drain_con),
        }
        self.input_variables = {
            'TA_Flow': (lambda m: m.FeedTaMFlow),
            'MA_Flow': (lambda m: m.FeedMaMFlow),
            'GAN_Flow': (lambda m: m.GNSplitterOutAMFlow),
            'LN_Flow': (lambda m: m.CondensorPrdtMFlow),
        }
        self.parameters = {
        }
        self.default_value={}
        self.parameter_scaling_factors={}
        self.output_noise = {
            'Purity_GAN': (lambda m: m.Purity_GAN_noise),
            'T_GAN': (lambda m: m.T_GAN_noise),
            'T_30thTray': (lambda m: m.T_30thTray_noise),
            'WN_Flow':(lambda m: m.WN_Flow_noise),
            'Drain_Flow':(lambda m: m.Drain_Flow_noise),
        }
        self.initial_value_file = os.path.join(os.path.dirname(__file__) + r"\NominalOpt_InitValue.txt")

    def build_body(self, hpc_model):
        hpc_model.Component = Set(initialize=['Oxygen', 'Nitrogen', 'Argon'])

        # ===================================
        #
        #      Process Parts
        #
        # ===================================
        hpc_model.PR_CrtcTemp = Param(hpc_model.Component,
                                       initialize={'Oxygen': 154.77, 'Argon': 150.71,
                                                   'Nitrogen': 126.15})  # Temp unit: K
        hpc_model.PR_CrtcPressure = Param(hpc_model.Component, initialize={'Oxygen': 5080000, 'Argon': 4864000,
                                                                             'Nitrogen': 3394000})  # Crtc pressure, pa
        hpc_model.PR_CompOmega = Param(hpc_model.Component,
                                        initialize={'Oxygen': 0.019, 'Argon': 0, 'Nitrogen': 0.045})
        hpc_model.PR_SatPresCorC = Param(hpc_model.Component, initialize={'Oxygen': -5.667, 'Nitrogen': -6.344})
        dict = {
            ('Oxygen', 'Oxygen'): 0,
            ('Oxygen', 'Argon'): 0.0265,
            ('Oxygen', 'Nitrogen'): -0.01238,
            ('Argon', 'Oxygen'): 0.0265,
            ('Argon', 'Argon'): 0,
            ('Argon', 'Nitrogen'): -0.004071,
            ('Nitrogen', 'Oxygen'): -0.01238,
            ('Nitrogen', 'Argon'): -0.004071,
            ('Nitrogen', 'Nitrogen'): 0
        }
        hpc_model.PR_Kij = Var(hpc_model.Component, hpc_model.Component, initialize=dict)
        hpc_model.PR_Kij['Oxygen', 'Oxygen'].fixed = True
        hpc_model.PR_Kij['Oxygen', 'Argon'].fixed = True
        hpc_model.PR_Kij['Argon', 'Oxygen'].fixed = True
        hpc_model.PR_Kij['Argon', 'Argon'].fixed = True
        hpc_model.PR_Kij['Argon', 'Nitrogen'].fixed = True
        hpc_model.PR_Kij['Nitrogen', 'Argon'].fixed = True
        hpc_model.PR_Kij['Nitrogen', 'Nitrogen'].fixed = True
        hpc_model.PR_RefTemp = Param(hpc_model.Component,
                                      initialize={'Oxygen': 298.15, 'Argon': 298.15, 'Nitrogen': 298.15})  # unit: K
        # Reference state enthalpy for ideal gas at reference temperature, cal/mol
        hpc_model.PR_DHFORM = Param(hpc_model.Component, initialize={'Oxygen': 0, 'Argon': 0, 'Nitrogen': 0})
        # Parameters to calculate the equation for the DIPPR ideal gas heat capacity model by Aly and Lee 1981,unit: cal/(mol*k)
        hpc_model.PR_CPIGDP1 = Param(hpc_model.Component,
                                      initialize={'Oxygen': 1 / 4180 * 29103, 'Argon': 1 / 4180 * 20786,
                                                  'Nitrogen': 1 / 4180 * 29105})
        hpc_model.PR_CPIGDP2 = Param(hpc_model.Component,
                                      initialize={'Oxygen': 1 / 4180 * 10040, 'Argon': 1 / 4180 * 0,
                                                  'Nitrogen': 1 / 4180 * 8614.9})
        hpc_model.PR_CPIGDP3 = Param(hpc_model.Component,
                                      initialize={'Oxygen': 2526.5, 'Argon': 0, 'Nitrogen': 1701.6})
        hpc_model.PR_CPIGDP4 = Param(hpc_model.Component,
                                      initialize={'Oxygen': 1 / 4180 * 9356, 'Argon': 1 / 4180 * 0,
                                                  'Nitrogen': 1 / 4180 * 103.47})
        hpc_model.PR_CPIGDP5 = Param(hpc_model.Component,
                                      initialize={'Oxygen': 1153.8, 'Argon': 0, 'Nitrogen': 909.79})
        hpc_model.PR_CPIGDP6 = Param(hpc_model.Component, initialize={'Oxygen': 50, 'Argon': -173.15, 'Nitrogen': 50})
        hpc_model.PR_CPIGDP7 = Param(hpc_model.Component,
                                      initialize={'Oxygen': 1500, 'Argon': 1226.85, 'Nitrogen': 1500})
        hpc_model.PR_RGas = Param(initialize=8.314)  # J/(mol*K)

        def PR_bi_Calc(m, comp):
            return 0.07780 * m.PR_RGas * m.PR_CrtcTemp[comp] / m.PR_CrtcPressure[comp]

        hpc_model.PR_bi = Param(hpc_model.Component, initialize=PR_bi_Calc)

        def PR_m_Calc(m, comp):
            return 0.37464 + 1.54226 * m.PR_CompOmega[comp] - 0.26992 * m.PR_CompOmega[comp] * m.PR_CompOmega[comp]

        hpc_model.PR_m = Param(hpc_model.Component, initialize=PR_m_Calc)

        # ===================================
        #
        #      Parameters & Sets
        #
        # ===================================
        # --------Feed---------
        hpc_model.FeedTaMFrac = Param(hpc_model.Component,
                                       initialize={'Nitrogen': 0.78118, 'Oxygen': 0.2095, 'Argon': 0.00932})
        hpc_model.FeedTaPressure = Param(initialize=100.000000)
        hpc_model.FeedTaTemp = Param(initialize=39.850000)
        hpc_model.FeedTaMEtlp = Param(initialize=424.928400)
        # --------FeedMa---------
        hpc_model.FeedMaMFrac = Param(hpc_model.Component,
                                       initialize={'Nitrogen': 0.78118, 'Oxygen': 0.2095, 'Argon': 0.00932})
        hpc_model.FeedMaPressure = Param(initialize=100.000000)
        hpc_model.FeedMaTemp = Param(initialize=39.850000)
        hpc_model.FeedMaMEtlp = Param(initialize=424.928400)

        # --------TaCooler---------
        hpc_model.TaCoolerPressure = Param(initialize=566.000000)
        hpc_model.TaCoolerVF = Param(initialize=0.400000)
        # --------MaCooler---------
        hpc_model.MaCoolerTemp = Param(initialize=-168.550000)
        hpc_model.MaCoolerPressure = Param(initialize=566.000000)
        # --------Column---------
        hpc_model.ColumnTrays = RangeSet(1, 60)
        hpc_model.ColumnTopPressure = Param(initialize=539.560000)
        hpc_model.ColumnBtmPressure = Param(initialize=564.000000)
        # --------Sump---------
        hpc_model.SumpLiqRho = Param(initialize=1000.000000)
        hpc_model.SumpSumpCSArea = Param(initialize=6.000000)
        # --------Condensor---------
        hpc_model.CondensorPressure = Param(initialize=538.600000)
        # --------Reboiler---------
        hpc_model.ReboilerLiqRho = Param(initialize=30282.000000)
        hpc_model.ReboilerRebCSArea = Param(initialize=6.000000)
        # --------Reboiler VLE & Enthalpy---------
        # --------GNSplitter---------
        hpc_model.GNSplitterOutAMFlow = Var(initialize=800, bounds=(0, 2000))
        hpc_model.GNSplitterOutAMFlow.fixed = True
        hpc_model.CondensorPrdtMFlow = Var(initialize=100, bounds=(0, 200))
        hpc_model.CondensorPrdtMFlow.fixed = True

        # --------HPCZeroReboiled---------
        hpc_model.HPCZeroReboiledNull = Param(initialize=0.000000)

        # =================================================================================================
        # =================================================================================================
        # ======================================
        #
        #        Model Parameters
        #
        # ======================================
        hpc_model.ColumnEffTray = Param(initialize=0.7)
        hpc_model.NitrogenOxygenKij = Param(initialize=-0.01238)
        hpc_model.ColumnUA = Param(initialize=2)
        # ======================================
        #
        #        Output noise
        #
        # ======================================
        hpc_model.Purity_GAN_noise = Var(initialize=0)
        hpc_model.Purity_GAN_noise.fixed=True
        hpc_model.T_GAN_noise = Var(initialize=0)
        hpc_model.T_GAN_noise.fixed=True
        hpc_model.T_30thTray_noise = Var(initialize=0)
        hpc_model.T_30thTray_noise.fixed=True
        hpc_model.WN_Flow_noise = Var(initialize=0)
        hpc_model.WN_Flow_noise.fixed=True
        hpc_model.Drain_Flow_noise = Var(initialize=0)
        hpc_model.Drain_Flow_noise.fixed=True
        # ======================================
        #
        #        Manipulated Variables
        #
        # ======================================
        hpc_model.FeedTaMFlow = Var(initialize=320.000000, bounds=(0, 2000))
        hpc_model.FeedTaMFlow.fixed = True
        hpc_model.FeedMaMFlow = Var(initialize=1680.000000, bounds=(0, 3000))
        hpc_model.FeedMaMFlow.fixed = True
        # =================================================================================================
        # =================================================================================================

        # ===================================
        #
        #         Variables
        #
        # ===================================
        # --------TaCooler---------
        hpc_model.TaCoolerTemp = Var()
        hpc_model.TaCoolerVapMFlow = Var(within=NonNegativeReals)
        hpc_model.TaCoolerLiqMFlow = Var(within=NonNegativeReals)
        hpc_model.TaCoolerMHeat = Var(within=NonNegativeReals)
        hpc_model.TaCoolerLiqMFrac = Var(hpc_model.Component, bounds=(0, 1))
        hpc_model.TaCoolerVapMFrac = Var(hpc_model.Component, bounds=(0, 1))
        hpc_model.TaCoolerSelfSqrtTr = Var(hpc_model.Component, initialize=0.8)
        hpc_model.TaCoolerSelfSqrt_ai = Var(hpc_model.Component, hpc_model.Component, initialize=0.17,
                                             bounds=(0.05, 0.4))
        hpc_model.TaCoolerSelfZL = Var(initialize=0.023966, bounds=(0, 0.05))
        hpc_model.TaCoolerSelfZV = Var(initialize=0.857334, bounds=(0.2, 3))
        hpc_model.TaCoolerSelfZLBL = Var(initialize=-4.878610, bounds=(-10, -1))
        hpc_model.TaCoolerSelfZLdivBL = Var(initialize=1.30898, bounds=(0.8, 2))
        hpc_model.TaCoolerSelfZVBV = Var(initialize=-0.1732, bounds=(-0.3, 0))
        hpc_model.TaCoolerSelfZVdivBV = Var(initialize=0.053666900, bounds=(0, 0.2))
        hpc_model.TaCoolerSelfVLE_K = Var(hpc_model.Component, bounds=(0, 8))
        # --------MaCooler---------
        hpc_model.MaCoolerMHeat = Var(within=NonNegativeReals)
        hpc_model.MaCoolerSelfSqrtTr = Var(hpc_model.Component, initialize=0.8)
        hpc_model.MaCoolerSelfSqrt_ai = Var(hpc_model.Component, hpc_model.Component, initialize=0.17,
                                             bounds=(0.05, 0.4))
        hpc_model.MaCoolerSelfZV = Var(initialize=0.866923293, bounds=(0.2, 3))
        hpc_model.MaCoolerSelfZVdivBV = Var(initialize=0.050441702, bounds=(0, 0.2))
        # --------AirMixer---------
        hpc_model.AirMixerOutMFlow = Var(within=NonNegativeReals)
        hpc_model.AirMixerOutMFrac = Var(hpc_model.Component, bounds=(0, 1))
        hpc_model.AirMixerOutMEtlp = Var()
        # --------Column---------
        hpc_model.ColumnVapLvMFlow = Var(hpc_model.ColumnTrays, within=NonNegativeReals)
        hpc_model.ColumnLiqLvMFlow = Var(hpc_model.ColumnTrays, within=NonNegativeReals)
        hpc_model.ColumnTrayTemp = Var(hpc_model.ColumnTrays)
        hpc_model.ColumnTrayPressure = Var(hpc_model.ColumnTrays)
        hpc_model.ColumnLiqLvMFrac = Var(hpc_model.ColumnTrays, hpc_model.Component, bounds=(0, 1))
        hpc_model.ColumnVapEqilMFrac = Var(hpc_model.ColumnTrays, hpc_model.Component, bounds=(0, 1))
        hpc_model.ColumnAllTraysSqrtTr = Var(hpc_model.ColumnTrays, hpc_model.Component, initialize=0.8)
        hpc_model.ColumnAllTraysSqrt_ai = Var(hpc_model.ColumnTrays, hpc_model.Component, hpc_model.Component,
                                               initialize=0.17, bounds=(0.05, 0.4))
        hpc_model.ColumnAllTraysZL = Var(hpc_model.ColumnTrays, initialize=0.023966, bounds=(0, 0.05))
        hpc_model.ColumnAllTraysZV = Var(hpc_model.ColumnTrays, initialize=0.857334, bounds=(0.2, 3))
        hpc_model.ColumnAllTraysZLBL = Var(hpc_model.ColumnTrays, initialize=-4.878610, bounds=(-10, -1))
        hpc_model.ColumnAllTraysZLdivBL = Var(hpc_model.ColumnTrays, initialize=1.30898, bounds=(0.8, 2))
        hpc_model.ColumnAllTraysZVBV = Var(hpc_model.ColumnTrays, initialize=-0.1732, bounds=(-0.3, 0))
        hpc_model.ColumnAllTraysZVdivBV = Var(hpc_model.ColumnTrays, initialize=0.053666900, bounds=(0, 0.2))
        hpc_model.ColumnAllTraysVLE_K = Var(hpc_model.ColumnTrays, hpc_model.Component, bounds=(0, 8))
        hpc_model.ColumnVapLvMFrac = Var(hpc_model.ColumnTrays, hpc_model.Component, bounds=(0, 1))
        # --------Sump---------
        hpc_model.SumpOutMFlow = Var(within=NonNegativeReals)
        hpc_model.SumpHldpMFrac = Var(hpc_model.Component, bounds=(0, 1))
        # --------Condensor---------
        hpc_model.CondensorRefMFlow = Var(within=NonNegativeReals)
        hpc_model.CondensorMHeatOut = Var(within=NonNegativeReals)
        hpc_model.CondensorOutletSqrtTr = Var(hpc_model.Component, initialize=0.8)
        hpc_model.CondensorOutletSqrt_ai = Var(hpc_model.Component, hpc_model.Component, initialize=0.17,
                                                bounds=(0.05, 0.4))
        hpc_model.CondensorOutletZL = Var(initialize=0.023916, bounds=(0, 0.05))
        hpc_model.CondensorOutletZLdivBL = Var(initialize=1.3089, bounds=(0.8, 2))
        # --------Reboiler---------
        hpc_model.ReboilerLiqLvMFlow = Var(within=NonNegativeReals)
        hpc_model.ReboilerVapLvMFlow = Var(within=NonNegativeReals)
        hpc_model.ReboilerTemp = Var()
        hpc_model.ReboilerLiqLvMFrac = Var(hpc_model.Component, bounds=(0, 1))
        hpc_model.ReboilerVapLvMFrac = Var(hpc_model.Component, bounds=(0, 1))
        hpc_model.ReboilerHdlpSqrtTr = Var(hpc_model.Component, initialize=0.8)
        hpc_model.ReboilerHdlpSqrt_ai = Var(hpc_model.Component, hpc_model.Component, initialize=0.17,
                                             bounds=(0.05, 0.4))
        hpc_model.ReboilerHdlpZL = Var(initialize=0.023966, bounds=(0, 0.05))
        hpc_model.ReboilerHdlpZV = Var(initialize=0.857334, bounds=(0.2, 3))
        hpc_model.ReboilerHdlpZLBL = Var(initialize=-4.878610, bounds=(-10, -1))
        hpc_model.ReboilerHdlpZLdivBL = Var(initialize=1.30898, bounds=(0.8, 2))
        hpc_model.ReboilerHdlpZVBV = Var(initialize=-0.1732, bounds=(-0.3, 0))
        hpc_model.ReboilerHdlpZVdivBV = Var(initialize=0.053666900, bounds=(0, 0.2))
        hpc_model.ReboilerHdlpVLE_K = Var(hpc_model.Component, bounds=(0, 8))
        # --------GNSplitter---------
        hpc_model.GNSplitterOutBMFlow = Var(within=NonNegativeReals)

        # ===================================
        #
        #         Expressions
        #
        # ===================================
        # --------TaCooler Enthalpy---------
        def TaCoolerSelfTr_Calc(m, comp):
            return (m.TaCoolerTemp + 273.15) / m.PR_CrtcTemp[comp]

        hpc_model.TaCoolerSelfTr = Expression(hpc_model.Component, rule=TaCoolerSelfTr_Calc)

        def TaCoolerSelfalpha_Calc(m, comp):
            return (1 + m.PR_m[comp] * (1 - m.TaCoolerSelfSqrtTr[comp])) ** 2

        hpc_model.TaCoolerSelfalpha = Expression(hpc_model.Component, rule=TaCoolerSelfalpha_Calc)

        def TaCoolerSelfai_Calc(m, comp):
            return 0.45724 * m.TaCoolerSelfalpha[comp] * (m.PR_RGas ** 2) * (m.PR_CrtcTemp[comp] ** 2) / \
                   m.PR_CrtcPressure[comp]

        hpc_model.TaCoolerSelfai = Expression(hpc_model.Component, rule=TaCoolerSelfai_Calc)

        def TaCoolerSelfaij_Calc(m, comp1, comp2):
            return m.TaCoolerSelfSqrt_ai[comp1, comp2] * (1 - m.PR_Kij[comp1, comp2])

        hpc_model.TaCoolerSelfaij = Expression(hpc_model.Component, hpc_model.Component, rule=TaCoolerSelfaij_Calc)

        def TaCoolerSelfaL_Calc(m):
            return sum([sum(
                [m.TaCoolerSelfaij[c1, c2] * m.TaCoolerLiqMFrac[c1] * m.TaCoolerLiqMFrac[c2] for c1 in m.Component]) for
                        c2 in m.Component])

        hpc_model.TaCoolerSelfaL = Expression(rule=TaCoolerSelfaL_Calc)

        def TaCoolerSelfbL_Calc(m):
            return sum([m.PR_bi[c] * m.TaCoolerLiqMFrac[c] for c in m.Component])

        hpc_model.TaCoolerSelfbL = Expression(rule=TaCoolerSelfbL_Calc)

        def TaCoolerSelfAL_Calc(m):
            return m.TaCoolerSelfaL * m.TaCoolerPressure * 1000 / ((m.PR_RGas * (m.TaCoolerTemp + 273.15)) ** 2)

        hpc_model.TaCoolerSelfAL = Expression(rule=TaCoolerSelfAL_Calc)

        def TaCoolerSelfBL_Calc(m):
            return m.TaCoolerSelfbL * m.TaCoolerPressure * 1000 / (m.PR_RGas * (m.TaCoolerTemp + 273.15))

        hpc_model.TaCoolerSelfBL = Expression(rule=TaCoolerSelfBL_Calc)

        def TaCoolerSelfa1_Calc(m):
            return m.TaCoolerSelfBL - 1

        hpc_model.TaCoolerSelfa1 = Expression(rule=TaCoolerSelfa1_Calc)

        def TaCoolerSelfa2_Calc(m):
            return m.TaCoolerSelfAL - 3 * m.TaCoolerSelfBL ** 2 - 2 * m.TaCoolerSelfBL

        hpc_model.TaCoolerSelfa2 = Expression(rule=TaCoolerSelfa2_Calc)

        def TaCoolerSelfa3_Calc(m):
            return m.TaCoolerSelfBL ** 2 + m.TaCoolerSelfBL ** 3 - m.TaCoolerSelfAL * m.TaCoolerSelfBL

        hpc_model.TaCoolerSelfa3 = Expression(rule=TaCoolerSelfa3_Calc)

        def TaCoolerSelfaV_Calc(m):
            return sum([sum(
                [m.TaCoolerSelfaij[c1, c2] * m.TaCoolerVapMFrac[c1] * m.TaCoolerVapMFrac[c2] for c1 in m.Component]) for
                        c2 in m.Component])

        hpc_model.TaCoolerSelfaV = Expression(rule=TaCoolerSelfaV_Calc)

        def TaCoolerSelfbV_Calc(m):
            return sum([m.PR_bi[c] * m.TaCoolerVapMFrac[c] for c in m.Component])

        hpc_model.TaCoolerSelfbV = Expression(rule=TaCoolerSelfbV_Calc)

        def TaCoolerSelfAV_Calc(m):
            return m.TaCoolerSelfaV * m.TaCoolerPressure * 1000 / ((m.PR_RGas * (m.TaCoolerTemp + 273.15)) ** 2)

        hpc_model.TaCoolerSelfAV = Expression(rule=TaCoolerSelfAV_Calc)

        def TaCoolerSelfBV_Calc(m):
            return m.TaCoolerSelfbV * m.TaCoolerPressure * 1000 / (m.PR_RGas * (m.TaCoolerTemp + 273.15))

        hpc_model.TaCoolerSelfBV = Expression(rule=TaCoolerSelfBV_Calc)

        def TaCoolerSelfb1_Calc(m):
            return m.TaCoolerSelfBV - 1

        hpc_model.TaCoolerSelfb1 = Expression(rule=TaCoolerSelfb1_Calc)

        def TaCoolerSelfb2_Calc(m):
            return m.TaCoolerSelfAV - 3 * m.TaCoolerSelfBV ** 2 - 2 * m.TaCoolerSelfBV

        hpc_model.TaCoolerSelfb2 = Expression(rule=TaCoolerSelfb2_Calc)

        def TaCoolerSelfb3_Calc(m):
            return m.TaCoolerSelfBV ** 2 + m.TaCoolerSelfBV ** 3 - m.TaCoolerSelfAV * m.TaCoolerSelfBV

        hpc_model.TaCoolerSelfb3 = Expression(rule=TaCoolerSelfb3_Calc)

        def TaCoolerSelfSL_Calc(m, comp):
            return sum([m.TaCoolerSelfaij[c, comp] * m.TaCoolerLiqMFrac[c] for c in m.Component])

        hpc_model.TaCoolerSelfSL = Expression(hpc_model.Component, rule=TaCoolerSelfSL_Calc)

        def TaCoolerSelfSV_Calc(m, comp):
            return sum([m.TaCoolerSelfaij[c, comp] * m.TaCoolerVapMFrac[c] for c in m.Component])

        hpc_model.TaCoolerSelfSV = Expression(hpc_model.Component, rule=TaCoolerSelfSV_Calc)

        def TaCoolerSelfPhiL_Calc(m, comp):
            return exp(m.PR_bi[comp] / m.TaCoolerSelfbL * (m.TaCoolerSelfZL - 1) - \
                       m.TaCoolerSelfZLBL - m.TaCoolerSelfAL / m.TaCoolerSelfBL / 2 / sqrt(2) * (2 * \
                                                                                                 m.TaCoolerSelfSL[
                                                                                                     comp] / m.TaCoolerSelfaL -
                                                                                                 m.PR_bi[
                                                                                                     comp] / m.TaCoolerSelfbL) * \
                       m.TaCoolerSelfZLdivBL)

        hpc_model.TaCoolerSelfPhiL = Expression(hpc_model.Component, rule=TaCoolerSelfPhiL_Calc)

        def TaCoolerSelfPhiV_Calc(m, comp):
            return exp(m.PR_bi[comp] / m.TaCoolerSelfbV * (m.TaCoolerSelfZV - 1) - \
                       m.TaCoolerSelfZVBV - m.TaCoolerSelfAV / m.TaCoolerSelfBV / 2 / sqrt(2) * (2 * \
                                                                                                 m.TaCoolerSelfSV[
                                                                                                     comp] / m.TaCoolerSelfaV -
                                                                                                 m.PR_bi[
                                                                                                     comp] / m.TaCoolerSelfbV) * \
                       m.TaCoolerSelfZVdivBV)

        hpc_model.TaCoolerSelfPhiV = Expression(hpc_model.Component, rule=TaCoolerSelfPhiV_Calc)

        def TaCoolerSelfH0_Calc(m, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.TaCoolerTemp + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (m.PR_DHFORM[comp] + m.PR_CPIGDP1[comp] * (
                            (m.TaCoolerTemp + 273.15) - m.PR_RefTemp[comp]) + 2 * m.PR_CPIGDP2[comp] * \
                               m.PR_CPIGDP3[comp] * (-exp(2 / (m.TaCoolerTemp + 273.15) * m.PR_CPIGDP3[comp]) + \
                                                     exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp])) / (
                                           exp(2 / (m.TaCoolerTemp + 273.15) * m.PR_CPIGDP3[comp]) - 1) / (
                                           exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp]) - 1) - \
                               2 * m.PR_CPIGDP4[comp] * m.PR_CPIGDP5[comp] * (
                                           exp(2 / (m.TaCoolerTemp + 273.15) * m.PR_CPIGDP5[comp]) - exp(
                                       2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp])) / \
                               (exp(2 / (m.TaCoolerTemp + 273.15) * m.PR_CPIGDP5[comp]) + 1) / (
                                           exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp]) + 1))

        hpc_model.TaCoolerSelfH0 = Expression(hpc_model.Component, rule=TaCoolerSelfH0_Calc)

        def TaCoolerSelfM_Calc(m, comp1, comp2):
            return m.PR_m[comp1] * m.TaCoolerSelfSqrtTr[comp1] / (2 * sqrt(m.TaCoolerSelfalpha[comp1]))

        hpc_model.TaCoolerSelfM = Expression(hpc_model.Component, hpc_model.Component,
                                              rule=TaCoolerSelfM_Calc)  # Calculate the enthalpy of liquid phase mixture %unit: J/mol

        def TaCoolerSelfHV_Calc(m):
            return sum([m.TaCoolerVapMFrac[c] * m.TaCoolerSelfH0[c] for c in m.Component]) + \
                   m.PR_RGas * (m.TaCoolerTemp + 273.15) * ((m.TaCoolerSelfZV - 1) - 1 / (2 ** 1.5) * (
                        1 / m.TaCoolerSelfBV) * m.TaCoolerSelfZVdivBV * \
                                                            sum([sum([m.TaCoolerVapMFrac[c1] * m.TaCoolerVapMFrac[
                                                                c2] * m.TaCoolerSelfAV * (1 + m.TaCoolerSelfM[c2, c1] + \
                                                                                          m.TaCoolerSelfM[c1, c2]) for
                                                                      c1 in m.Component]) for c2 in m.Component]))

        hpc_model.TaCoolerVapMEtlp = Expression(rule=TaCoolerSelfHV_Calc)

        # Calculate the enthalpy of liquid phase mixture %unit: J/mol
        def TaCoolerSelfHL_Calc(m):
            return sum([m.TaCoolerLiqMFrac[c] * m.TaCoolerSelfH0[c] for c in m.Component]) + \
                   m.PR_RGas * (m.TaCoolerTemp + 273.15) * ((m.TaCoolerSelfZL - 1) - 1 / (2 ** 1.5) * (
                        1 / m.TaCoolerSelfBL) * m.TaCoolerSelfZLdivBL * \
                                                            sum([sum([m.TaCoolerLiqMFrac[c1] * m.TaCoolerLiqMFrac[
                                                                c2] * m.TaCoolerSelfAL * (1 + m.TaCoolerSelfM[c2, c1] + \
                                                                                          m.TaCoolerSelfM[c1, c2]) for
                                                                      c1 in m.Component]) for c2 in m.Component]))

        hpc_model.TaCoolerLiqMEtlp = Expression(rule=TaCoolerSelfHL_Calc)

        def TaCoolerMEtlp(m):
            return (
                               m.TaCoolerLiqMEtlp * m.TaCoolerLiqMFlow + m.TaCoolerVapMEtlp * m.TaCoolerVapMFlow) / m.FeedTaMFlow

        hpc_model.TaCoolerMEtlp = Expression(rule=TaCoolerMEtlp)

        # --------MaCooler Enthalpy---------
        def MaCoolerSelfTr_Calc(m, comp):
            return (m.MaCoolerTemp + 273.15) / m.PR_CrtcTemp[comp]

        hpc_model.MaCoolerSelfTr = Expression(hpc_model.Component, rule=MaCoolerSelfTr_Calc)

        def MaCoolerSelfalpha_Calc(m, comp):
            return (1 + m.PR_m[comp] * (1 - m.MaCoolerSelfSqrtTr[comp])) ** 2

        hpc_model.MaCoolerSelfalpha = Expression(hpc_model.Component, rule=MaCoolerSelfalpha_Calc)

        def MaCoolerSelfai_Calc(m, comp):
            return 0.45724 * m.MaCoolerSelfalpha[comp] * (m.PR_RGas ** 2) * (m.PR_CrtcTemp[comp] ** 2) / \
                   m.PR_CrtcPressure[comp]

        hpc_model.MaCoolerSelfai = Expression(hpc_model.Component, rule=MaCoolerSelfai_Calc)

        def MaCoolerSelfaij_Calc(m, comp1, comp2):
            return m.MaCoolerSelfSqrt_ai[comp1, comp2] * (1 - m.PR_Kij[comp1, comp2])

        hpc_model.MaCoolerSelfaij = Expression(hpc_model.Component, hpc_model.Component, rule=MaCoolerSelfaij_Calc)

        def MaCoolerSelfaV_Calc(m):
            return sum(
                [sum([m.MaCoolerSelfaij[c1, c2] * m.FeedMaMFrac[c1] * m.FeedMaMFrac[c2] for c1 in m.Component]) for c2 in
                 m.Component])

        hpc_model.MaCoolerSelfaV = Expression(rule=MaCoolerSelfaV_Calc)

        def MaCoolerSelfbV_Calc(m):
            return sum([m.PR_bi[c] * m.FeedMaMFrac[c] for c in m.Component])

        hpc_model.MaCoolerSelfbV = Expression(rule=MaCoolerSelfbV_Calc)

        def MaCoolerSelfAV_Calc(m):
            return m.MaCoolerSelfaV * m.MaCoolerPressure * 1000 / ((m.PR_RGas * (m.MaCoolerTemp + 273.15)) ** 2)

        hpc_model.MaCoolerSelfAV = Expression(rule=MaCoolerSelfAV_Calc)

        def MaCoolerSelfBV_Calc(m):
            return m.MaCoolerSelfbV * m.MaCoolerPressure * 1000 / (m.PR_RGas * (m.MaCoolerTemp + 273.15))

        hpc_model.MaCoolerSelfBV = Expression(rule=MaCoolerSelfBV_Calc)

        def MaCoolerSelfb1_Calc(m):
            return m.MaCoolerSelfBV - 1

        hpc_model.MaCoolerSelfb1 = Expression(rule=MaCoolerSelfb1_Calc)

        def MaCoolerSelfb2_Calc(m):
            return m.MaCoolerSelfAV - 3 * m.MaCoolerSelfBV ** 2 - 2 * m.MaCoolerSelfBV

        hpc_model.MaCoolerSelfb2 = Expression(rule=MaCoolerSelfb2_Calc)

        def MaCoolerSelfb3_Calc(m):
            return m.MaCoolerSelfBV ** 2 + m.MaCoolerSelfBV ** 3 - m.MaCoolerSelfAV * m.MaCoolerSelfBV

        hpc_model.MaCoolerSelfb3 = Expression(rule=MaCoolerSelfb3_Calc)

        def MaCoolerSelfH0_Calc(m, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.MaCoolerTemp + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (m.PR_DHFORM[comp] + m.PR_CPIGDP1[comp] * (
                        (m.MaCoolerTemp + 273.15) - m.PR_RefTemp[comp]) + 2 * m.PR_CPIGDP2[comp] * \
                               m.PR_CPIGDP3[comp] * (-exp(2 / (m.MaCoolerTemp + 273.15) * m.PR_CPIGDP3[comp]) + \
                                                     exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp])) / (
                                           exp(2 / (m.MaCoolerTemp + 273.15) * m.PR_CPIGDP3[comp]) - 1) / (
                                           exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp]) - 1) - \
                               2 * m.PR_CPIGDP4[comp] * m.PR_CPIGDP5[comp] * (
                                           exp(2 / (m.MaCoolerTemp + 273.15) * m.PR_CPIGDP5[comp]) - exp(
                                       2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp])) / \
                               (exp(2 / (m.MaCoolerTemp + 273.15) * m.PR_CPIGDP5[comp]) + 1) / (
                                           exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp]) + 1))

        hpc_model.MaCoolerSelfH0 = Expression(hpc_model.Component, rule=MaCoolerSelfH0_Calc)

        def MaCoolerSelfM_Calc(m, comp1, comp2):
            return m.PR_m[comp1] * m.MaCoolerSelfSqrtTr[comp1] / (2 * sqrt(m.MaCoolerSelfalpha[comp1]))

        hpc_model.MaCoolerSelfM = Expression(hpc_model.Component, hpc_model.Component,
                                              rule=MaCoolerSelfM_Calc)  # Calculate the enthalpy of liquid phase mixture %unit: J/mol

        def MaCoolerSelfHV_Calc(m):
            return sum([m.FeedMaMFrac[c] * m.MaCoolerSelfH0[c] for c in m.Component]) + \
                   m.PR_RGas * (m.MaCoolerTemp + 273.15) * ((m.MaCoolerSelfZV - 1) - 1 / (2 ** 1.5) * (
                        1 / m.MaCoolerSelfBV) * m.MaCoolerSelfZVdivBV * \
                                                            sum([sum([m.FeedMaMFrac[c1] * m.FeedMaMFrac[
                                                                c2] * m.MaCoolerSelfAV * (1 + m.MaCoolerSelfM[c2, c1] + \
                                                                                          m.MaCoolerSelfM[c1, c2]) for
                                                                      c1 in m.Component]) for c2 in m.Component]))

        hpc_model.MaCoolerMEtlp = Expression(rule=MaCoolerSelfHV_Calc)

        # --------Column Component Flowrate---------
        def ColumnVapLvMCompFlow(m, tray, comp):
            return m.ColumnVapLvMFrac[tray, comp] * m.ColumnVapLvMFlow[tray]

        hpc_model.ColumnVapLvMCompFlow = Expression(hpc_model.ColumnTrays, hpc_model.Component,
                                                     rule=ColumnVapLvMCompFlow)

        def ColumnLiqLvMCompFlow(m, tray, comp):
            return m.ColumnLiqLvMFrac[tray, comp] * m.ColumnLiqLvMFlow[tray]

        hpc_model.ColumnLiqLvMCompFlow = Expression(hpc_model.ColumnTrays, hpc_model.Component,
                                                     rule=ColumnLiqLvMCompFlow)

        # --------Column Thermo and Enthalpy---------
        def ColumnAllTraysTr_Calc(m, tray, comp):
            return (m.ColumnTrayTemp[tray] + 273.15) / m.PR_CrtcTemp[comp]

        hpc_model.ColumnAllTraysTr = Expression(hpc_model.ColumnTrays, hpc_model.Component,
                                                 rule=ColumnAllTraysTr_Calc)

        def ColumnAllTraysalpha_Calc(m, tray, comp):
            return (1 + m.PR_m[comp] * (1 - m.ColumnAllTraysSqrtTr[tray, comp])) ** 2

        hpc_model.ColumnAllTraysalpha = Expression(hpc_model.ColumnTrays, hpc_model.Component,
                                                    rule=ColumnAllTraysalpha_Calc)

        def ColumnAllTraysai_Calc(m, tray, comp):
            return 0.45724 * m.ColumnAllTraysalpha[tray, comp] * (m.PR_RGas ** 2) * (m.PR_CrtcTemp[comp] ** 2) / \
                   m.PR_CrtcPressure[comp]

        hpc_model.ColumnAllTraysai = Expression(hpc_model.ColumnTrays, hpc_model.Component,
                                                 rule=ColumnAllTraysai_Calc)

        def ColumnAllTraysaij_Calc(m, tray, comp1, comp2):
            return m.ColumnAllTraysSqrt_ai[tray, comp1, comp2] * (1 - m.PR_Kij[comp1, comp2])

        hpc_model.ColumnAllTraysaij = Expression(hpc_model.ColumnTrays, hpc_model.Component, hpc_model.Component,
                                                  rule=ColumnAllTraysaij_Calc)

        def ColumnAllTraysaL_Calc(m, tray):
            return sum([sum(
                [m.ColumnAllTraysaij[tray, c1, c2] * m.ColumnLiqLvMFrac[tray, c1] * m.ColumnLiqLvMFrac[tray, c2] for c1
                 in m.Component]) for c2 in m.Component])

        hpc_model.ColumnAllTraysaL = Expression(hpc_model.ColumnTrays, rule=ColumnAllTraysaL_Calc)

        def ColumnAllTraysbL_Calc(m, tray):
            return sum([m.PR_bi[c] * m.ColumnLiqLvMFrac[tray, c] for c in m.Component])

        hpc_model.ColumnAllTraysbL = Expression(hpc_model.ColumnTrays, rule=ColumnAllTraysbL_Calc)

        def ColumnAllTraysAL_Calc(m, tray):
            return m.ColumnAllTraysaL[tray] * m.ColumnTrayPressure[tray] * 1000 / (
                        (m.PR_RGas * (m.ColumnTrayTemp[tray] + 273.15)) ** 2)

        hpc_model.ColumnAllTraysAL = Expression(hpc_model.ColumnTrays, rule=ColumnAllTraysAL_Calc)

        def ColumnAllTraysBL_Calc(m, tray):
            return m.ColumnAllTraysbL[tray] * m.ColumnTrayPressure[tray] * 1000 / (
                        m.PR_RGas * (m.ColumnTrayTemp[tray] + 273.15))

        hpc_model.ColumnAllTraysBL = Expression(hpc_model.ColumnTrays, rule=ColumnAllTraysBL_Calc)

        def ColumnAllTraysa1_Calc(m, tray):
            return m.ColumnAllTraysBL[tray] - 1

        hpc_model.ColumnAllTraysa1 = Expression(hpc_model.ColumnTrays, rule=ColumnAllTraysa1_Calc)

        def ColumnAllTraysa2_Calc(m, tray):
            return m.ColumnAllTraysAL[tray] - 3 * m.ColumnAllTraysBL[tray] ** 2 - 2 * m.ColumnAllTraysBL[tray]

        hpc_model.ColumnAllTraysa2 = Expression(hpc_model.ColumnTrays, rule=ColumnAllTraysa2_Calc)

        def ColumnAllTraysa3_Calc(m, tray):
            return m.ColumnAllTraysBL[tray] ** 2 + m.ColumnAllTraysBL[tray] ** 3 - m.ColumnAllTraysAL[tray] * \
                   m.ColumnAllTraysBL[tray]

        hpc_model.ColumnAllTraysa3 = Expression(hpc_model.ColumnTrays, rule=ColumnAllTraysa3_Calc)

        def ColumnAllTraysaV_Calc(m, tray):
            return sum([sum(
                [m.ColumnAllTraysaij[tray, c1, c2] * m.ColumnVapEqilMFrac[tray, c1] * m.ColumnVapEqilMFrac[tray, c2] for
                 c1 in m.Component]) for c2 in m.Component])

        hpc_model.ColumnAllTraysaV = Expression(hpc_model.ColumnTrays, rule=ColumnAllTraysaV_Calc)

        def ColumnAllTraysbV_Calc(m, tray):
            return sum([m.PR_bi[c] * m.ColumnVapEqilMFrac[tray, c] for c in m.Component])

        hpc_model.ColumnAllTraysbV = Expression(hpc_model.ColumnTrays, rule=ColumnAllTraysbV_Calc)

        def ColumnAllTraysAV_Calc(m, tray):
            return m.ColumnAllTraysaV[tray] * m.ColumnTrayPressure[tray] * 1000 / (
                        (m.PR_RGas * (m.ColumnTrayTemp[tray] + 273.15)) ** 2)

        hpc_model.ColumnAllTraysAV = Expression(hpc_model.ColumnTrays, rule=ColumnAllTraysAV_Calc)

        def ColumnAllTraysBV_Calc(m, tray):
            return m.ColumnAllTraysbV[tray] * m.ColumnTrayPressure[tray] * 1000 / (
                        m.PR_RGas * (m.ColumnTrayTemp[tray] + 273.15))

        hpc_model.ColumnAllTraysBV = Expression(hpc_model.ColumnTrays, rule=ColumnAllTraysBV_Calc)

        def ColumnAllTraysb1_Calc(m, tray):
            return m.ColumnAllTraysBV[tray] - 1

        hpc_model.ColumnAllTraysb1 = Expression(hpc_model.ColumnTrays, rule=ColumnAllTraysb1_Calc)

        def ColumnAllTraysb2_Calc(m, tray):
            return m.ColumnAllTraysAV[tray] - 3 * m.ColumnAllTraysBV[tray] ** 2 - 2 * m.ColumnAllTraysBV[tray]

        hpc_model.ColumnAllTraysb2 = Expression(hpc_model.ColumnTrays, rule=ColumnAllTraysb2_Calc)

        def ColumnAllTraysb3_Calc(m, tray):
            return m.ColumnAllTraysBV[tray] ** 2 + m.ColumnAllTraysBV[tray] ** 3 - m.ColumnAllTraysAV[tray] * \
                   m.ColumnAllTraysBV[tray]

        hpc_model.ColumnAllTraysb3 = Expression(hpc_model.ColumnTrays, rule=ColumnAllTraysb3_Calc)

        def ColumnAllTraysSL_Calc(m, tray, comp):
            return sum([m.ColumnAllTraysaij[tray, c, comp] * m.ColumnLiqLvMFrac[tray, c] for c in m.Component])

        hpc_model.ColumnAllTraysSL = Expression(hpc_model.ColumnTrays, hpc_model.Component,
                                                 rule=ColumnAllTraysSL_Calc)

        def ColumnAllTraysSV_Calc(m, tray, comp):
            return sum([m.ColumnAllTraysaij[tray, c, comp] * m.ColumnVapEqilMFrac[tray, c] for c in m.Component])

        hpc_model.ColumnAllTraysSV = Expression(hpc_model.ColumnTrays, hpc_model.Component,
                                                 rule=ColumnAllTraysSV_Calc)

        def ColumnAllTraysPhiL_Calc(m, tray, comp):
            return exp(m.PR_bi[comp] / m.ColumnAllTraysbL[tray] * (m.ColumnAllTraysZL[tray] - 1) - \
                       m.ColumnAllTraysZLBL[tray] - m.ColumnAllTraysAL[tray] / m.ColumnAllTraysBL[tray] / 2 / sqrt(
                2) * (2 * \
                      m.ColumnAllTraysSL[tray, comp] / m.ColumnAllTraysaL[tray] - m.PR_bi[comp] / m.ColumnAllTraysbL[
                          tray]) * \
                       m.ColumnAllTraysZLdivBL[tray])

        hpc_model.ColumnAllTraysPhiL = Expression(hpc_model.ColumnTrays, hpc_model.Component,
                                                   rule=ColumnAllTraysPhiL_Calc)

        def ColumnAllTraysPhiV_Calc(m, tray, comp):
            return exp(m.PR_bi[comp] / m.ColumnAllTraysbV[tray] * (m.ColumnAllTraysZV[tray] - 1) - \
                       m.ColumnAllTraysZVBV[tray] - m.ColumnAllTraysAV[tray] / m.ColumnAllTraysBV[tray] / 2 / sqrt(
                2) * (2 * \
                      m.ColumnAllTraysSV[tray, comp] / m.ColumnAllTraysaV[tray] - m.PR_bi[comp] / m.ColumnAllTraysbV[
                          tray]) * \
                       m.ColumnAllTraysZVdivBV[tray])

        hpc_model.ColumnAllTraysPhiV = Expression(hpc_model.ColumnTrays, hpc_model.Component,
                                                   rule=ColumnAllTraysPhiV_Calc)

        def ColumnAllTraysH0_Calc(m, tray, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.ColumnTrayTemp[tray] + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (m.PR_DHFORM[comp] + m.PR_CPIGDP1[comp] * (
                            (m.ColumnTrayTemp[tray] + 273.15) - m.PR_RefTemp[comp]) + 2 * m.PR_CPIGDP2[comp] * \
                               m.PR_CPIGDP3[comp] * (-exp(2 / (m.ColumnTrayTemp[tray] + 273.15) * m.PR_CPIGDP3[comp]) + \
                                                     exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp])) / (
                                           exp(2 / (m.ColumnTrayTemp[tray] + 273.15) * m.PR_CPIGDP3[comp]) - 1) / (
                                           exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp]) - 1) - \
                               2 * m.PR_CPIGDP4[comp] * m.PR_CPIGDP5[comp] * (
                                           exp(2 / (m.ColumnTrayTemp[tray] + 273.15) * m.PR_CPIGDP5[comp]) - exp(
                                       2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp])) / \
                               (exp(2 / (m.ColumnTrayTemp[tray] + 273.15) * m.PR_CPIGDP5[comp]) + 1) / (
                                           exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp]) + 1))

        hpc_model.ColumnAllTraysH0 = Expression(hpc_model.ColumnTrays, hpc_model.Component,
                                                 rule=ColumnAllTraysH0_Calc)

        def ColumnAllTraysM_Calc(m, tray, comp1, comp2):
            return m.PR_m[comp1] * m.ColumnAllTraysSqrtTr[tray, comp1] / (2 * sqrt(m.ColumnAllTraysalpha[tray, comp1]))

        hpc_model.ColumnAllTraysM = Expression(hpc_model.ColumnTrays, hpc_model.Component, hpc_model.Component,
                                                rule=ColumnAllTraysM_Calc)  # Calculate the enthalpy of liquid phase mixture %unit: J/mol

        def ColumnAllTraysHV_Calc(m, tray):
            return sum([m.ColumnVapLvMFrac[tray, c] * m.ColumnAllTraysH0[tray, c] for c in m.Component]) + \
                   m.PR_RGas * (m.ColumnTrayTemp[tray] + 273.15) * (
                               (m.ColumnAllTraysZV[tray] - 1) - 1 / (2 ** 1.5) * (1 / m.ColumnAllTraysBV[tray]) *
                               m.ColumnAllTraysZVdivBV[tray] * \
                               sum([sum([m.ColumnVapLvMFrac[tray, c1] * m.ColumnVapLvMFrac[tray, c2] *
                                         m.ColumnAllTraysAV[tray] * (1 + m.ColumnAllTraysM[tray, c2, c1] + \
                                                                     m.ColumnAllTraysM[tray, c1, c2]) for c1 in
                                         m.Component]) for c2 in m.Component]))

        hpc_model.ColumnVapLvMEtlp = Expression(hpc_model.ColumnTrays, rule=ColumnAllTraysHV_Calc)

        # Calculate the enthalpy of liquid phase mixture %unit: J/mol
        def ColumnAllTraysHL_Calc(m, tray):
            return sum([m.ColumnLiqLvMFrac[tray, c] * m.ColumnAllTraysH0[tray, c] for c in m.Component]) + \
                   m.PR_RGas * (m.ColumnTrayTemp[tray] + 273.15) * (
                               (m.ColumnAllTraysZL[tray] - 1) - 1 / (2 ** 1.5) * (1 / m.ColumnAllTraysBL[tray]) *
                               m.ColumnAllTraysZLdivBL[tray] * \
                               sum([sum([m.ColumnLiqLvMFrac[tray, c1] * m.ColumnLiqLvMFrac[tray, c2] *
                                         m.ColumnAllTraysAL[tray] * (1 + m.ColumnAllTraysM[tray, c2, c1] + \
                                                                     m.ColumnAllTraysM[tray, c1, c2]) for c1 in
                                         m.Component]) for c2 in m.Component]))

        hpc_model.ColumnLiqLvMEtlp = Expression(hpc_model.ColumnTrays, rule=ColumnAllTraysHL_Calc)

        # --------Condensor Enthalpy---------
        def CondensorOutletTr_Calc(m, comp):
            return (m.ColumnTrayTemp[60] + 273.15) / m.PR_CrtcTemp[comp]

        hpc_model.CondensorOutletTr = Expression(hpc_model.Component, rule=CondensorOutletTr_Calc)

        def CondensorOutletalpha_Calc(m, comp):
            return (1 + m.PR_m[comp] * (1 - m.CondensorOutletSqrtTr[comp])) ** 2

        hpc_model.CondensorOutletalpha = Expression(hpc_model.Component, rule=CondensorOutletalpha_Calc)

        def CondensorOutletai_Calc(m, comp):
            return 0.45724 * m.CondensorOutletalpha[comp] * (m.PR_RGas ** 2) * (m.PR_CrtcTemp[comp] ** 2) / \
                   m.PR_CrtcPressure[comp]

        hpc_model.CondensorOutletai = Expression(hpc_model.Component, rule=CondensorOutletai_Calc)

        def CondensorOutletaij_Calc(m, comp1, comp2):
            return m.CondensorOutletSqrt_ai[comp1, comp2] * (1 - m.PR_Kij[comp1, comp2])

        hpc_model.CondensorOutletaij = Expression(hpc_model.Component, hpc_model.Component,
                                                   rule=CondensorOutletaij_Calc)

        def CondensorOutletaL_Calc(m):
            return sum([sum(
                [m.CondensorOutletaij[c1, c2] * m.ColumnVapLvMFrac[60, c1] * m.ColumnVapLvMFrac[60, c2] for c1 in
                 m.Component]) for c2 in m.Component])

        hpc_model.CondensorOutletaL = Expression(rule=CondensorOutletaL_Calc)

        def CondensorOutletbL_Calc(m):
            return sum([m.PR_bi[c] * m.ColumnVapLvMFrac[60, c] for c in m.Component])

        hpc_model.CondensorOutletbL = Expression(rule=CondensorOutletbL_Calc)

        def CondensorOutletAL_Calc(m):
            return m.CondensorOutletaL * m.CondensorPressure * 1000 / (
                        (m.PR_RGas * (m.ColumnTrayTemp[60] + 273.15)) ** 2)

        hpc_model.CondensorOutletAL = Expression(rule=CondensorOutletAL_Calc)

        def CondensorOutletBL_Calc(m):
            return m.CondensorOutletbL * m.CondensorPressure * 1000 / (m.PR_RGas * (m.ColumnTrayTemp[60] + 273.15))

        hpc_model.CondensorOutletBL = Expression(rule=CondensorOutletBL_Calc)

        def CondensorOutleta1_Calc(m):
            return m.CondensorOutletBL - 1

        hpc_model.CondensorOutleta1 = Expression(rule=CondensorOutleta1_Calc)

        def CondensorOutleta2_Calc(m):
            return m.CondensorOutletAL - 3 * m.CondensorOutletBL ** 2 - 2 * m.CondensorOutletBL

        hpc_model.CondensorOutleta2 = Expression(rule=CondensorOutleta2_Calc)

        def CondensorOutleta3_Calc(m):
            return m.CondensorOutletBL ** 2 + m.CondensorOutletBL ** 3 - m.CondensorOutletAL * m.CondensorOutletBL

        hpc_model.CondensorOutleta3 = Expression(rule=CondensorOutleta3_Calc)

        def CondensorOutletH0_Calc(m, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.ColumnTrayTemp[60] + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (m.PR_DHFORM[comp] + m.PR_CPIGDP1[comp] * (
                            (m.ColumnTrayTemp[60] + 273.15) - m.PR_RefTemp[comp]) + 2 * m.PR_CPIGDP2[comp] * \
                               m.PR_CPIGDP3[comp] * (-exp(2 / (m.ColumnTrayTemp[60] + 273.15) * m.PR_CPIGDP3[comp]) + \
                                                     exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp])) / (
                                           exp(2 / (m.ColumnTrayTemp[60] + 273.15) * m.PR_CPIGDP3[comp]) - 1) / (
                                           exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp]) - 1) - \
                               2 * m.PR_CPIGDP4[comp] * m.PR_CPIGDP5[comp] * (
                                           exp(2 / (m.ColumnTrayTemp[60] + 273.15) * m.PR_CPIGDP5[comp]) - exp(
                                       2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp])) / \
                               (exp(2 / (m.ColumnTrayTemp[60] + 273.15) * m.PR_CPIGDP5[comp]) + 1) / (
                                           exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp]) + 1))

        hpc_model.CondensorOutletH0 = Expression(hpc_model.Component, rule=CondensorOutletH0_Calc)

        def CondensorOutletM_Calc(m, comp1, comp2):
            return m.PR_m[comp1] * m.CondensorOutletSqrtTr[comp1] / (2 * sqrt(m.CondensorOutletalpha[comp1]))

        hpc_model.CondensorOutletM = Expression(hpc_model.Component, hpc_model.Component,
                                                 rule=CondensorOutletM_Calc)  # Calculate the enthalpy of liquid phase mixture %unit: J/mol

        def CondensorOutletHL_Calc(m):
            return sum([m.ColumnVapLvMFrac[60, c] * m.CondensorOutletH0[c] for c in m.Component]) + \
                   m.PR_RGas * (m.ColumnTrayTemp[60] + 273.15) * ((m.CondensorOutletZL - 1) - 1 / (2 ** 1.5) * (
                        1 / m.CondensorOutletBL) * m.CondensorOutletZLdivBL * \
                                                                  sum([sum([m.ColumnVapLvMFrac[60, c1] *
                                                                            m.ColumnVapLvMFrac[
                                                                                60, c2] * m.CondensorOutletAL * (
                                                                                        1 + m.CondensorOutletM[c2, c1] + \
                                                                                        m.CondensorOutletM[c1, c2]) for
                                                                            c1 in m.Component]) for c2 in m.Component]))

        hpc_model.CondensorOutMEtlp = Expression(rule=CondensorOutletHL_Calc)

        # --------Reboiler VLE & Enthalpy---------
        def ReboilerHdlpTr_Calc(m, comp):
            return (m.ReboilerTemp + 273.15) / m.PR_CrtcTemp[comp]

        hpc_model.ReboilerHdlpTr = Expression(hpc_model.Component, rule=ReboilerHdlpTr_Calc)

        def ReboilerHdlpalpha_Calc(m, comp):
            return (1 + m.PR_m[comp] * (1 - m.ReboilerHdlpSqrtTr[comp])) ** 2

        hpc_model.ReboilerHdlpalpha = Expression(hpc_model.Component, rule=ReboilerHdlpalpha_Calc)

        def ReboilerHdlpai_Calc(m, comp):
            return 0.45724 * m.ReboilerHdlpalpha[comp] * (m.PR_RGas ** 2) * (m.PR_CrtcTemp[comp] ** 2) / \
                   m.PR_CrtcPressure[comp]

        hpc_model.ReboilerHdlpai = Expression(hpc_model.Component, rule=ReboilerHdlpai_Calc)

        def ReboilerHdlpaij_Calc(m, comp1, comp2):
            return m.ReboilerHdlpSqrt_ai[comp1, comp2] * (1 - m.PR_Kij[comp1, comp2])

        hpc_model.ReboilerHdlpaij = Expression(hpc_model.Component, hpc_model.Component, rule=ReboilerHdlpaij_Calc)

        def ReboilerHdlpaL_Calc(m):
            return sum([sum(
                [m.ReboilerHdlpaij[c1, c2] * m.ReboilerLiqLvMFrac[c1] * m.ReboilerLiqLvMFrac[c2] for c1 in m.Component])
                        for c2 in m.Component])

        hpc_model.ReboilerHdlpaL = Expression(rule=ReboilerHdlpaL_Calc)

        def ReboilerHdlpbL_Calc(m):
            return sum([m.PR_bi[c] * m.ReboilerLiqLvMFrac[c] for c in m.Component])

        hpc_model.ReboilerHdlpbL = Expression(rule=ReboilerHdlpbL_Calc)

        def ReboilerHdlpAL_Calc(m):
            return m.ReboilerHdlpaL * m.ColumnTrayPressure[1] * 1000 / ((m.PR_RGas * (m.ReboilerTemp + 273.15)) ** 2)

        hpc_model.ReboilerHdlpAL = Expression(rule=ReboilerHdlpAL_Calc)

        def ReboilerHdlpBL_Calc(m):
            return m.ReboilerHdlpbL * m.ColumnTrayPressure[1] * 1000 / (m.PR_RGas * (m.ReboilerTemp + 273.15))

        hpc_model.ReboilerHdlpBL = Expression(rule=ReboilerHdlpBL_Calc)

        def ReboilerHdlpa1_Calc(m):
            return m.ReboilerHdlpBL - 1

        hpc_model.ReboilerHdlpa1 = Expression(rule=ReboilerHdlpa1_Calc)

        def ReboilerHdlpa2_Calc(m):
            return m.ReboilerHdlpAL - 3 * m.ReboilerHdlpBL ** 2 - 2 * m.ReboilerHdlpBL

        hpc_model.ReboilerHdlpa2 = Expression(rule=ReboilerHdlpa2_Calc)

        def ReboilerHdlpa3_Calc(m):
            return m.ReboilerHdlpBL ** 2 + m.ReboilerHdlpBL ** 3 - m.ReboilerHdlpAL * m.ReboilerHdlpBL

        hpc_model.ReboilerHdlpa3 = Expression(rule=ReboilerHdlpa3_Calc)

        def ReboilerHdlpaV_Calc(m):
            return sum([sum(
                [m.ReboilerHdlpaij[c1, c2] * m.ReboilerVapLvMFrac[c1] * m.ReboilerVapLvMFrac[c2] for c1 in m.Component])
                        for c2 in m.Component])

        hpc_model.ReboilerHdlpaV = Expression(rule=ReboilerHdlpaV_Calc)

        def ReboilerHdlpbV_Calc(m):
            return sum([m.PR_bi[c] * m.ReboilerVapLvMFrac[c] for c in m.Component])

        hpc_model.ReboilerHdlpbV = Expression(rule=ReboilerHdlpbV_Calc)

        def ReboilerHdlpAV_Calc(m):
            return m.ReboilerHdlpaV * m.ColumnTrayPressure[1] * 1000 / ((m.PR_RGas * (m.ReboilerTemp + 273.15)) ** 2)

        hpc_model.ReboilerHdlpAV = Expression(rule=ReboilerHdlpAV_Calc)

        def ReboilerHdlpBV_Calc(m):
            return m.ReboilerHdlpbV * m.ColumnTrayPressure[1] * 1000 / (m.PR_RGas * (m.ReboilerTemp + 273.15))

        hpc_model.ReboilerHdlpBV = Expression(rule=ReboilerHdlpBV_Calc)

        def ReboilerHdlpb1_Calc(m):
            return m.ReboilerHdlpBV - 1

        hpc_model.ReboilerHdlpb1 = Expression(rule=ReboilerHdlpb1_Calc)

        def ReboilerHdlpb2_Calc(m):
            return m.ReboilerHdlpAV - 3 * m.ReboilerHdlpBV ** 2 - 2 * m.ReboilerHdlpBV

        hpc_model.ReboilerHdlpb2 = Expression(rule=ReboilerHdlpb2_Calc)

        def ReboilerHdlpb3_Calc(m):
            return m.ReboilerHdlpBV ** 2 + m.ReboilerHdlpBV ** 3 - m.ReboilerHdlpAV * m.ReboilerHdlpBV

        hpc_model.ReboilerHdlpb3 = Expression(rule=ReboilerHdlpb3_Calc)

        def ReboilerHdlpSL_Calc(m, comp):
            return sum([m.ReboilerHdlpaij[c, comp] * m.ReboilerLiqLvMFrac[c] for c in m.Component])

        hpc_model.ReboilerHdlpSL = Expression(hpc_model.Component, rule=ReboilerHdlpSL_Calc)

        def ReboilerHdlpSV_Calc(m, comp):
            return sum([m.ReboilerHdlpaij[c, comp] * m.ReboilerVapLvMFrac[c] for c in m.Component])

        hpc_model.ReboilerHdlpSV = Expression(hpc_model.Component, rule=ReboilerHdlpSV_Calc)

        def ReboilerHdlpPhiL_Calc(m, comp):
            return exp(m.PR_bi[comp] / m.ReboilerHdlpbL * (m.ReboilerHdlpZL - 1) - \
                       m.ReboilerHdlpZLBL - m.ReboilerHdlpAL / m.ReboilerHdlpBL / 2 / sqrt(2) * (2 * \
                                                                                                 m.ReboilerHdlpSL[
                                                                                                     comp] / m.ReboilerHdlpaL -
                                                                                                 m.PR_bi[
                                                                                                     comp] / m.ReboilerHdlpbL) * \
                       m.ReboilerHdlpZLdivBL)

        hpc_model.ReboilerHdlpPhiL = Expression(hpc_model.Component, rule=ReboilerHdlpPhiL_Calc)

        def ReboilerHdlpPhiV_Calc(m, comp):
            return exp(m.PR_bi[comp] / m.ReboilerHdlpbV * (m.ReboilerHdlpZV - 1) - \
                       m.ReboilerHdlpZVBV - m.ReboilerHdlpAV / m.ReboilerHdlpBV / 2 / sqrt(2) * (2 * \
                                                                                                 m.ReboilerHdlpSV[
                                                                                                     comp] / m.ReboilerHdlpaV -
                                                                                                 m.PR_bi[
                                                                                                     comp] / m.ReboilerHdlpbV) * \
                       m.ReboilerHdlpZVdivBV)

        hpc_model.ReboilerHdlpPhiV = Expression(hpc_model.Component, rule=ReboilerHdlpPhiV_Calc)

        def ReboilerHdlpH0_Calc(m, comp):  # unit: J/mol
            if comp == 'Argon':
                return 20786 * ((m.ReboilerTemp + 273.15) - 298.15) / 1000
            else:
                return 4.18 * (m.PR_DHFORM[comp] + m.PR_CPIGDP1[comp] * (
                            (m.ReboilerTemp + 273.15) - m.PR_RefTemp[comp]) + 2 * m.PR_CPIGDP2[comp] * \
                               m.PR_CPIGDP3[comp] * (-exp(2 / (m.ReboilerTemp + 273.15) * m.PR_CPIGDP3[comp]) + \
                                                     exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp])) / (
                                           exp(2 / (m.ReboilerTemp + 273.15) * m.PR_CPIGDP3[comp]) - 1) / (
                                           exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP3[comp]) - 1) - \
                               2 * m.PR_CPIGDP4[comp] * m.PR_CPIGDP5[comp] * (
                                           exp(2 / (m.ReboilerTemp + 273.15) * m.PR_CPIGDP5[comp]) - exp(
                                       2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp])) / \
                               (exp(2 / (m.ReboilerTemp + 273.15) * m.PR_CPIGDP5[comp]) + 1) / (
                                           exp(2 / m.PR_RefTemp[comp] * m.PR_CPIGDP5[comp]) + 1))

        hpc_model.ReboilerHdlpH0 = Expression(hpc_model.Component, rule=ReboilerHdlpH0_Calc)

        def ReboilerHdlpM_Calc(m, comp1, comp2):
            return m.PR_m[comp1] * m.ReboilerHdlpSqrtTr[comp1] / (2 * sqrt(m.ReboilerHdlpalpha[comp1]))

        hpc_model.ReboilerHdlpM = Expression(hpc_model.Component, hpc_model.Component,
                                              rule=ReboilerHdlpM_Calc)  # Calculate the enthalpy of liquid phase mixture %unit: J/mol

        def ReboilerHdlpHV_Calc(m):
            return sum([m.ReboilerVapLvMFrac[c] * m.ReboilerHdlpH0[c] for c in m.Component]) + \
                   m.PR_RGas * (m.ReboilerTemp + 273.15) * ((m.ReboilerHdlpZV - 1) - 1 / (2 ** 1.5) * (
                        1 / m.ReboilerHdlpBV) * m.ReboilerHdlpZVdivBV * \
                                                            sum([sum([m.ReboilerVapLvMFrac[c1] * m.ReboilerVapLvMFrac[
                                                                c2] * m.ReboilerHdlpAV * (1 + m.ReboilerHdlpM[c2, c1] + \
                                                                                          m.ReboilerHdlpM[c1, c2]) for
                                                                      c1 in m.Component]) for c2 in m.Component]))

        hpc_model.ReboilerVapLvMEtlp = Expression(rule=ReboilerHdlpHV_Calc)

        # Calculate the enthalpy of liquid phase mixture %unit: J/mol
        def ReboilerHdlpHL_Calc(m):
            return sum([m.ReboilerLiqLvMFrac[c] * m.ReboilerHdlpH0[c] for c in m.Component]) + \
                   m.PR_RGas * (m.ReboilerTemp + 273.15) * ((m.ReboilerHdlpZL - 1) - 1 / (2 ** 1.5) * (
                        1 / m.ReboilerHdlpBL) * m.ReboilerHdlpZLdivBL * \
                                                            sum([sum([m.ReboilerLiqLvMFrac[c1] * m.ReboilerLiqLvMFrac[
                                                                c2] * m.ReboilerHdlpAL * (1 + m.ReboilerHdlpM[c2, c1] + \
                                                                                          m.ReboilerHdlpM[c1, c2]) for
                                                                      c1 in m.Component]) for c2 in m.Component]))

        hpc_model.ReboilerLiqLvMEtlp = Expression(rule=ReboilerHdlpHL_Calc)

        # ===================================
        #
        #         Constraints
        #
        # ===================================
        def NitrogenOxygenKijEqu(m):
            yield hpc_model.PR_Kij['Nitrogen', 'Oxygen'] == hpc_model.NitrogenOxygenKij
            yield hpc_model.PR_Kij['Oxygen', 'Nitrogen'] == hpc_model.NitrogenOxygenKij

        hpc_model.NitrogenOxygenKijEqu = ConstraintList(rule=NitrogenOxygenKijEqu)

        # ----------------------------------
        #           Feed
        # ----------------------------------

        # ----------------------------------
        #           TaCooler
        # ----------------------------------
        # --------TaCooler Mass Balance---------
        def TaCoolerMassBlnc(m, comp):
            return m.TaCoolerVapMFlow * m.TaCoolerVapMFrac[comp] + \
                   m.TaCoolerLiqMFlow * m.TaCoolerLiqMFrac[comp] - m.FeedTaMFlow * m.FeedTaMFrac[comp] == 0

        hpc_model.TaCoolerMassBlnc = Constraint(hpc_model.Component, rule=TaCoolerMassBlnc)

        def TaCoolerVF_Spec(m):
            return (m.TaCoolerVF * m.FeedTaMFlow - m.TaCoolerVapMFlow) * 1e-2 == 0

        hpc_model.TaCoolerVF_Spec = Constraint(rule=TaCoolerVF_Spec)

        # --------TaCooler Energy Balance---------
        def TaCoolerEngBlnc(m):
            return ((m.FeedTaMEtlp - m.TaCoolerMEtlp) * m.FeedTaMFlow - m.TaCoolerMHeat) * 1.000000e-05 == 0

        hpc_model.TaCoolerEngBlnc = Constraint(rule=TaCoolerEngBlnc)

        # --------TaCooler Phase Equilibrium---------
        def TaCoolerSelfSqrtTr_Calc(m, comp):
            return m.TaCoolerSelfSqrtTr[comp] ** 2 - m.TaCoolerSelfTr[comp] == 0

        hpc_model.TaCoolerSelfSqrtTr_Calc = Constraint(hpc_model.Component, rule=TaCoolerSelfSqrtTr_Calc)

        def TaCoolerSelfSqrt_ai_Calc(m, comp1, comp2):
            return m.TaCoolerSelfSqrt_ai[comp1, comp2] ** 2 - m.TaCoolerSelfai[comp1] * m.TaCoolerSelfai[comp2] == 0

        hpc_model.TaCoolerSelfSqrt_ai_Calc = Constraint(hpc_model.Component, hpc_model.Component,
                                                         rule=TaCoolerSelfSqrt_ai_Calc)

        def TaCoolerSelfZL_Calc(m):
            return m.TaCoolerSelfZL ** 3 + m.TaCoolerSelfa1 * m.TaCoolerSelfZL ** 2 + \
                   m.TaCoolerSelfa2 * m.TaCoolerSelfZL + m.TaCoolerSelfa3 == 0

        hpc_model.TaCoolerSelfZL_Calc = Constraint(rule=TaCoolerSelfZL_Calc)

        def TaCoolerSelfZL_Con1(m):
            return 3 * m.TaCoolerSelfZL ** 2 + 2 * m.TaCoolerSelfa1 * m.TaCoolerSelfZL + m.TaCoolerSelfa2 >= 0

        hpc_model.TaCoolerSelfZL_Con1 = Constraint(rule=TaCoolerSelfZL_Con1)

        def TaCoolerSelfZL_Con2(m):
            return 6 * m.TaCoolerSelfZL + 2 * m.TaCoolerSelfa1 <= 0

        hpc_model.TaCoolerSelfZL_Con2 = Constraint(rule=TaCoolerSelfZL_Con2)

        def TaCoolerSelfZV_Calc(m):
            return m.TaCoolerSelfZV ** 3 + m.TaCoolerSelfb1 * m.TaCoolerSelfZV ** 2 + \
                   m.TaCoolerSelfb2 * m.TaCoolerSelfZV + m.TaCoolerSelfb3 == 0

        hpc_model.TaCoolerSelfZV_Calc = Constraint(rule=TaCoolerSelfZV_Calc)

        def TaCoolerSelfZV_Con1(m):
            return 3 * m.TaCoolerSelfZV ** 2 + 2 * m.TaCoolerSelfb1 * m.TaCoolerSelfZV + m.TaCoolerSelfb2 >= 0

        hpc_model.TaCoolerSelfZV_Con1 = Constraint(rule=TaCoolerSelfZV_Con1)

        def TaCoolerSelfZV_Con2(m):
            return 6 * m.TaCoolerSelfZV + 2 * m.TaCoolerSelfb1 >= 0

        hpc_model.TaCoolerSelfZV_Con2 = Constraint(rule=TaCoolerSelfZV_Con2)

        def TaCoolerSelfZLBL_Calc(m):
            return exp(m.TaCoolerSelfZLBL) - (m.TaCoolerSelfZL - m.TaCoolerSelfBL) == 0

        hpc_model.TaCoolerSelfZLBL_Calc = Constraint(rule=TaCoolerSelfZLBL_Calc)

        def TaCoolerSelfZLdivBL_Calc(m):
            return exp(m.TaCoolerSelfZLdivBL) - (m.TaCoolerSelfZL + 2.414 * m.TaCoolerSelfBL) / (
                        m.TaCoolerSelfZL - 0.414 * m.TaCoolerSelfBL) == 0

        hpc_model.TaCoolerSelfZLdivBL_Calc = Constraint(rule=TaCoolerSelfZLdivBL_Calc)

        def TaCoolerSelfZVBV_Calc(m):
            return exp(m.TaCoolerSelfZVBV) - (m.TaCoolerSelfZV - m.TaCoolerSelfBV) == 0

        hpc_model.TaCoolerSelfZVBV_Calc = Constraint(rule=TaCoolerSelfZVBV_Calc)

        def TaCoolerSelfZVdivBV_Calc(m):
            return exp(m.TaCoolerSelfZVdivBV) - (m.TaCoolerSelfZV + 2.414 * m.TaCoolerSelfBV) / (
                        m.TaCoolerSelfZV - 0.414 * m.TaCoolerSelfBV) == 0

        hpc_model.TaCoolerSelfZVdivBV_Calc = Constraint(rule=TaCoolerSelfZVdivBV_Calc)

        def TaCoolerSelfK_Calc(m, comp):
            return m.TaCoolerSelfPhiL[comp] / m.TaCoolerSelfPhiV[comp] - m.TaCoolerSelfVLE_K[comp] == 0

        hpc_model.TaCoolerSelfK_Calc = Constraint(hpc_model.Component, rule=TaCoolerSelfK_Calc)

        def TaCoolerSelfVLEEqu(m, comp):
            return m.TaCoolerLiqMFrac[comp] * m.TaCoolerSelfVLE_K[comp] - m.TaCoolerVapMFrac[comp] == 0

        hpc_model.TaCoolerSelfVLEEqu = Constraint(hpc_model.Component, rule=TaCoolerSelfVLEEqu)

        # --------TaCooler Summation---------
        def TaCoolerLiqSum(m):
            return sum([m.TaCoolerLiqMFrac[c] for c in m.Component]) == 1

        hpc_model.TaCoolerLiqSum = Constraint(rule=TaCoolerLiqSum)

        def TaCoolerVapSum(m):
            return sum([m.TaCoolerVapMFrac[c] for c in m.Component]) == 1

        hpc_model.TaCoolerVapSum = Constraint(rule=TaCoolerVapSum)

        # ----------------------------------
        #           MaCooler
        # ----------------------------------
        # --------MaCooler Enthalpy---------
        def MaCoolerSelfSqrtTr_Calc(m, comp):
            return m.MaCoolerSelfSqrtTr[comp] ** 2 - m.MaCoolerSelfTr[comp] == 0

        hpc_model.MaCoolerSelfSqrtTr_Calc = Constraint(hpc_model.Component, rule=MaCoolerSelfSqrtTr_Calc)

        def MaCoolerSelfSqrt_ai_Calc(m, comp1, comp2):
            return m.MaCoolerSelfSqrt_ai[comp1, comp2] ** 2 - m.MaCoolerSelfai[comp1] * m.MaCoolerSelfai[comp2] == 0

        hpc_model.MaCoolerSelfSqrt_ai_Calc = Constraint(hpc_model.Component, hpc_model.Component,
                                                         rule=MaCoolerSelfSqrt_ai_Calc)

        def MaCoolerSelfZV_Calc(m):
            return m.MaCoolerSelfZV ** 3 + m.MaCoolerSelfb1 * m.MaCoolerSelfZV ** 2 + \
                   m.MaCoolerSelfb2 * m.MaCoolerSelfZV + m.MaCoolerSelfb3 == 0

        hpc_model.MaCoolerSelfZV_Calc = Constraint(rule=MaCoolerSelfZV_Calc)

        def MaCoolerSelfZV_Con1(m):
            return 3 * m.MaCoolerSelfZV ** 2 + 2 * m.MaCoolerSelfb1 * m.MaCoolerSelfZV + m.MaCoolerSelfb2 >= 0

        hpc_model.MaCoolerSelfZV_Con1 = Constraint(rule=MaCoolerSelfZV_Con1)

        def MaCoolerSelfZV_Con2(m):
            return 6 * m.MaCoolerSelfZV + 2 * m.MaCoolerSelfb1 >= 0

        hpc_model.MaCoolerSelfZV_Con2 = Constraint(rule=MaCoolerSelfZV_Con2)

        def MaCoolerSelfZVdivBV_Calc(m):
            return exp(m.MaCoolerSelfZVdivBV) - (m.MaCoolerSelfZV + 2.414 * m.MaCoolerSelfBV) / (
                        m.MaCoolerSelfZV - 0.414 * m.MaCoolerSelfBV) == 0

        hpc_model.MaCoolerSelfZVdivBV_Calc = Constraint(rule=MaCoolerSelfZVdivBV_Calc)

        # --------MaCooler Energy Balance---------
        def MaCoolerEngBlnc(m):
            return ((m.FeedMaMEtlp - m.MaCoolerMEtlp) * m.FeedMaMFlow - m.MaCoolerMHeat) * 1.000000e-04 == 0

        hpc_model.MaCoolerEngBlnc = Constraint(rule=MaCoolerEngBlnc)

        # ----------------------------------
        #           AirMixer
        # ----------------------------------
        # --------AirMixer Mass Balance---------
        def AirMixerOutMassBlnc(m):
            return (m.AirMixerOutMFlow - m.FeedMaMFlow - m.FeedTaMFlow) * 0.001000 == 0

        hpc_model.AirMixerOutMassBlnc = Constraint(rule=AirMixerOutMassBlnc)

        def AirMixerOutCompMassBlnc(m, comp):
            return (m.AirMixerOutMFrac[comp] * m.AirMixerOutMFlow - m.FeedMaMFrac[comp] * m.FeedMaMFlow -
                    m.FeedTaMFrac[comp] * m.FeedTaMFlow) * 0.001000 == 0

        hpc_model.AirMixerOutCompMassBlnc = Constraint(hpc_model.Component, rule=AirMixerOutCompMassBlnc)

        # --------AirMixer Energy Balance---------
        def AirMixerEnergyBlnc(m):
            return (
                               m.AirMixerOutMEtlp * m.AirMixerOutMFlow - m.MaCoolerMEtlp * m.FeedMaMFlow - m.TaCoolerMEtlp * m.FeedTaMFlow) * 1.000000e-08 == 0

        hpc_model.AirMixerEnergyBlnc = Constraint(rule=AirMixerEnergyBlnc)

        # ----------------------------------
        #           Column
        # ----------------------------------
        # --------Column Mass Balance---------
        def ColumnMassBlnc(m, tray, comp):
            if tray == 1:
                return (m.HPCZeroReboiledNull * m.HPCZeroReboiledNull + m.ColumnLiqLvMCompFlow[tray + 1, comp] -
                        m.ColumnLiqLvMCompFlow[tray, comp] - m.ColumnVapLvMCompFlow[tray, comp] + m.AirMixerOutMFlow *
                        m.AirMixerOutMFrac[comp]) * 0.001000 == 0
            elif tray == 60:
                return (m.CondensorRefMFlow * m.ColumnVapLvMFrac[60, comp] + m.ColumnVapLvMCompFlow[tray - 1, comp] -
                        m.ColumnLiqLvMCompFlow[tray, comp] - m.ColumnVapLvMCompFlow[tray, comp]) * 0.001000 == 0
            else:
                return (m.ColumnLiqLvMCompFlow[tray + 1, comp] + m.ColumnVapLvMCompFlow[tray - 1, comp] -
                        m.ColumnLiqLvMCompFlow[tray, comp] - m.ColumnVapLvMCompFlow[tray, comp]) * 0.001000 == 0

        hpc_model.ColumnMassBlnc = Constraint(hpc_model.ColumnTrays, hpc_model.Component, rule=ColumnMassBlnc)

        # --------Column Energy Balance---------
        def ColumnEngBlnc(m, tray):
            if tray == 1:
                return (m.HPCZeroReboiledNull * m.HPCZeroReboiledNull \
                        + m.ColumnLiqLvMFlow[tray + 1] * m.ColumnLiqLvMEtlp[tray + 1] \
                        - m.ColumnLiqLvMFlow[tray] * m.ColumnLiqLvMEtlp[tray] \
                        - m.ColumnVapLvMFlow[tray] * m.ColumnVapLvMEtlp[tray] \
                        + m.AirMixerOutMFlow * m.AirMixerOutMEtlp) * 0.000010 == 0
            elif tray == 60:
                return (m.CondensorRefMFlow * m.CondensorOutMEtlp \
                        + m.ColumnVapLvMFlow[tray - 1] * m.ColumnVapLvMEtlp[tray - 1] \
                        - m.ColumnLiqLvMFlow[tray] * m.ColumnLiqLvMEtlp[tray] \
                        - m.ColumnVapLvMFlow[tray] * m.ColumnVapLvMEtlp[tray]) * 0.000010 == 0
            else:
                return (m.ColumnLiqLvMFlow[tray + 1] * m.ColumnLiqLvMEtlp[tray + 1] \
                        + m.ColumnVapLvMFlow[tray - 1] * m.ColumnVapLvMEtlp[tray - 1] \
                        - m.ColumnLiqLvMFlow[tray] * m.ColumnLiqLvMEtlp[tray] \
                        - m.ColumnVapLvMFlow[tray] * m.ColumnVapLvMEtlp[tray] \
                - m.ColumnUA * (m.ColumnTrayTemp[tray]-25))*0.000010 == 0

        hpc_model.ColumnEngBlnc = Constraint(hpc_model.ColumnTrays, rule=ColumnEngBlnc)

        # --------Column Phase Equilibrium & System Parts---------
        def ColumnAllTraysSqrtTr_Calc(m, tray, comp):
            return m.ColumnAllTraysSqrtTr[tray, comp] ** 2 - m.ColumnAllTraysTr[tray, comp] == 0

        hpc_model.ColumnAllTraysSqrtTr_Calc = Constraint(hpc_model.ColumnTrays, hpc_model.Component,
                                                          rule=ColumnAllTraysSqrtTr_Calc)

        def ColumnAllTraysSqrt_ai_Calc(m, tray, comp1, comp2):
            return m.ColumnAllTraysSqrt_ai[tray, comp1, comp2] ** 2 - m.ColumnAllTraysai[tray, comp1] * \
                   m.ColumnAllTraysai[tray, comp2] == 0

        hpc_model.ColumnAllTraysSqrt_ai_Calc = Constraint(hpc_model.ColumnTrays, hpc_model.Component,
                                                           hpc_model.Component, rule=ColumnAllTraysSqrt_ai_Calc)

        def ColumnAllTraysZL_Calc(m, tray):
            return m.ColumnAllTraysZL[tray] ** 3 + m.ColumnAllTraysa1[tray] * m.ColumnAllTraysZL[tray] ** 2 + \
                   m.ColumnAllTraysa2[tray] * m.ColumnAllTraysZL[tray] + m.ColumnAllTraysa3[tray] == 0

        hpc_model.ColumnAllTraysZL_Calc = Constraint(hpc_model.ColumnTrays, rule=ColumnAllTraysZL_Calc)

        def ColumnAllTraysZL_Con1(m, tray):
            return 3 * m.ColumnAllTraysZL[tray] ** 2 + 2 * m.ColumnAllTraysa1[tray] * m.ColumnAllTraysZL[tray] + \
                   m.ColumnAllTraysa2[tray] >= 0

        hpc_model.ColumnAllTraysZL_Con1 = Constraint(hpc_model.ColumnTrays, rule=ColumnAllTraysZL_Con1)

        def ColumnAllTraysZL_Con2(m, tray):
            return 6 * m.ColumnAllTraysZL[tray] + 2 * m.ColumnAllTraysa1[tray] <= 0

        hpc_model.ColumnAllTraysZL_Con2 = Constraint(hpc_model.ColumnTrays, rule=ColumnAllTraysZL_Con2)

        def ColumnAllTraysZV_Calc(m, tray):
            return m.ColumnAllTraysZV[tray] ** 3 + m.ColumnAllTraysb1[tray] * m.ColumnAllTraysZV[tray] ** 2 + \
                   m.ColumnAllTraysb2[tray] * m.ColumnAllTraysZV[tray] + m.ColumnAllTraysb3[tray] == 0

        hpc_model.ColumnAllTraysZV_Calc = Constraint(hpc_model.ColumnTrays, rule=ColumnAllTraysZV_Calc)

        def ColumnAllTraysZV_Con1(m, tray):
            return 3 * m.ColumnAllTraysZV[tray] ** 2 + 2 * m.ColumnAllTraysb1[tray] * m.ColumnAllTraysZV[tray] + \
                   m.ColumnAllTraysb2[tray] >= 0

        hpc_model.ColumnAllTraysZV_Con1 = Constraint(hpc_model.ColumnTrays, rule=ColumnAllTraysZV_Con1)

        def ColumnAllTraysZV_Con2(m, tray):
            return 6 * m.ColumnAllTraysZV[tray] + 2 * m.ColumnAllTraysb1[tray] >= 0

        hpc_model.ColumnAllTraysZV_Con2 = Constraint(hpc_model.ColumnTrays, rule=ColumnAllTraysZV_Con2)

        def ColumnAllTraysZLBL_Calc(m, tray):
            return exp(m.ColumnAllTraysZLBL[tray]) - (m.ColumnAllTraysZL[tray] - m.ColumnAllTraysBL[tray]) == 0

        hpc_model.ColumnAllTraysZLBL_Calc = Constraint(hpc_model.ColumnTrays, rule=ColumnAllTraysZLBL_Calc)

        def ColumnAllTraysZLdivBL_Calc(m, tray):
            return exp(m.ColumnAllTraysZLdivBL[tray]) - (
                        m.ColumnAllTraysZL[tray] + 2.414 * m.ColumnAllTraysBL[tray]) / (
                               m.ColumnAllTraysZL[tray] - 0.414 * m.ColumnAllTraysBL[tray]) == 0

        hpc_model.ColumnAllTraysZLdivBL_Calc = Constraint(hpc_model.ColumnTrays, rule=ColumnAllTraysZLdivBL_Calc)

        def ColumnAllTraysZVBV_Calc(m, tray):
            return exp(m.ColumnAllTraysZVBV[tray]) - (m.ColumnAllTraysZV[tray] - m.ColumnAllTraysBV[tray]) == 0

        hpc_model.ColumnAllTraysZVBV_Calc = Constraint(hpc_model.ColumnTrays, rule=ColumnAllTraysZVBV_Calc)

        def ColumnAllTraysZVdivBV_Calc(m, tray):
            return exp(m.ColumnAllTraysZVdivBV[tray]) - (
                        m.ColumnAllTraysZV[tray] + 2.414 * m.ColumnAllTraysBV[tray]) / (
                               m.ColumnAllTraysZV[tray] - 0.414 * m.ColumnAllTraysBV[tray]) == 0

        hpc_model.ColumnAllTraysZVdivBV_Calc = Constraint(hpc_model.ColumnTrays, rule=ColumnAllTraysZVdivBV_Calc)

        def ColumnAllTraysK_Calc(m, tray, comp):
            return m.ColumnAllTraysPhiL[tray, comp] / m.ColumnAllTraysPhiV[tray, comp] - m.ColumnAllTraysVLE_K[
                tray, comp] == 0

        hpc_model.ColumnAllTraysK_Calc = Constraint(hpc_model.ColumnTrays, hpc_model.Component,
                                                     rule=ColumnAllTraysK_Calc)

        def ColumnAllTraysVLEEqu(m, tray, comp):
            return m.ColumnLiqLvMFrac[tray, comp] * m.ColumnAllTraysVLE_K[tray, comp] - m.ColumnVapEqilMFrac[
                tray, comp] == 0

        hpc_model.ColumnAllTraysVLEEqu = Constraint(hpc_model.ColumnTrays, hpc_model.Component,
                                                     rule=ColumnAllTraysVLEEqu)

        def ColumnTrayEffEqu(m, tray, comp):
            if tray == 1:
                return m.ColumnVapLvMFrac[tray, comp] == m.HPCZeroReboiledNull + 1 * (
                            m.ColumnVapEqilMFrac[tray, comp] - m.HPCZeroReboiledNull)
            else:
                return m.ColumnVapLvMFrac[tray, comp] == m.ColumnVapLvMFrac[tray - 1, comp] + m.ColumnEffTray * (
                        m.ColumnVapEqilMFrac[tray, comp] - m.ColumnVapLvMFrac[tray - 1, comp])

        hpc_model.ColumnTrayEffEqu = Constraint(hpc_model.ColumnTrays, hpc_model.Component, rule=ColumnTrayEffEqu)

        # --------Column Summation---------
        def ColumnLiqSum(m, tray):
            return sum([m.ColumnLiqLvMFrac[tray, c] for c in m.Component]) == 1

        hpc_model.ColumnLiqSum = Constraint(hpc_model.ColumnTrays, rule=ColumnLiqSum)

        def ColumnVapSum(m, tray):
            return sum([m.ColumnVapLvMFrac[tray, c] for c in m.Component]) == 1

        hpc_model.ColumnVapSum = Constraint(hpc_model.ColumnTrays, rule=ColumnVapSum)

        # --------Column Pressure Profile---------
        def ColumnPProf(m, tray):
            return ((m.ColumnTopPressure - m.ColumnBtmPressure) / 59 * (tray - 1) + m.ColumnBtmPressure -
                    m.ColumnTrayPressure[tray]) * 0.010000 == 0

        hpc_model.ColumnPProf = Constraint(hpc_model.ColumnTrays, rule=ColumnPProf)

        # ----------------------------------
        #           Sump
        # ----------------------------------
        # --------Sump Mass Balance---------
        def SumpMassBlnc(m):
            return (m.ColumnLiqLvMFlow[1] - m.SumpOutMFlow) * 0.001000 == 0

        hpc_model.SumpMassBlnc = Constraint(rule=SumpMassBlnc)

        def SumpHldpMFracSpec(m, comp):
            return m.SumpHldpMFrac[comp] - m.ColumnLiqLvMFrac[1, comp] == 0

        hpc_model.SumpHldpMFracSpec = Constraint(hpc_model.Component, rule=SumpHldpMFracSpec)

        # ----------------------------------
        #           Condensor
        # ----------------------------------
        # --------Condensor Mass Balance---------
        def CondensorMassBlnc(m):
            return (m.CondensorRefMFlow + m.CondensorPrdtMFlow - m.GNSplitterOutBMFlow) * 0.001000 == 0

        hpc_model.CondensorMassBlnc = Constraint(rule=CondensorMassBlnc)

        # --------Condensor Energy Balance---------
        def CondensorEngBlnc(m):
            return (m.GNSplitterOutBMFlow * (
                        m.ColumnVapLvMEtlp[60] - m.CondensorOutMEtlp) - m.CondensorMHeatOut) * 0.000010 == 0

        hpc_model.CondensorEngBlnc = Constraint(rule=CondensorEngBlnc)

        # --------Condensor Enthalpy---------
        def CondensorOutletSqrtTr_Calc(m, comp):
            return m.CondensorOutletSqrtTr[comp] ** 2 - m.CondensorOutletTr[comp] == 0

        hpc_model.CondensorOutletSqrtTr_Calc = Constraint(hpc_model.Component, rule=CondensorOutletSqrtTr_Calc)

        def CondensorOutletSqrt_ai_Calc(m, comp1, comp2):
            return m.CondensorOutletSqrt_ai[comp1, comp2] ** 2 - m.CondensorOutletai[comp1] * m.CondensorOutletai[
                comp2] == 0

        hpc_model.CondensorOutletSqrt_ai_Calc = Constraint(hpc_model.Component, hpc_model.Component,
                                                            rule=CondensorOutletSqrt_ai_Calc)

        def CondensorOutletZL_Calc(m):
            return m.CondensorOutletZL ** 3 + m.CondensorOutleta1 * m.CondensorOutletZL ** 2 + \
                   m.CondensorOutleta2 * m.CondensorOutletZL + m.CondensorOutleta3 == 0

        hpc_model.CondensorOutletZL_Calc = Constraint(rule=CondensorOutletZL_Calc)

        def CondensorOutletZL_Con1(m):
            return 3 * m.CondensorOutletZL ** 2 + 2 * m.CondensorOutleta1 * m.CondensorOutletZL + m.CondensorOutleta2 >= 0

        hpc_model.CondensorOutletZL_Con1 = Constraint(rule=CondensorOutletZL_Con1)

        def CondensorOutletZL_Con2(m):
            return 6 * m.CondensorOutletZL + 2 * m.CondensorOutleta1 <= 0

        hpc_model.CondensorOutletZL_Con2 = Constraint(rule=CondensorOutletZL_Con2)

        def CondensorOutletZLdivBL_Calc(m):
            return exp(m.CondensorOutletZLdivBL) - (m.CondensorOutletZL + 2.414 * m.CondensorOutletBL) / (
                        m.CondensorOutletZL - 0.414 * m.CondensorOutletBL) == 0

        hpc_model.CondensorOutletZLdivBL_Calc = Constraint(rule=CondensorOutletZLdivBL_Calc)

        # ----------------------------------
        #           Reboiler
        # ----------------------------------
        # --------Reboiler Mass Balance---------
        def ReboilerMassBlnc(m, comp):
            return (m.SumpOutMFlow * m.ColumnLiqLvMFrac[1, comp] - \
                    m.ReboilerLiqLvMFlow * m.ReboilerLiqLvMFrac[comp] - \
                    m.ReboilerVapLvMFlow * m.ReboilerVapLvMFrac[comp]) * 0.001000 == 0

        hpc_model.ReboilerMassBlnc = Constraint(hpc_model.Component, rule=ReboilerMassBlnc)

        # --------Reboiler Energy Balance---------
        def ReboilerEngBlnc(m):
            return (m.SumpOutMFlow * m.ColumnLiqLvMEtlp[1] + m.CondensorMHeatOut - \
                    m.ReboilerLiqLvMFlow * m.ReboilerLiqLvMEtlp - \
                    m.ReboilerVapLvMFlow * m.ReboilerVapLvMEtlp) * 0.000010 == 0

        hpc_model.ReboilerEngBlnc = Constraint(rule=ReboilerEngBlnc)

        # --------Reboiler Phase Equilibrium---------
        def ReboilerHdlpSqrtTr_Calc(m, comp):
            return m.ReboilerHdlpSqrtTr[comp] ** 2 - m.ReboilerHdlpTr[comp] == 0

        hpc_model.ReboilerHdlpSqrtTr_Calc = Constraint(hpc_model.Component, rule=ReboilerHdlpSqrtTr_Calc)

        def ReboilerHdlpSqrt_ai_Calc(m, comp1, comp2):
            return m.ReboilerHdlpSqrt_ai[comp1, comp2] ** 2 - m.ReboilerHdlpai[comp1] * m.ReboilerHdlpai[comp2] == 0

        hpc_model.ReboilerHdlpSqrt_ai_Calc = Constraint(hpc_model.Component, hpc_model.Component,
                                                         rule=ReboilerHdlpSqrt_ai_Calc)

        def ReboilerHdlpZL_Calc(m):
            return m.ReboilerHdlpZL ** 3 + m.ReboilerHdlpa1 * m.ReboilerHdlpZL ** 2 + \
                   m.ReboilerHdlpa2 * m.ReboilerHdlpZL + m.ReboilerHdlpa3 == 0

        hpc_model.ReboilerHdlpZL_Calc = Constraint(rule=ReboilerHdlpZL_Calc)

        def ReboilerHdlpZL_Con1(m):
            return 3 * m.ReboilerHdlpZL ** 2 + 2 * m.ReboilerHdlpa1 * m.ReboilerHdlpZL + m.ReboilerHdlpa2 >= 0

        hpc_model.ReboilerHdlpZL_Con1 = Constraint(rule=ReboilerHdlpZL_Con1)

        def ReboilerHdlpZL_Con2(m):
            return 6 * m.ReboilerHdlpZL + 2 * m.ReboilerHdlpa1 <= 0

        hpc_model.ReboilerHdlpZL_Con2 = Constraint(rule=ReboilerHdlpZL_Con2)

        def ReboilerHdlpZV_Calc(m):
            return m.ReboilerHdlpZV ** 3 + m.ReboilerHdlpb1 * m.ReboilerHdlpZV ** 2 + \
                   m.ReboilerHdlpb2 * m.ReboilerHdlpZV + m.ReboilerHdlpb3 == 0

        hpc_model.ReboilerHdlpZV_Calc = Constraint(rule=ReboilerHdlpZV_Calc)

        def ReboilerHdlpZV_Con1(m):
            return 3 * m.ReboilerHdlpZV ** 2 + 2 * m.ReboilerHdlpb1 * m.ReboilerHdlpZV + m.ReboilerHdlpb2 >= 0

        hpc_model.ReboilerHdlpZV_Con1 = Constraint(rule=ReboilerHdlpZV_Con1)

        def ReboilerHdlpZV_Con2(m):
            return 6 * m.ReboilerHdlpZV + 2 * m.ReboilerHdlpb1 >= 0

        hpc_model.ReboilerHdlpZV_Con2 = Constraint(rule=ReboilerHdlpZV_Con2)

        def ReboilerHdlpZLBL_Calc(m):
            return exp(m.ReboilerHdlpZLBL) - (m.ReboilerHdlpZL - m.ReboilerHdlpBL) == 0

        hpc_model.ReboilerHdlpZLBL_Calc = Constraint(rule=ReboilerHdlpZLBL_Calc)

        def ReboilerHdlpZLdivBL_Calc(m):
            return exp(m.ReboilerHdlpZLdivBL) - (m.ReboilerHdlpZL + 2.414 * m.ReboilerHdlpBL) / (
                        m.ReboilerHdlpZL - 0.414 * m.ReboilerHdlpBL) == 0

        hpc_model.ReboilerHdlpZLdivBL_Calc = Constraint(rule=ReboilerHdlpZLdivBL_Calc)

        def ReboilerHdlpZVBV_Calc(m):
            return exp(m.ReboilerHdlpZVBV) - (m.ReboilerHdlpZV - m.ReboilerHdlpBV) == 0

        hpc_model.ReboilerHdlpZVBV_Calc = Constraint(rule=ReboilerHdlpZVBV_Calc)

        def ReboilerHdlpZVdivBV_Calc(m):
            return exp(m.ReboilerHdlpZVdivBV) - (m.ReboilerHdlpZV + 2.414 * m.ReboilerHdlpBV) / (
                        m.ReboilerHdlpZV - 0.414 * m.ReboilerHdlpBV) == 0

        hpc_model.ReboilerHdlpZVdivBV_Calc = Constraint(rule=ReboilerHdlpZVdivBV_Calc)

        def ReboilerHdlpK_Calc(m, comp):
            return m.ReboilerHdlpPhiL[comp] / m.ReboilerHdlpPhiV[comp] - m.ReboilerHdlpVLE_K[comp] == 0

        hpc_model.ReboilerHdlpK_Calc = Constraint(hpc_model.Component, rule=ReboilerHdlpK_Calc)

        def ReboilerHdlpVLEEqu(m, comp):
            return m.ReboilerLiqLvMFrac[comp] * m.ReboilerHdlpVLE_K[comp] - m.ReboilerVapLvMFrac[comp] == 0

        hpc_model.ReboilerHdlpVLEEqu = Constraint(hpc_model.Component, rule=ReboilerHdlpVLEEqu)

        # --------Reboiler Summation---------
        def ReboilerLiqSum(m):
            return sum([m.ReboilerLiqLvMFrac[c] for c in m.Component]) == 1

        hpc_model.ReboilerLiqSum = Constraint(rule=ReboilerLiqSum)

        def ReboilerVapSum(m):
            return sum([m.ReboilerVapLvMFrac[c] for c in m.Component]) == 1

        hpc_model.ReboilerVapSum = Constraint(rule=ReboilerVapSum)

        # ----------------------------------
        #           GNSplitter
        # ----------------------------------
        # --------GNSplitter Mass Balance---------
        def GNSplitterMassBlnc(m):
            return (m.GNSplitterOutAMFlow + m.GNSplitterOutBMFlow - m.ColumnVapLvMFlow[60]) * 0.001000 == 0

        hpc_model.GNSplitterMassBlnc = Constraint(rule=GNSplitterMassBlnc)


    def build_rto(self, hpc_model, cv_func):
        def obj(m):
            return 4.07 * m.FeedTaMFlow + 1.71 * m.FeedMaMFlow
        hpc_model.obj = Expression(rule=obj)

        def feed_con1(m):
            return -(m.FeedTaMFlow - (m.FeedTaMFlow + m.FeedMaMFlow) * 0.1) * 0.001000
        hpc_model.feed_con1 = Expression(rule=feed_con1)

        def feed_con2(m):
            return (m.FeedTaMFlow - (m.FeedTaMFlow + m.FeedMaMFlow) * 0.4) * 0.001000
        hpc_model.feed_con2 = Expression(rule=feed_con2)

        def purity_con_noise_free(m):
            return 0.999- cv_func['Purity_GAN'].__call__(m)
        hpc_model.purity_con_noise_free = Expression(rule=purity_con_noise_free)

        def purity_con(m):
            return 0.999- cv_func['Purity_GAN'].__call__(m)-m.Purity_GAN_noise
        hpc_model.purity_con = Expression(rule=purity_con)

        def drain_con_noise_free(m):
            return 10 - cv_func['Drain_Flow'].__call__(m)
        hpc_model.drain_con_noise_free = Expression(rule=drain_con_noise_free)

        def drain_con(m):
            return 10 - cv_func['Drain_Flow'].__call__(m)-m.Drain_Flow_noise
        hpc_model.drain_con = Expression(rule=drain_con)