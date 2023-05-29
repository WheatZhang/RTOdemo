from rtolib.core import ProblemDescription


symbol_list = {
                'MV': ("TA_Flow","MA_Flow"),
                'CV': ('Purity_GAN','T_GAN','T_30thTray','WN_Flow','Drain_Flow'),# output
                'OBJ': "obj",# objective function
                'CON': ("feed_con1","feed_con2","purity_con","drain_con"),  # inequality constraint
                'SPEC': ('GAN_Flow', 'LN_Flow'),  # specification
            }
bounds = {
        "TA_Flow":(200, 400),
        "MA_Flow":(1600, 2000),
        'GAN_Flow':(700,900),
        'LN_Flow':(0,200),
        "Purity_GAN":(0, 1),
        "T_GAN":(None, None),
        "T_30thTray":(None, None),
        "WN_Flow":(0, None),
        "Drain_Flow":(0, None),
        "obj":(None, None),
        "feed_con1":(None, None),
        "feed_con2":(None, None),
        "purity_con":(None, None),
        "drain_con":(None, None),
    }
scaling_factors = {
    "TA_Flow": 100,
    "MA_Flow": 100,
    'GAN_Flow':100,
    'LN_Flow':100,
    "Purity_GAN": 0.01,
    "T_GAN": 10,
    "T_30thTray": 10,
    "WN_Flow": 100,
    "Drain_Flow": 10,
    "obj": 1000,
    "feed_con1": 100,
    "feed_con2": 100,
    "purity_con": 0.01,
    "drain_con": 1,
    }
default_values = {
    "TA_Flow": 500,
    "MA_Flow": 1600,
    'GAN_Flow':800,
    'LN_Flow':100,
    "Purity_GAN": 0.9989999900000225,
    "T_GAN": -178.05452596582163,
    "T_30thTray": -177.2581337959447,
    "WN_Flow": 1127.5420914732154,
    "Drain_Flow": 0,
    "obj": 4400,
    "feed_con1": 200,
    "feed_con2": 200,
    "purity_con": 0.001,
    "drain_con": 0,
    }
default_HPC_2d_description=ProblemDescription(symbol_list, bounds, scaling_factors, default_values)

symbol_list = {
                'MV': ("TA_Flow","MA_Flow", 'GAN_Flow'),
                'CV': ('Purity_GAN','T_GAN','T_30thTray','WN_Flow','Drain_Flow'),# output
                'OBJ': "obj",# objective function
                'CON': ("feed_con1","feed_con2","purity_con","drain_con"),  # inequality constraint
                'SPEC': ('LN_Flow'),  # specification
            }
default_HPC_3d_description=ProblemDescription(symbol_list, bounds, scaling_factors, default_values)

symbol_list = {
                'MV': ("TA_Flow","MA_Flow", 'GAN_Flow', 'LN_Flow'),
                'CV': ('Purity_GAN','T_GAN','T_30thTray','WN_Flow','Drain_Flow'),# output
                'OBJ': "obj",# objective function
                'CON': ("feed_con1","feed_con2","purity_con","drain_con"),  # inequality constraint
                'SPEC': None,  # specification
            }
default_HPC_4d_description=ProblemDescription(symbol_list, bounds, scaling_factors, default_values)
