from rtolib.core import ProblemDescription


symbol_list = {
                'MV': ("Fb","Tr"),# input
                'CV': ('XFr_A','XFr_B','XFr_E','XFr_P','XFr_G'),# output
                'OBJ': "cost",# objective function
                'CON': (),  # inequality constraint
                'SPEC': None,  # specification
            }
bounds = {
        "Fb":(3, 6),
        "Tr":(70, 100),
        "cost":(None, None),
        "XFr_A":(0, 1),
        "XFr_B":(0, 1),
        "XFr_E":(0, 1),
        "XFr_P":(0, 1),
        "XFr_G":(0, 1),
    }
scaling_factors = {
        "Fb":10,
        "Tr":100,
        "cost":100,
        "XFr_A":1,
        "XFr_B":1,
        "XFr_E":1,
        "XFr_P":1,
        "XFr_G":1,
    }
default_values = {
        "Fb":5,
        "Tr":80,
        "cost":170,
        "XFr_A":0.109,
        "XFr_B":0.447,
        "XFr_E":0.252,
        "XFr_P":0.063,
        "XFr_G":0.105,
    }
default_WOR_description=ProblemDescription(symbol_list, bounds, scaling_factors, default_values)
