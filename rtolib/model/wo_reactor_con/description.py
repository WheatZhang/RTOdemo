from rtolib.core import ProblemDescription


symbol_list = {
                'MV': ("Fa","Fb","Tr"),# input
                'CV': ('XFr_A','XFr_B','XFr_E','XFr_P','XFr_G'),# output
                'OBJ': "cost",# objective function
                'CON': ("con",),  # inequality constraint
                'SPEC': None,  # specification
            }
bounds = {
        "Fa":(3, 4.5),
        "Fb":(6, 11),
        "Tr":(80, 105),
        "cost":(None, None),
        "con":(None, None),
        "XFr_A":(0, 1),
        "XFr_B":(0, 1),
        "XFr_E":(0, 1),
        "XFr_P":(0, 1),
        "XFr_G":(0, 1),
    }
scaling_factors = {
        "Fa":1.5,#10
        "Fb":5,#10
        "Tr":25,#100
        "cost":10,#10
        "con":0.1,
        "XFr_A":1,
        "XFr_B":1,
        "XFr_E":1,
        "XFr_P":1,
        "XFr_G":1,
    }
default_values = {
        "Fa":4,
        "Fb":9,
        "Tr":90,
        "cost":210,
        "con":0.08,
        "XFr_A":0.109,
        "XFr_B":0.447,
        "XFr_E":0.252,
        "XFr_P":0.063,
        "XFr_G":0.105,
    }
default_WOR_description=ProblemDescription(symbol_list, bounds, scaling_factors, default_values)
