from rtolib.core import ProblemDescription


symbol_list = {
                'MV': ("S0","F"),# input
                'CV': [],# output
                'OBJ': "profit",# objective function
                'CON': (),  # inequality constraint
                'SPEC': None,  # specification
            }

bounds={}
scaling_factors={}
for var in ['X', 'P', 'S']:
    for time in range(32):
        if time == 31:
            index=192
            name=str(192)
        else:
            index=time*6.0+6.0
            name=str(time*6+6)
        bounds[var+"_"+name]=(0,None)
        if var == "X":
            scaling_factors[var+"_"+name]=10
        elif var == "P":
            scaling_factors[var+"_"+name]=0.1
        elif var == "S":
            scaling_factors[var+"_"+name]=1



bounds = {
        "Fb":(3, 6),
        "Tr":(70, 100),
        "profit":(None, None),
        "XFr_A":(0, 1),
        "XFr_B":(0, 1),
        "XFr_E":(0, 1),
        "XFr_P":(0, 1),
        "XFr_G":(0, 1),
    }
scaling_factors = {
        "Fb":10,
        "Tr":100,
        "profit":100,
        "XFr_A":1,
        "XFr_B":1,
        "XFr_E":1,
        "XFr_P":1,
        "XFr_G":1,
    }
default_values = {
        "Fb":5,
        "Tr":80,
        "profit":170,
        "XFr_A":0.109,
        "XFr_B":0.447,
        "XFr_E":0.252,
        "XFr_P":0.063,
        "XFr_G":0.105,
    }
default_WOR_description=ProblemDescription(symbol_list, bounds, scaling_factors, default_values)
