from rtolib.core import ProblemDescription


symbol_list = {
                'MV': ("Fb1","Fb2","Fb3","Tr1","Tr2","Tr3"),# input
                'CV': (),# output
                'OBJ': "cost",# objective function
                'CON': ("con",),  # inequality constraint
                'SPEC': None,  # specification
            }
bounds = {
        "Fb1":(3, 7),
        "Fb2":(3, 7),
        "Fb3":(3, 7),
        "Tr1":(60, 120),
        "Tr2":(60, 120),
        "Tr3":(60, 120),
        "cost":(None, None),
        "con":(None, None),
    }
scaling_factors = {
        "Fb1":5,
        "Fb2":5,
        "Fb3":5,
        "Tr1":20,
        "Tr2":20,
        "Tr3":20,
        "cost":10,
        "con":1,
    }
default_values = {
        "Fb1":5,
        "Fb2":5,
        "Fb3":5,
        "Tr1":80,
        "Tr2":80,
        "Tr3":80,
    }
default_WOR_description=ProblemDescription(symbol_list, bounds, scaling_factors, default_values)
