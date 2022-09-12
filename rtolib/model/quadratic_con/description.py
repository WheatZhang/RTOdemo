from rtolib.core import ProblemDescription


symbol_list = {
                'MV': ("u1","u2"),# input
                'CV': ("cost","con"),# output
                'OBJ': "cost",# objective function
                'CON': ("con",),  # inequality constraint
                'SPEC': None,  # specification
            }
bounds = {
        "u1":(None, None),
        "u2":(None, None),
        "cost":(None, None),
        "con":(None, None),
    }
scaling_factors = {
        "u1":1,
        "u2":1,
        "cost":1,
        "con":1,
    }
default_values = {
        "u1":1,
        "u2":1,
        "cost":0,
        "con":0,
    }
quadratic_con_problem_description=ProblemDescription(symbol_list, bounds, scaling_factors, default_values)
