from rtolib.core import ProblemDescription

symbol_list = {
                'MV': ("u1",),# input
                'CV': ("cost",),# output
                'OBJ': "cost",# objective function
                'CON': (),  # inequality constraint
                'SPEC': None,  # specification
            }
bounds = {
        "u1":(-5, 5),
        "cost":(None, None),
    }
scaling_factors = {
        "u1":1,
        "cost":1,
    }
default_values = {
        "u1":1,
        "cost":0,
    }
unconstrained_quadratic_problem_description=ProblemDescription(symbol_list, bounds, scaling_factors, default_values)

symbol_list = {
                'MV': ("u1",),# input
                'CV': ("cost",),# output
                'OBJ': "cost",# objective function
                'CON': ("con",),  # inequality constraint
                'SPEC': None,  # specification
            }
bounds = {
        "u1":(-5, 5),
        "cost":(None, None),
        "con":(None, None),
    }
scaling_factors = {
        "u1":1,
        "cost":1,
        "con":1,
    }
default_values = {
        "u1":1,
        "cost":0,
        "con":0,
    }
constrained_quadratic_problem_description=ProblemDescription(symbol_list, bounds, scaling_factors, default_values)