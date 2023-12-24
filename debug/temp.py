from rtolib.util.misc import square_circle_mapping
import numpy as np

x = np.array([1,-1,1])
radius = 3
circle_scaling = np.array([1,1,1])
result = square_circle_mapping(x, radius, circle_scaling)
print(np.linalg.norm(result))
print(result)