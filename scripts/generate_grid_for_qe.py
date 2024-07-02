import numpy as np


step_width = 0.005

grid_size = 0.2


k_points = []


for z in np.arange(-grid_size,grid_size,step_width):
    for x in np.arange(-grid_size,grid_size,step_width):
        k = [x, 0, z]
        k_points.append(k)


out_file = open("Parameters/k_points.dat", "w")


for k in k_points:
    out_file.write(" ".join(map(str, k)) + " 1 \n")