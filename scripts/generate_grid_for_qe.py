import numpy as np


step_width = 0.002

grid_size_z = 0.2
grid_size_x = 0.1


k_points = []


for z in np.arange(-grid_size_z,grid_size_z,step_width):
    for x in np.arange(-grid_size_x,grid_size_x,step_width):
        if x == 0 and z == 0:
            continue
        k = [x, 0, z]
        k_points.append(k)


out_file = open("Parameters/k_points.dat", "w")


for k in k_points:
    out_file.write(" ".join(map(str, k)) + " 1 \n")