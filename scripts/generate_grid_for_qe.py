import numpy as np
from mpl_toolkits.mplot3d import Axes3D


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




# diagonal plane


vec1 = [1, 0, .5]
vec2 = [0, 1, .5]

scalar_range = 0.2
step_width = 0.004

k_plane = []

for scalar1 in np.arange(-scalar_range,scalar_range,step_width):
    for scalar2 in np.arange(-scalar_range,scalar_range,step_width):
        k = [scalar1*vec1[0] + scalar2*vec2[0], scalar1*vec1[1] + scalar2*vec2[1], scalar1*vec1[2] + scalar2*vec2[2]]
        k_plane.append(k)

out_file_plane = open("Parameters/k_points_diag_plane.dat", "w")

for k in k_plane:
    out_file_plane.write(" ".join(map(str, k)) + " 1 \n")


import matplotlib.pyplot as plt

# ... existing code ...

# Plotting the points of the plane in 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Extracting x, y, and z coordinates from k_plane
x = [k[0] for k in k_plane]
y = [k[1] for k in k_plane]
z = [k[2] for k in k_plane]

# Plotting the points
ax.scatter(x, y, z)

# Setting labels and title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Points of the Plane')

# Displaying the plot
plt.show()
