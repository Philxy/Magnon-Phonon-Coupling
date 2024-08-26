import matplotlib.pyplot as plt
import numpy as np


kx,ky,kz = [],[],[]
ks = []
with open('Parameters/20x20x20_ph_disp.txt', 'r') as file:

    for line in file:
        line = line.strip().split(" ")
        kx.append(float(line[0]))
        ky.append(float(line[1]))
        kz.append(float(line[2]))
        ks.append([float(line[0]),float(line[1]),float(line[2])])


gen_k_points = []

N_i = 20

a_1 = np.array([1.0 / 2.0, 1.0 / 2.0, -1.0 / 2.0])
a_2 = np.array([-1.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0])
a_3 = np.array([1.0 / 2.0, -1.0 / 2.0, 1.0 / 2.0])

b_1 = 2 * np.pi * np.cross(a_2, a_3) / np.dot(a_3, np.cross(a_1, a_2))
b_2 = 2 * np.pi * np.cross(a_3, a_1) / np.dot(a_3, np.cross(a_1, a_2))
b_3 = 2 * np.pi * np.cross(a_1, a_2) / np.dot(a_3, np.cross(a_1, a_2))


for n_i in range(N_i):
    for n_j in range(N_i):
        for n_k in range(N_i):
            k = np.array([0.0,0.0,0.0])
            k += (n_i + (1.0-N_i)/2.0)/N_i*b_1
            k += (n_j + (1.0-N_i)/2.0)/N_i*b_2
            k += (n_k + (1.0-N_i)/2.0)/N_i*b_3

            gen_k_points.append(k/(2*np.pi))


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x = [gen_k_points[i][0] for i in range(len(gen_k_points))]
y = [gen_k_points[i][1] for i in range(len(gen_k_points))]
z = [gen_k_points[i][2] for i in range(len(gen_k_points))]

gen_kpts_out_file = open("Parameters/gen_kpts.txt", "w")

for i in range(len(gen_k_points)):
    gen_kpts_out_file.write(str(gen_k_points[i][0]) + " " + str(gen_k_points[i][1]) + " " + str(gen_k_points[i][2]) + " 1" +  "\n")

ax.scatter(x, y, z, alpha=.25)

plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x = [ks[i][0] for i in range(len(ks))]
y = [ks[i][1] for i in range(len(ks))]
z = [ks[i][2] for i in range(len(ks))]

ax.scatter(x, y, z, alpha=1)

plt.show()


fig = plt.figure()
plt.scatter(x,y)
plt.show()

import csv

irr_pts = []

# Open and read the file
with open("Parameters/irr_20x20x20.txt", "r") as f:
    reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
    for line in reader:
        # Convert each line into a list of floats and append to irr_pts
        irr_pts.append([float(value) for value in line])

symm_file = open("Parameters/tu_graz_symm_bcc.txt", "r")

symm_ops = symm_file.readlines()

import re

refolded_points = []

for pt_count in range(len(irr_pts)):
    for operation in symm_ops:
        operation = operation.strip().replace("\n", "")
        x,y,z = irr_pts[pt_count][0], irr_pts[pt_count][1], irr_pts[pt_count][2]
        operation = operation.replace("x", str(x))
        operation = operation.replace("y", str(y))
        operation = operation.replace("z", str(z))
        new_point = [eval(coord) for coord in operation.split(',')]
        refolded_points.append(new_point)


x_refolded = [refolded_points[i][0] for i in range(len(refolded_points))]
y_refolded = [refolded_points[i][1] for i in range(len(refolded_points))]
z_refolded = [refolded_points[i][2] for i in range(len(refolded_points))]

plt.clf()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(x_refolded, y_refolded, z_refolded, alpha=1)

plt.show()
plt.clf()


matrix = np.array([b_1,b_2,b_3]).transpose()
matrix_inv = np.linalg.inv(matrix)


for i in range(len(refolded_points)):
    refolded_points[i] = np.dot(matrix_inv, refolded_points[i])


x_refolded = [refolded_points[i][0] for i in range(len(refolded_points))]
y_refolded = [refolded_points[i][1] for i in range(len(refolded_points))]
z_refolded = [refolded_points[i][2] for i in range(len(refolded_points))]


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(x_refolded, y_refolded, z_refolded, alpha=1)

plt.show()
plt.clf()
