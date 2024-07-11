import matplotlib.pyplot as plt
import numpy as np
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

mi,ma = -2*np.pi,2*np.pi

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

            #k+= n_i*b_1/N_i
            #k+= n_j*b_2/N_i
            #k+= n_k*b_3/N_i

            gen_k_points.append(k/(2*np.pi))


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x = [gen_k_points[i][0] for i in range(len(gen_k_points))]
y = [gen_k_points[i][1] for i in range(len(gen_k_points))]
z = [gen_k_points[i][2] for i in range(len(gen_k_points))]

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