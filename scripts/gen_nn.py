import numpy as np

# lattice vectors for bcc
a1 = [0.5, 0.5, -0.5]
a2 = [-0.5, 0.5, 0.5]
a3 = [0.5, -0.5, 0.5]

# contains points of the lattice relative to the atom located at (0,0,0)
lattice = []
distances = []  # contains the distances of the atom from the origin in the order given in 'lattice'

n = 8

for i in range(-n, n+1, 1):
    for j in range(-n, n+1, 1):
        for k in range(-n, n+1, 1):
            #if i == 0 and j == 0 and k == 0:  # skip origin
            #    continue
            x = i*a1[0] + j*a2[0] + k*a3[0]
            y = i*a1[1] + j*a2[1] + k*a3[1]
            z = i*a1[2] + j*a2[2] + k*a3[2]
            dist = np.sqrt(x*x + y*y + z*z)
            distances.append(dist)
            lattice.append([x, y, z])


# sort the list containing the atoms of the generated lattice based on the relative distance from the atom at the origin
combined = sorted(zip(distances, lattice))
a_sorted, b_sorted = zip(*combined)
a_sorted = list(a_sorted)
b_sorted = list(b_sorted)


# dictionary given as: 'distance: [[lattice_point],[lattice_point],...]
dict = {}

for i in range(len(distances)):
    d = distances[i]
    dict[d] = []

for i in range(len(distances)):
    d = distances[i]
    dict[d].append(lattice[i])


# sort dict based on the distance from the origin
sorted_dict = {k: dict[k] for k in sorted(dict)}



file = open('Parameters/nn{}.txt'.format(n), 'w')
file.write('x,y,z\n')

for n, d in enumerate(sorted_dict):
    for v in dict[d]:
        file.write(str(v[0]) + "," + str(v[1]) + "," + str(v[2]) + "\n")