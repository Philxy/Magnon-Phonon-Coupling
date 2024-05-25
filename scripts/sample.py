import matplotlib.pyplot as plt
import numpy as np


kVectors = []

num = 50
range = .5

#for kx in np.linspace(-range, range, num):
#        for kz in np.linspace(-range, range, num):
#            kVectors.append([kx, 0, kz])


#with open("Outputs/wholeBZ/grid_sample.txt", "w") as file:
#
#    for k in kVectors:
#        file.write(f'{k[0]} {k[1]} {k[2]} 1\n')


# Gamma- H
range = 1
for kx,ky,kz in zip(np.linspace(0, range, num),np.linspace(0, range, num), np.linspace(0, range, num)):
            kVectors.append([0, 0, kz])


# H - N 
for kx,ky,kz in zip(np.linspace(0, .5, num),np.linspace(0, .5, num), np.linspace(1, .5, num)):
            kVectors.append([0, ky, kz])

# N - Gamma
for kx,ky,kz in zip(np.linspace(0, 0.5, num),np.linspace(.5, 0, num), np.linspace(.5, 0, num)):
            kVectors.append([0, ky, kz])


# Gamma - P
for kx,ky,kz in zip(np.linspace(0, .5, num),np.linspace(0, .5, num), np.linspace(0, .5, num)):
            kVectors.append([kx, ky, kz])

# P - H 
for kx,ky,kz in zip(np.linspace(.5, 0, num),np.linspace(.5, 0, num), np.linspace(.5, 1, num)):
            kVectors.append([kx, ky, kz])


with open("Outputs/full_path_new/full_path.txt", "w") as file:

    for k in kVectors:
        file.write(f'{k[0]} {k[1]} {k[2]} 1\n')


