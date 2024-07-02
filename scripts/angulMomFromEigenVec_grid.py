import numpy as np
import matplotlib.pyplot as plt


data_ev = "Outputs/Grid/ang_eig_fromEV.txt"
data_phonon = "Outputs/Grid/grid_formatted.txt"

L = []
kVec = []


with open(data_ev, 'r') as f:
    for line_count, line in enumerate(f):
        if line == "\n":
            continue
        line = line.strip()
        splitted_line = line.split(" ")

        Lx = []
        Ly = []
        Lz = []

        for j in range(3):
            for i in range(j,8*3,3):
                if j == 0:
                    Lx.append(float(splitted_line[i]))
                if j == 1:
                    Ly.append(float(splitted_line[i]))
                if j == 2:
                    Lz.append(float(splitted_line[i]))

        L.append([Lx, Ly, Lz])


with open(data_phonon, 'r') as f:
    f.readline()
    f.readline()
    for line_count, line in enumerate(f):
        if line == "\n":
            continue
        line = line.strip()
        splitted_line = line.split(" ")

        k = [float(k) for k in splitted_line[:3]]
        kVec.append(k)


print(len(L), len(kVec))
assert len(L) == len(kVec)



fig = plt.figure()
ax = fig.add_subplot(projection='3d')


kx = [k[0] for k in kVec]
kz = [k[2] for k in kVec]

Lx_all_modes = []
Ly_all_modes = []
Lz_all_modes = []

for i in range(len(L)):

    Lx = L[i][0]
    Ly = L[i][1]
    Lz = L[i][2]

    Lx_all_modes.append([Lx[j] for j in range(8)])
    Ly_all_modes.append([Ly[j] for j in range(8)])
    Lz_all_modes.append([Lz[j] for j in range(8)])


for j in range(8):
    Lx_point = [Lx_all_modes[i][j] for i in range(len(L))]
    Ly_point = [Ly_all_modes[i][j] for i in range(len(L))]
    Lz_point = [Lz_all_modes[i][j] for i in range(len(L))]

    ax.scatter(kx, kz, Lz_point, c=Lz_point, cmap='coolwarm', s=10)

ax.set_xlim(-.5, .5)
ax.set_xlabel('$k_x$')
ax.set_zlabel('$k_z$')

plt.show()