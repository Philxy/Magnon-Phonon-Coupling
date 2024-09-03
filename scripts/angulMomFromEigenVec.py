import matplotlib.pyplot as plt
import numpy as np

import scienceplots


filename = "Outputs/Angular_Gxyz/GHz_ang_eig_fromEV.txt"
filename = "Outputs/AngularMomentum/ang_eig_fromEV_GHz.txt"
filename = "Outputs/test_GH/ang_eig_fromEV.txt"
#filename = "Outputs/GP_path/ang_eig_fromEV.txt"
#filename = 'Outputs/test_GH copy/ang_eig_fromEV.txt'

plt.style.use('science')


with open(filename, 'r') as f:
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

        plt.scatter([line_count for _ in range(8)], Lz, color='tab:blue', s=3)

xmin = 10
xmax = 490

plt.figsize=(12/2.52, 2/2.52)


plt.xlim(xmin,xmax)

plt.xticks([xmin,xmax], labels=["$\Gamma$", "$H$"])

#plt.plot(x_axis, Lx, label='Lx')
#plt.plot(x_axis, Ly, label='Ly')
#plt.scatter(x_axis, Lz, label='Lz')

plt.ylim(-1.1,1.1)

plt.ylabel("angular momentum ($\hbar$)")
plt.tight_layout()
#plt.savefig("scripts/Figures/ang_eig_fromEV_GH_colpa.png", dpi=500)
plt.show()