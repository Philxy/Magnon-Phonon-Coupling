import numpy as np
import matplotlib.pyplot as plt
import scienceplots
from scipy.interpolate import griddata

plt.style.use('science')

data_ev = "Outputs/Diag_Grid/ang_eig_fromEV_acc.txt"
data_phonon = "Outputs/Diag_Grid/diag_grid_ph_disp_acc.txt"

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






kx = [k[0]/(2*np.pi) for k in kVec]
ky = [k[1]/(2*np.pi) for k in kVec]
kz = [k[2]/(2*np.pi) for k in kVec]

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

"""
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

for j in [2]:
    Lx_point = [Lx_all_modes[i][j] for i in range(len(L))]
    Ly_point = [Ly_all_modes[i][j] for i in range(len(L))]
    Lz_point = [Lz_all_modes[i][j] for i in range(len(L))]

    if j == 1:
        Lz_point = [np.abs(Lz_point[i]) for i in range(len(L))]

    if j == 2:
        Lz_point = [-np.abs(Lz_point[i]) for i in range(len(L))]
        pass
    if j == 3:
        Lz_point = [np.abs(Lz_point[i]) for i in range(len(L))]


    kx = np.array(kx)
    ky = np.array(ky)
    kz = np.array(kz)
    Lz_point = np.array(Lz_point)

    

    # Create grid data for contour plot
    xi = np.linspace(kx.min(), kx.max(), 100)
    yi = np.linspace(ky.min(), ky.max(), 100)
    zi = np.linspace(kz.min(), kz.max(), 100)
    xi, zi = np.meshgrid(xi, zi)

    # Interpolate the data
    Lz_interpolated = griddata((kx, kz), Lz_point, (xi, zi), method='cubic')

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel('$k_x$ ($2\pi/a$)')
    ax.set_ylabel('$k_z$ ($2\pi/a$)')
    ax.set_zlabel(r'angular momentum ($\hbar$)')
    ax.set_zlim(-1, 1)

    # Plotting the surface
    surf = ax.plot_surface(xi, zi, Lz_interpolated, cmap='coolwarm', edgecolor='none', norm=plt.Normalize(-1,1))

    # Adding a color bar
    cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
    cbar.set_label('angular momentum ($\hbar$)')

    plt.show()
"""

fig, axs = plt.subplots(2, 2, figsize=(14/2.52,14/2.52), sharex=True, sharey=True)

for j in [0,1,2,3]:
    Lx_point = [Lx_all_modes[i][j] for i in range(len(L))]
    Ly_point = [Ly_all_modes[i][j] for i in range(len(L))]
    Lz_point = [Lz_all_modes[i][j] for i in range(len(L))]

    kx = np.array(kx)
    kz = np.array(kz)
    Lz_point = np.array(Lz_point)

    if j == 3:
        Lz_point = [np.abs(Lz_point[i]) for i in range(len(L))]
        axs[0][0].scatter(kx, kz, c=Lz_point, cmap='coolwarm', s=10, marker='s', norm=plt.Normalize(-1,1))
    if j == 1:
        Lz_point = [np.abs(Lz_point[i]) for i in range(len(L))]
        axs[0][1].scatter(kx, kz, c=Lz_point, cmap='coolwarm', s=10, marker='s', norm=plt.Normalize(-1,1))
    if j == 2:
        Lz_point = [-np.abs(Lz_point[i]) for i in range(len(L))]
        axs[1][0].scatter(kx, kz, c=Lz_point, cmap='coolwarm', s=10, marker='s', norm=plt.Normalize(-1,1))
    if j == 0:
        Lz_point = [np.abs(Lz_point[i]) for i in range(len(L))]
        axs[1][1].scatter(kx, kz, c=Lz_point, cmap='coolwarm', s=10, marker='s', norm=plt.Normalize(-1,1))



axs[1][0].set_xlabel('$k_x$ ($2\pi/a$)')
axs[1][0].set_ylabel('$k_z$ ($2\pi/a$)')

axs[0][0].set_ylabel('$k_z$ ($2\pi/a$)')
axs[1][1].set_xlabel('$k_x$ ($2\pi/a$)')

for i in range(2):
    for j in range(2):
        axs[i][j].set_xlim(-0.1,0.1)
        axs[i][j].set_ylim(-0.1,0.1)
        axs[i][j].set_xticks([-0.05,0,0.05])
        axs[i][j].set_yticks([-0.1,0,0.1])


plt.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-1,1), cmap='coolwarm'), ax=axs, label='angular momentum ($\hbar$)', pad=-.45)

plt.tight_layout()
plt.savefig("scripts/Figures/ang_eig_fromEV_diag_acc.png", dpi=500)
plt.show()
plt.clf()


fig = plt.figure()
ax = fig.add_subplot(projection='3d')



for j in [3]:
    Lx_point = [Lx_all_modes[i][j] for i in range(len(L))]
    Ly_point = [Ly_all_modes[i][j] for i in range(len(L))]
    Lz_point = [Lz_all_modes[i][j] for i in range(len(L))]

    kx = np.array(kx)
    ky = np.array(ky)
    kz = np.array(kz)


    ax.set_xlim(-0.1,0.1)
    ax.set_ylim(-0.1,0.1)
    ax.set_zlim(-0.1,0.1)

    ax.scatter(kx, ky, kz, c=Lz_point, cmap='coolwarm', marker='s', norm=plt.Normalize(-1,1))


    plt.show()
