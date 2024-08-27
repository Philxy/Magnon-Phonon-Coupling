import numpy as np
import matplotlib.pyplot as plt



import numpy as np
import matplotlib.pyplot as plt
import scienceplots
from scipy.interpolate import griddata

plt.style.use('science')

data_ev = "Outputs/colpa_grid_diag/ang_eig_fromEV.txt"
data_phonon = "Outputs/colpa_grid_diag/diag_grid_ph_disp_acc.txt"

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


#Lz = [L[i][2][7] for i in range(len(L))]
#plt.scatter(kx, kz, c=Lz, cmap='coolwarm', s=10, marker='s', norm=plt.Normalize(-1,1))
#plt.show()



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

for j in [4,5,6,7]:
    Lx_point = [Lx_all_modes[i][j] for i in range(len(L))]
    Ly_point = [Ly_all_modes[i][j] for i in range(len(L))]
    Lz_point = [Lz_all_modes[i][j] for i in range(len(L))]

    kx = np.array(kx)
    kz = np.array(kz)
    Lz_point = np.array(Lz_point)

    if j == 7:
        Lz_point = [np.abs(Lz_point[i]) for i in range(len(L))]
        axs[0][0].scatter(kx, kz, c=Lz_point, cmap='coolwarm', s=10, marker='s', norm=plt.Normalize(-1,1))
    if j == 5:
        Lz_point = [np.abs(Lz_point[i]) for i in range(len(L))]
        axs[0][1].scatter(kx, kz, c=Lz_point, cmap='coolwarm', s=10, marker='s', norm=plt.Normalize(-1,1))
    if j == 6:
        Lz_point = [-np.abs(Lz_point[i]) for i in range(len(L))]
        axs[1][0].scatter(kx, kz, c=Lz_point, cmap='coolwarm', s=10, marker='s', norm=plt.Normalize(-1,1))
    if j == 4:
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
#plt.savefig("scripts/Figures/ang_eig_fromEV_diag_acc.png", dpi=500)
plt.show()
plt.clf()


from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

#fig = plt.figure(figsize=(10, 7))
#ax = fig.add_subplot(111, projection='3d')
fig, axs = plt.subplots(2, 2, figsize=(22, 7), subplot_kw={'projection': '3d'})

kx_small = []
ky_small = []
kz_small = []
Lz_small_all_modes = []

for i in range(len(L)):
    kx_temp = kVec[i][0]/(2*np.pi)
    ky_temp = kVec[i][1]/(2*np.pi)
    kz_temp = kVec[i][2]/(2*np.pi)
    if np.abs(kx_temp) < 0.1 and np.abs(ky_temp) < 0.1 and np.abs(kz_temp) < 0.1:
        kx_small.append(kx_temp)
        ky_small.append(ky_temp)
        kz_small.append(kz_temp)
        Lz_small_all_modes.append(Lz_all_modes[i])

kx = kx_small
ky = ky_small
kz = kz_small
Lz_all_modes = Lz_small_all_modes


for j in [4,5,6,7]:
    Lx_point = [Lx_all_modes[i][j] for i in range(len(kx))]
    Ly_point = [Ly_all_modes[i][j] for i in range(len(kx))]
    Lz_point = [Lz_all_modes[i][j] for i in range(len(kx))]

    kx = np.array(kx)
    ky = np.array(ky)
    kz = np.array(kz)
    Lz_point = np.array(Lz_point)

    if j == 7:
        Lz_point = [-np.abs(Lz_point[i]) for i in range(len(kx))]
        axs[0][0].scatter(xs=kx, ys=ky, zs=kz,  c=Lz_point, cmap='coolwarm', s=10, marker='s', norm=plt.Normalize(-1,1))
    if j == 5:
        Lz_point = [-np.abs(Lz_point[i]) for i in range(len(kx))]
        axs[0][1].scatter(xs=kx, ys=ky, zs=kz, c=Lz_point, cmap='coolwarm', s=10, marker='s', norm=plt.Normalize(-1,1))
    if j == 6:
        Lz_point = [np.abs(Lz_point[i]) for i in range(len(kx))]
        axs[1][0].scatter(xs=kx, ys=ky, zs=kz, c=Lz_point, cmap='coolwarm', s=10, marker='s', norm=plt.Normalize(-1,1))
    if j == 4:
        Lz_point = [np.abs(Lz_point[i]) for i in range(len(kx))]
        axs[1][1].scatter(xs=kx, ys=ky, zs=kz, c=Lz_point, cmap='coolwarm', s=10, marker='s', norm=plt.Normalize(-1,1))


axs[1][1].set_zlabel('$k_z$ ($2\pi/a$)')
axs[1][1].set_ylabel('$k_y$ ($2\pi/a$)')
axs[1][1].set_xlabel('$k_x$ ($2\pi/a$)')

for i in range(2):
    for j in range(2):
        axs[i][j].set_xlim(-0.1,0.1)
        axs[i][j].set_ylim(-0.1,0.1)
        axs[i][j].set_xticks([-0.05,0,0.05])
        axs[i][j].set_yticks([-0.1,0,0.1])

cbar = plt.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-1,1), cmap='coolwarm'), ax=axs, label='angular momentum', pad=-.45)
cbar.set_ticklabels([r'$\hbar$','0',r'$-\hbar$'])
cbar.set_ticks([-1,0,1])
plt.tight_layout()
plt.savefig("scripts/Figures/ang_eig_fromEV_diag_acc_colpa.png", dpi=700)
plt.show()
plt.clf()

exit()

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
