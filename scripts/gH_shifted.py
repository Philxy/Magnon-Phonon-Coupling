import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
import scienceplots


kVectors = []
energies = []
mag_disp = []
ph_disp = []

with open("Outputs/wholeBZ/eigenenergies_intercept_acc.txt", "r") as file:
    file.readline()
    for line in file:
        line = line.strip()
        data = line.split(',')

        energies.append([float(ev) for ev in data])


with open("Outputs/intercept_acc.txt", "r") as file:
    file.readline()
    for line in file:
        data = line.split(',')

        mag_disp.append(float(data[3].replace('\n','').replace('(','')))
        


with open("Outputs/wholeBZ/grid_formatted.txt", "r") as file:
    file.readline()
    for line in file:
        data = line.split(' ')

        ph_disp.append([float(d) for d in data[3:6]])



with open("Outputs/wholeBZ/grid_formatted.txt", "r") as file:
    file.readline()
    for line in file:
        data = line.split(' ')

        kVectors.append([float(k) for k in data[:3]])



def find_intersection_kVec(kVectors, mag_disp, ph_disp):

    threshold = 1
    kVecs_close_to_interception = []

    for i in range(len(mag_disp)):
        for branch in range(3):
            if abs(mag_disp[i] - ph_disp[i][branch]) < threshold:
                kVecs_close_to_interception.append(kVectors[i])

    return kVecs_close_to_interception


def sample_close_to_some_kVec(kVec, range, num_setps):

    kVectors = []

    for i in np.linspace(-range, range, num_setps):
        for j in np.linspace(-range, range, num_setps):
            for k in np.linspace(-range, range, num_setps):
                kVectors.append([kVec[0]+i, kVec[1]+j, kVec[2]+k])

    return kVectors


def get_kVecs_close_to_intercep(kVecs, threshold, range, num_steps, mag_disp, ph_disp):
    
        kVecs_intersec = find_intersection_kVec(kVecs, mag_disp, ph_disp)
    
        res = []

        for kVec in kVecs_intersec:
            res.extend(sample_close_to_some_kVec(kVec, range, num_steps))
    
        return res


def print_kVec_for_QE_to_file(kVectors, filepath):

    with open(filepath, "w") as file:
        for k in kVectors:
            file.write(f'{k[0]/(2*np.pi)} {k[1]/(2*np.pi)} {k[2]/(2*np.pi)} 1\n')

        

# range near intersecting points: hybridization range is roughly 0.04 * 2pi/a.
# 
#kVecs_close_to_all_intercep = get_kVecs_close_to_intercep(kVectors, threshold=1, range=0.03*2*np.pi, num_steps=5, mag_disp=mag_disp, ph_disp=ph_disp)

#print('num kVecs close to all intersections: ', len(kVecs_close_to_all_intercep))


#print_kVec_for_QE_to_file(kVecs_close_to_all_intercep, "Outputs/wholeBZ/kVecs_close_to_all_intercep.txt")

#interception_kVec = find_intersection_kVec(kVectors, mag_disp, ph_disp)

#print('num intersecting k Points: ' ,len(interception_kVec))

epsilon = 0.01

indices = []

for i in range(len(kVectors)):
    if abs(kVectors[i][1]) < epsilon :
        indices.append(i)

y_coords = []
x_coords = []
z_coords = []

mag = []
ph1 = []
ph2 = []
ph3 = []
ph_disp_diff = []


#en = [energies[i][ind] for i in range(len(energies))]
en = []

energies1, energies2, energies3, energies4 = [], [], [], []


kxs = np.linspace(0, 0.6, 14)
dk = 0.008

def skip(a, kxs, dk):
    for k in kxs:
        if  (a < k+dk) & (a > k-dk):
            return False
    return True

for i in indices:

    kx,ky,kz = kVectors[i]

    if kx < 0: 
        continue

    if skip(kx, kxs, dk):
        continue
            
    y_coords.append(kVectors[i][1])
    x_coords.append(kVectors[i][0])
    z_coords.append(kVectors[i][2])
    en.append(energies[i][6])
    energies1.append(energies[i][4])
    energies2.append(energies[i][5])
    energies3.append(energies[i][6])
    energies4.append(energies[i][7])
    

    mag_dist = np.min([np.abs(mag_disp[i]-energies[i][j]) for j in range(4,8)])
    ph_dist1 = np.min([np.abs(ph_disp[i][0]-energies[i][j]) for j in range(4,8)])
    ph_dist2 = np.min([np.abs(ph_disp[i][1]-energies[i][j]) for j in range(4,8)])
    ph_dist3 = np.min([np.abs(ph_disp[i][2]-energies[i][j]) for j in range(4,8)])

    ph_disp_diff.append(  np.abs(ph_disp[i][0]-ph_disp[i][1]))

    mag.append(mag_dist)
    ph1.append(ph_dist1)
    ph2.append(ph_dist2)
    ph3.append(ph_dist3)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


# Normalize the z coordinates for coloring
norm = Normalize(vmin=0, vmax=6)
cmap = plt.get_cmap('viridis')

# Obtain colors from the colormap
colors = cmap(norm(en))

#ax.scatter(x_coords, z_coords, en, c=colors)

ax.scatter(x_coords, z_coords, energies1, c='r', s=1)
ax.scatter(x_coords, z_coords, energies2, c='g', s=1)
ax.scatter(x_coords, z_coords, energies3, c='b', s=1)
ax.scatter(x_coords, z_coords, energies4, c='y', s=1)


#ax.plot(x_coords, z_coords, en)


'''
x, z = np.array(x_coords).reshape(100,100), np.array(z_coords).reshape(100,100)
en = np.array(en).reshape(100,100)
'''

ax.set_xlabel('$x$')
ax.set_ylabel('$z$')
ax.set_zlabel('dispersion relation (meV)')

plt.show()
plt.clf()

# 0.3

plt.style.use('science')
fig, ax = plt.subplots(figsize=(8/2.52, 8/2.52) )
aim_kx = 0.3

idxs = []
for i in range(len(x_coords)):
    kx = x_coords[i]
    ky = y_coords[i]
    kz = z_coords[i]

    if np.abs(kx - aim_kx) < 0.025:
        idxs.append(i)
        

ax.set_title(r'$k_x = $ '  + str(aim_kx) + '$\cdot 2\pi/a$')
ax.plot([z_coords[i] for i in idxs], [energies1[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.plot([z_coords[i] for i in idxs], [energies2[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.plot([z_coords[i] for i in idxs], [energies3[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.plot([z_coords[i] for i in idxs], [energies4[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.set_xlabel('$k_z$ ($2\pi/a$)')
ax.set_ylabel('dispersion relation (meV)')

# Add inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

x1, x2, y1, y2 = 0.0, 0.22, 1.9, 3
axins = inset_axes(ax, width="50%", height="50%", loc='upper left')
axins.plot([z_coords[i] for i in idxs], [energies1[i] for i in idxs], color='tab:blue', linewidth=1.5)
axins.plot([z_coords[i] for i in idxs], [energies2[i] for i in idxs], color='tab:blue', linewidth=1.5)
axins.plot([z_coords[i] for i in idxs], [energies3[i] for i in idxs], color='tab:blue', linewidth=1.5)
axins.plot([z_coords[i] for i in idxs], [energies4[i] for i in idxs], color='tab:blue', linewidth=1.5)
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.set_xticklabels([])
axins.set_yticklabels([])
mark_inset(ax, axins, loc1=4, loc2=3, fc="none", ec="0.5")

ax.set_xlim(0,.6)
plt.tight_layout()
plt.savefig('scripts/Figures/1d_kz_disp'+str(aim_kx)+'.png', dpi=600)
plt.show()

# 0.2

plt.style.use('science')
fig, ax = plt.subplots(figsize=(8/2.52, 8/2.52) )
aim_kx = 0.2

idxs = []
for i in range(len(x_coords)):
    kx = x_coords[i]
    ky = y_coords[i]
    kz = z_coords[i]

    if np.abs(kx - aim_kx) < 0.025:
        idxs.append(i)
        

ax.set_title(r'$k_x = $ '  + str(aim_kx) + '$\cdot 2\pi/a$')
ax.plot([z_coords[i] for i in idxs], [energies1[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.plot([z_coords[i] for i in idxs], [energies2[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.plot([z_coords[i] for i in idxs], [energies3[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.plot([z_coords[i] for i in idxs], [energies4[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.set_xlabel('$k_z$ ($2\pi/a$)')
ax.set_ylabel('dispersion relation (meV)')

ax.set_xlim(0,.6)
plt.tight_layout()
plt.savefig('scripts/Figures/1d_kz_disp'+str(aim_kx)+'.png', dpi=600)
plt.show()


# 0.1

plt.style.use('science')
fig, ax = plt.subplots(figsize=(8/2.52, 8/2.52) )
aim_kx = 0.1

idxs = []
for i in range(len(x_coords)):
    kx = x_coords[i]
    ky = y_coords[i]
    kz = z_coords[i]

    if np.abs(kx - aim_kx) < 0.025:
        idxs.append(i)
        

ax.set_title(r'$k_x = $ '  + str(aim_kx) + '$\cdot 2\pi/a$')
ax.plot([z_coords[i] for i in idxs], [energies1[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.plot([z_coords[i] for i in idxs], [energies2[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.plot([z_coords[i] for i in idxs], [energies3[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.plot([z_coords[i] for i in idxs], [energies4[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.set_xlabel('$k_z$ ($2\pi/a$)')
ax.set_ylabel('dispersion relation (meV)')

ax.set_xlim(0,.6)
plt.tight_layout()
plt.savefig('scripts/Figures/1d_kz_disp'+str(aim_kx)+'.png', dpi=600)
plt.show()


# 0.05

plt.style.use('science')
fig, ax = plt.subplots(figsize=(8/2.52, 8/2.52) )
aim_kx = 0.05

idxs = []
for i in range(len(x_coords)):
    kx = x_coords[i]
    ky = y_coords[i]
    kz = z_coords[i]

    if np.abs(kx - aim_kx) < 0.025:
        idxs.append(i)
        

ax.set_title(r'$k_x = $ '  + str(aim_kx) + '$\cdot 2\pi/a$')
ax.plot([z_coords[i] for i in idxs], [energies1[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.plot([z_coords[i] for i in idxs], [energies2[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.plot([z_coords[i] for i in idxs], [energies3[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.plot([z_coords[i] for i in idxs], [energies4[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.set_xlabel('$k_z$ ($2\pi/a$)')
ax.set_ylabel('dispersion relation (meV)')

ax.set_xlim(0,.6)
plt.tight_layout()
plt.savefig('scripts/Figures/1d_kz_disp'+str(aim_kx)+'.png', dpi=600)
#plt.show()
plt.clf()


# 0.275

plt.style.use('science')
fig, ax = plt.subplots(figsize=(8/2.52, 8/2.52) )
aim_kx = 0.275

idxs = []
for i in range(len(x_coords)):
    kx = x_coords[i]
    ky = y_coords[i]
    kz = z_coords[i]

    if np.abs(kx - aim_kx) < 0.02:
        idxs.append(i)
        

ax.set_title(r'$k_x = $ '  + str(aim_kx) + '$\cdot 2\pi/a$')
ax.plot([z_coords[i] for i in idxs], [energies1[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.plot([z_coords[i] for i in idxs], [energies2[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.plot([z_coords[i] for i in idxs], [energies3[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.plot([z_coords[i] for i in idxs], [energies4[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.set_xlabel('$k_z$ ($2\pi/a$)')
ax.set_ylabel('dispersion relation (meV)')

# Add inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

x1, x2, y1, y2 = 0.0, 0.3, 1.4, 3
axins = inset_axes(ax, width="60%", height="50%", loc='upper left')
axins.plot([z_coords[i] for i in idxs], [energies1[i] for i in idxs], color='tab:blue', linewidth=1.5)
axins.plot([z_coords[i] for i in idxs], [energies2[i] for i in idxs], color='tab:blue', linewidth=1.5)
axins.plot([z_coords[i] for i in idxs], [energies3[i] for i in idxs], color='tab:blue', linewidth=1.5)
axins.plot([z_coords[i] for i in idxs], [energies4[i] for i in idxs], color='tab:blue', linewidth=1.5)
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.set_xticklabels([])
axins.set_yticklabels([])
mark_inset(ax, axins, loc1=4, loc2=3, fc="none", ec="0.5")



ax.set_xlim(0,.6)
plt.tight_layout()
plt.savefig('scripts/Figures/1d_kz_disp'+str(aim_kx)+'.png', dpi=600)
plt.show()


# negative

plt.style.use('science')
fig, ax = plt.subplots(figsize=(8/2.52, 8/2.52) )
aim_kx = 0.0

idxs = []
for i in range(len(x_coords)):
    kx = x_coords[i]
    ky = y_coords[i]
    kz = z_coords[i]

    if np.abs(kx - aim_kx) < 0.009:
        idxs.append(i)
        

ax.plot([z_coords[i] for i in idxs], [energies1[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.plot([z_coords[i] for i in idxs], [energies2[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.plot([z_coords[i] for i in idxs], [energies3[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.plot([z_coords[i] for i in idxs], [energies4[i] for i in idxs], color='tab:blue', linewidth=1.5)
ax.set_xlabel('$k_z$ ($2\pi/a$)')
ax.set_ylabel('dispersion relation (meV)')

ax.set_xlim(-.6,.6)
plt.tight_layout()
plt.savefig('scripts/Figures/1d_kz_disp'+str(aim_kx)+'.png', dpi=600)
#plt.show()
plt.clf()
