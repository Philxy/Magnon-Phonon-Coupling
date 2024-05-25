import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
import scienceplots


kVectors = []
energies = []
mag_disp = []
ph_disp = []
angular_momentum = []

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


with open("Outputs/wholeBZ/angularMomentum_intercept_acc.txt", "r") as file:
    file.readline()
    for line in file:
        data = line.split(',')
        data = [d.replace('(','') for d in data]
        data = [d.replace(')','') for d in data]
        data = [d.replace('\n','') for d in data]

        angular_momentum.append([float(l) for i,l in enumerate(data) if l != '' and i%2 == 0])


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

indices = [0,1]

for i in range(len(kVectors)):
    if abs(kVectors[i][1]) < epsilon :
        indices.append(i)


y_coords = [0,0]
x_coords = [0,0]
z_coords = [0,0]

mag = [0,0]
ph1 = [0,0]
ph2 = [0,0]
ph3 = [0,0]
ph_disp_diff = [0,0]
angul1 = [-1.,1.]
angul2 = [-1.,1.]
angul3 = [-1.,1.]
angul4 = [-1.,1.]



#en = [energies[i][ind] for i in range(len(energies))]
en = [0,0]

energies1, energies2, energies3, energies4 = [0,0], [0,0], [0,0], [0,0]



for i in indices[:len(indices)-1]:

    kx,ky,kz = kVectors[i]

  
    y_coords.append(kVectors[i][1]/(2*np.pi))
    x_coords.append(kVectors[i][0]/(2*np.pi))
    z_coords.append(kVectors[i][2]/(2*np.pi))
    #en.append(energies[i][6])
    #energies1.append(energies[i][4])
    #energies2.append(energies[i][5])
    #energies3.append(energies[i][6])
    #energies4.append(energies[i][7])
    

    #mag_dist = np.min([np.abs(mag_disp[i]-energies[i][j]) for j in range(4,8)])
    #ph_dist1 = np.min([np.abs(ph_disp[i][0]-energies[i][j]) for j in range(4,8)])
    #ph_dist2 = np.min([np.abs(ph_disp[i][1]-energies[i][j]) for j in range(4,8)])
    #ph_dist3 = np.min([np.abs(ph_disp[i][2]-energies[i][j]) for j in range(4,8)])

    angul1.append(angular_momentum[i][4]*2)
    angul2.append(angular_momentum[i][5]*2)
    angul3.append(angular_momentum[i][6]*2)
    angul4.append(angular_momentum[i][7]*2)

    #ph_disp_diff.append(  np.abs(ph_disp[i][0]-ph_disp[i][1]))

    #mag.append(mag_dist)
    #ph1.append(ph_dist1)
    #ph2.append(ph_dist2)
    #ph3.append(ph_dist3)



# angular momentum 2d colormap 

fig, ax = plt.subplots(figsize=(10, 8))


from scipy.interpolate import griddata

# Convert lists to numpy arrays for easier handling
y_coords = np.array(y_coords)
z_coords = np.array(z_coords)
x_coords = np.array(x_coords)

angul1 = np.array(angul1)
angul2 = np.array(angul2)
angul3 = np.array(angul3)
angul4 = np.array(angul4)


# Create grid coordinates
yi = np.linspace(min(y_coords), max(y_coords), 102)
xi = np.linspace(min(x_coords), max(x_coords), 102)
zi = np.linspace(min(z_coords), max(z_coords), 102)

# X - Z grid
xi, zi = np.meshgrid(xi, zi)

# Interpolate mag values onto the grid

an1 = griddata((x_coords, z_coords), angul1, (xi, zi), method='nearest')
an2 = griddata((x_coords, z_coords), angul2, (xi, zi), method='nearest')
an3 = griddata((x_coords, z_coords), angul3, (xi, zi), method='nearest')
an4 = griddata((x_coords, z_coords), angul4, (xi, zi), method='nearest')

# Plotting
plt.style.use('science')


fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(10/2.52, 12/2.52), sharex=True, sharey=True)

from matplotlib.colors import Normalize

cmap = plt.get_cmap('coolwarm')
norm = Normalize(vmin=-1.0, vmax=1.0)

# mag
#axs[0][0].contourf(xi, zi, an1, levels=100, cmap=cmap, norm=norm)
#axs[1][0].contourf(xi, zi, an2, levels=100, cmap=cmap, norm=norm)
#axs[0][1].contourf(xi, zi, an3, levels=100, cmap=cmap, norm=norm)
#axs[1][1].contourf(xi, zi, an4, levels=100, cmap=cmap, norm=norm)


# Attempt to plot the intersection points between mag and ph disp
intersec = []

for i in range(len(ph1)):
    if abs(ph1[i] - mag[i]) < 0.001:
        intersec.append(i)
#axs[0][0].scatter([x_coords[i] for i in intersec], [z_coords[i] for i in intersec], c='black', s=4)




axs[0][0].scatter(x_coords, z_coords, c=angul1,  cmap=cmap, norm=norm)
axs[1][0].scatter(x_coords, z_coords, c=angul2,  cmap=cmap, norm=norm)
axs[0][1].scatter(x_coords, z_coords, c=angul3,  cmap=cmap, norm=norm)
axs[1][1].scatter(x_coords, z_coords, c=angul4,  cmap=cmap, norm=norm)



# Create colorbar with specified ticks using the dummy contour
#cbar = fig.colorbar(contour, ax=ax, ticks=[-.97, 0, 0.97*np.max(angul3)])
#cbar.set_label('angular momentum')
#
#cbar.set_ticklabels([r'$-\hbar$', '0', r'$\hbar$'])

xmin = -.6/(2*np.pi)
ymin = -max(z_coords)
axs[0][0].set_xlim(xmin, -xmin)
axs[0][0].set_ylim(ymin, -ymin)

axs[0][1].set_xlim(xmin, -xmin)
axs[0][1].set_ylim(ymin, -ymin)

axs[1][0].set_xlim(xmin, -xmin)
axs[1][0].set_ylim(ymin, -ymin)

axs[1][1].set_xlim(xmin, -xmin)
axs[1][1].set_ylim(ymin, -ymin)

axs[1][0].set_xlabel(r'$k_x$ ($2\pi/a$)')
axs[1][1].set_xlabel(r'$k_x$ ($2\pi/a$)')
axs[0][0].set_ylabel(r'$k_z$ ($2\pi/a$)')
axs[1][0].set_ylabel(r'$k_z$ ($2\pi/a$)')


dummy = axs[0][0].scatter([0,0],[0,0], s=0, cmap=cmap, norm=norm, c=[-1,1])
#cbar = fig.colorbar(dummy, ticks=[-1, 0, 1])

fig.subplots_adjust(right=.98)
cbar_ax = fig.add_axes([0.98, 0.11, 0.04, .853])
cbar = fig.colorbar(dummy, cax=cbar_ax)

axs[0][0].text(xmin*.8, -ymin*.7, s='a)')
axs[0][1].text(xmin*.8, -ymin*.7, s='b)')
axs[1][0].text(xmin*.8, -ymin*.7, s='c)')
axs[1][1].text(xmin*.8, -ymin*.7, s='d)')

cbar.set_label('angular momentum')
cbar.set_ticks([-1,0,1])
cbar.set_ticklabels([r'$+\hbar$', '0', r'$-\hbar$'])
plt.tight_layout()
plt.savefig('scripts/Figures/2d_angular.png',dpi=600)
plt.show()


'''
fig = plt.figure()
ax = fig.add_subplot()

norm = Normalize(vmin=0, vmax=.11)
cmap = plt.get_cmap('viridis')

ax.scatter(x_coords, z_coords, label='band 1', color=cmap(norm(ph1)))

    #ax.scatter(kVectors[i][0],kVectors[i][2], label='band 2', color=cmap(norm(ph1[i])))
    #ax.scatter(kVectors[i][0],kVectors[i][2], label='band 3', color=cmap(norm(ph2[i])))
    #ax.scatter(kVectors[i][0],kVectors[i][2], label='band 4', color=cmap(norm(ph3[i])))

    #print(kVectors[i][0],kVectors[i][1], kVectors[i][2])

    #print(kVectors[i][0])

plt.show()
'''

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


from scipy.interpolate import griddata

# Convert lists to numpy arrays for easier handling
y_coords = np.array(y_coords)
z_coords = np.array(z_coords)
x_coords = np.array(x_coords)

mag = np.array(mag) 
ph1 = np.array(ph1)
ph2 = np.array(ph2)
ph3 = np.array(ph3)

# Create grid coordinates
yi = np.linspace(min(y_coords), max(y_coords), 100)
xi = np.linspace(min(x_coords), max(x_coords), 100)
zi = np.linspace(min(z_coords), max(z_coords), 100)

# X - Z grid
xi, zi = np.meshgrid(xi, zi)

# Interpolate mag values onto the grid
mag_zi = griddata((x_coords, z_coords), mag, (xi, zi), method='cubic')
ph1_zi = griddata((x_coords, z_coords), ph1, (xi, zi), method='cubic')
ph2_zi = griddata((x_coords, z_coords), ph2, (xi, zi), method='cubic')
ph3_zi = griddata((x_coords, z_coords), ph3, (xi, zi), method='cubic')
diff = griddata((x_coords, z_coords), ph_disp_diff, (xi, zi), method='nearest')

# Plotting
plt.figure(figsize=(10, 8))
plt.style.use('science')

# mag
contour = plt.contourf(xi, zi, mag_zi, levels=50, cmap='viridis')
plt.colorbar(contour, label='energy difference (meV)')
plt.xlabel('X Coordinate')
plt.ylabel('Z Coordinate')
plt.savefig('scripts/Figures/2d_xz_mag.png')
plt.show()

# ph1
contour = plt.contourf(xi, zi, ph1_zi, levels=50, cmap='viridis')
plt.colorbar(contour, label='energy difference (meV)')
plt.xlabel('X Coordinate')
plt.ylabel('Z Coordinate')
plt.savefig('scripts/Figures/2d_xz_ph1.png')
plt.show()


# ph2
contour = plt.contourf(xi, zi, ph2_zi, levels=50, cmap='viridis')
plt.colorbar(contour, label='energy difference (meV)')
plt.xlabel('X Coordinate')
plt.ylabel('Z Coordinate')
plt.savefig('scripts/Figures/2d_xz_ph2.png')
plt.show()


# ph3
contour = plt.contourf(xi, zi, ph3_zi, levels=50, cmap='viridis')
plt.colorbar(contour, label='energy difference (meV)')
plt.xlabel('X Coordinate')
plt.ylabel('Z Coordinate')
plt.savefig('scripts/Figures/2d_xz_ph3.png')
plt.show()

#diff 
contour = plt.contourf(xi, zi, diff, levels=50, cmap='viridis')
plt.colorbar(contour, label='energy difference (meV)')
plt.xlabel('X Coordinate')
plt.ylabel('Z Coordinate')
plt.savefig('scripts/Figures/2d_xz_diff.png')
plt.show()

