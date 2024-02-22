import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import re
import matplotlib

# Enable LaTeX text rendering globally
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amsfonts}'


def split_into_n_parts(lst, n):
    # Calculate the size of each chunk
    chunk_size = len(lst) // n
    # Calculate the number of elements that will be left after evenly splitting
    remainder = len(lst) % n
    
    # Initialize variables to keep track of the current index
    start_idx = 0
    parts = []
    
    for i in range(n):
        # If there are leftovers, distribute them one by one to the first 'remainder' chunks
        end_idx = start_idx + chunk_size + (1 if i < remainder else 0)
        # Append the current chunk to the parts list
        parts.append(lst[start_idx:end_idx])
        # Update the start index for the next chunk
        start_idx = end_idx
        
    return parts


# Define a function to convert tuple strings to complex numbers
def parse_complex(tuple_str):
    tuple_str = tuple_str.replace(')', '')
    tuple_str = tuple_str.replace('(', '')
    tuple_str = tuple_str.split(',')
    return complex(float(tuple_str[0]), float(tuple_str[1]))


file_path = 'Outputs/numbersD.txt'

lattice_constant = 1

# Initialize a list to hold the parsed data
datak = []
dataD = []
# Open the file and read the data
with open(file_path, 'r') as file:
    file.readline()
    for line in file:
        line = line.strip('\n')
        parts = line.split(',')

        kVector = [float(part) for part in parts[:3]]

        # Use regular expression to find all tuples
        matches = re.findall(r'\([-\d.e]+,[-\d.e]+\)', line)
        Dxx, Dxy, Dxz, Dyx, Dyy, Dyz = [parse_complex(z) for z in matches][:6]

        Ds = [Dxx, Dxy, Dxz, Dyx, Dyy, Dyz]
        datak.append(kVector)
        dataD.append(Ds)

x, y, z = zip(*datak)


# Display the path in k space
fig = plt.figure()

# Add a 3D subplot
ax = fig.add_subplot(111, projection='3d')

# Plot the points
ax.plot(x, y, z, marker='o')

ax.scatter(0,0,0, label=r'$\Gamma$', color='r') # Gamma Point
ax.text(0,0,0, r'$\Gamma$', color='r')

ax.scatter(0,0, 2*np.pi/lattice_constant, label='H', color='r')
ax.text(0,0,2*np.pi/lattice_constant, 'H', color='r')

ax.scatter(np.pi/lattice_constant,np.pi/lattice_constant, np.pi/lattice_constant, label='P', color='r')
ax.text(np.pi/lattice_constant,np.pi/lattice_constant, np.pi/lattice_constant, 'P', color='r')

ax.scatter(0,np.pi/lattice_constant, np.pi/lattice_constant, label='N', color='r')
ax.text(0,np.pi/lattice_constant, np.pi/lattice_constant, 'N', color='r')

# Adding labels for clarity
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')

# Set title
ax.set_title('3D Path Visualization')

# Show the plot
plt.show()
plt.clf()

fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2,
                               sharex=True, figsize=(16/2.52, 8/2.52))

x = np.linspace(0, 1, len(datak))

xx = [d[0].imag for d in dataD]
xy = [d[1].imag for d in dataD]
xz = [d[2].imag for d in dataD]

yx = [d[3].imag for d in dataD]
yy = [d[4].imag for d in dataD]
yz = [d[5].imag for d in dataD]

# Labeling the x axis
# number_of_labels = int(len(x) / 6.0)+1  # Interval for special labels
# labels = [r'$\Gamma$', 'H', 'N' ,r'$\Gamma$', 'P', 'H']
# label_values = [x[i] for i in range(0, len(x), number_of_labels)]
# plt.xticks(x[::number_of_labels], labels=labels)

# tick_labels = [r'$\Gamma$', 'H', 'N' ,r'$\Gamma$', 'P', 'H', 'P', 'N']
# ax2.set_xticks(ticks=np.append(np.linspace(0,1,6),), labels=tick_labels )

for i in range(0, 7, 1):
    ax1.axvline(i/6.0, color='black', alpha=0.3)
    ax2.axvline(i/6.0, color='black', alpha=0.3)

ax1.plot(x, xx, label=r'$\mu = x$', linestyle='dotted', color='tab:orange')
ax1.plot(x, xy, label=r'$\mu = y$', color='tab:green', linestyle='dashdot')
ax1.plot(x, xz, label=r'$\mu = z$', color='tab:red', linestyle='dashed')

ax1.set_ylabel(r'$\mathrm{Im} \mathcal{D}^{x\mu}_\mathbf{k}$')

ax2.plot(x, yx, label=r'$\mu = x$', linestyle='dotted', color='tab:orange')
ax2.plot(x, yy, label=r'$\mu = y$', color='tab:green', linestyle='dashdot')
ax2.plot(x, yz, label=r'$\mu = z$', color='tab:red', linestyle='dashed')

ax2.set_ylabel(r'$\mathrm{Im}\mathcal{D}^{y\mu}_\mathbf{k}$')

ax1.set_yticks([0, 5])
ax2.set_xticks([])
ax2.legend()
ax1.legend()
plt.savefig('scripts/Figures/bcc_Fe_DSLC_parameter.pdf')
plt.show()
plt.clf()


#from aquarel import load_theme
#
#theme = load_theme("scientific")
#theme.apply()
## ... plotting code here
#theme.apply_transforms()

plt.style.use("seaborn-v0_8-bright")


#plt.plot(x, [Dx-Dy for Dx, Dy in zip(xx,yx)])
#plt.plot(x, [Dx-Dy for Dx, Dy in zip(yx,yy)])
#plt.plot(x, [Dx-Dy for Dx, Dy in zip(xz,yz)])
#plt.show()


paths_xx = split_into_n_parts(xx,6)
paths_xy = split_into_n_parts(xy,6)
paths_xz = split_into_n_parts(xz,6)

paths_yx = split_into_n_parts(yx,6)
paths_yy = split_into_n_parts(yy,6)
paths_yz = split_into_n_parts(yz,6)

path_data = [[paths_xx, paths_xy, paths_xz], [paths_yx,paths_yy,paths_yz]]

paths = np.linspace(0,1,len(paths_xx[0]))

fig, axs = plt.subplots(ncols=6,nrows=2,figsize=(16/2.25,10/2.52))

plt.subplots_adjust(wspace=0.05, hspace=0.06)

lw=2.5

for col in range(6):
    for row in range(2):
        if row == 0:
            linex, = axs[row][col].plot(paths, paths_xx[col], label=r'$x$', linestyle='dotted', linewidth=lw)
            cx = linex.get_color()
            axs[row][col].plot(paths, paths_xx[col], linestyle='solid', linewidth=lw/8.0, color=cx)
            
            liney, = axs[row][col].plot(paths, paths_xy[col], label=r'$y$', linestyle='dashdot', linewidth=lw)
            cy = liney.get_color()
            axs[row][col].plot(paths, paths_xy[col], linestyle='solid', linewidth=lw/8.0, color=cy)

            linez, = axs[row][col].plot(paths, paths_xz[col], label=r'$z$', linestyle='dashed', linewidth=lw)
            cz = linez.get_color()
            axs[row][col].plot(paths, paths_xz[col], linestyle='solid', linewidth=lw/8.0, color=cz)


            axs[row][col].set_xticks([])
            axs[row][col].set_yticks([])
            axs[row][col].set_ylim([-3,10])
        else:
            linex, = axs[row][col].plot(paths, paths_yx[col], label=r'$x$', linestyle='dotted', linewidth=lw)
            cx = linex.get_color()
            axs[row][col].plot(paths, paths_yx[col], linestyle='solid', linewidth=lw/8.0, color=cx)

            liney, = axs[row][col].plot(paths, paths_yy[col], label=r'$y$', linestyle='dashdot', linewidth=lw)
            cy = liney.get_color()
            axs[row][col].plot(paths, paths_yy[col], linestyle='solid', linewidth=lw/8.0, color=cy)

            linez, = axs[row][col].plot(paths, paths_yz[col], label=r'$z$', linestyle='dashed', linewidth=lw)
            cz = linez.get_color()
            axs[row][col].plot(paths, paths_yz[col], linestyle='solid', linewidth=lw/8.0, color=cz)

            axs[row][col].set_yticks([])
            axs[row][col].set_ylim([-9,3])


axs[1][0].set_xticks([0,1], [r'$\Gamma$', '$H$'])
axs[1][1].set_xticks([0,1], ['$H$', '$N$'])
axs[1][2].set_xticks([0,1], ['$N$',r'$\Gamma$'])
axs[1][3].set_xticks([0,1], [r'$\Gamma$', '$P$'])
axs[1][4].set_xticks([0,1], ['$P$', '$H$'])
axs[1][5].set_xticks([0,1], ['$P$', '$N$'])

axs[0][0].set_yticks([0,5])
axs[1][0].set_yticks([-5,0])

axs[0][0].set_ylabel(r'$\mathrm{Im}\mathcal{D}^{x\mu}_\mathbf{k}$')
axs[1][0].set_ylabel(r'$\mathrm{Im}\mathcal{D}^{y\mu}_\mathbf{k}$')


axs[0][5].legend(title=r'$\mu$')

#plt.savefig('scripts/Figures/ImD_bccFe.pdf')
plt.show()
