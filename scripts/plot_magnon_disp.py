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


file_path = 'Outputs/numbersIso.txt'

lattice_constant = 1

# Initialize a list to hold the parsed data
datak = []
dataJ = []
# Open the file and read the data
with open(file_path, 'r') as file:
    file.readline()
    for line in file:
        line = line.strip('\n')
        parts = line.split(',')

        kVector = [float(part) for part in parts[:3]]

        # Use regular expression to find all tuples
        matches = re.findall(r'\([-\d.e]+,[-\d.e]+\)', line)
        J = [parse_complex(z) for z in matches]

        datak.append(kVector)
        dataJ.append(J[0])




fig, axs = plt.subplots(ncols=1,nrows=1,figsize=(16/2.25,8/2.52), sharey=True)
plt.style.use("seaborn-v0_8-bright")

num_paths = 6
k_points_per_path = int(len(datak)/num_paths)

x = np.arange(0, len(datak))

split = 5 * k_points_per_path

# the J are given in mRy, convert them to THz by dividing by hbar

ein_mRy_in_eV = 1E-3 * 13.605693
ein_mRy_in_meV = 13.605693

ein_mRy_durch_hbar_in_Hz = 2.067068666810E13 
ein_mRy_durch_hbar_in_THz = 20.67068666810


#axs.plot(x[:split],[J.real for J in dataJ][:split],linewidth=2.5) # plot in mRy
axs.plot(x[:split],[J.real*ein_mRy_durch_hbar_in_THz for J in dataJ][:split],linewidth=2.5) #plot in THz
#axs.plot(x[:split],[J.real*ein_mRy_in_meV for J in dataJ][:split],linewidth=2.5) #plot in meV

labels = [r'$\Gamma$', '$H$', '$N$', r'$\Gamma$', '$P$', '$H$']

for i in range(len(labels)):
    plt.axvline(x=i*k_points_per_path, color='black',linewidth=0.5)


axs.set_xlim(0,(len(labels)-1)*k_points_per_path)


axs.set_xticks([i * k_points_per_path for i in range(len(labels))], labels )



axs.set_ylabel(r'magnon dispersion (THz)')

plt.savefig('scripts/Figures/mag_disp_bccFe.pdf')
plt.show()
plt.clf()


exit()

plt.style.use("seaborn-v0_8-bright")





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
            cy = linex.get_color()
            axs[row][col].plot(paths, paths_yy[col], linestyle='solid', linewidth=lw/8.0, color=cx)

            linez, = axs[row][col].plot(paths, paths_yz[col], label=r'$z$', linestyle='dashed', linewidth=lw)
            cz = linex.get_color()
            axs[row][col].plot(paths, paths_yz[col], linestyle='solid', linewidth=lw/8.0, color=cx)

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

plt.show()
