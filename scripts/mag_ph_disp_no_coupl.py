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


file_path = 'scripts/Data/set_pol_vec_input_output/magDisp.txt'

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
        dataJ.append(J[0].real)




fig, axs = plt.subplots(ncols=1,nrows=1,figsize=(16/2.52,8/2.52), sharey=True)
plt.style.use("seaborn-v0_8-bright")

num_paths = 5
k_points_per_path = int(len(datak)/num_paths)


split = 5 * k_points_per_path

x = np.linspace(0, 1, split)


# the J are given in mRy, convert them to THz by dividing by hbar

ein_mRy_in_eV = 1E-3 * 13.605693
ein_mRy_in_meV = 13.605693

ein_mRy_durch_hbar_in_Hz = 2.067068666810E13 
ein_mRy_durch_hbar_in_THz = 20.67068666810

axs.plot(x[:split],[J.real for J in dataJ][:split],linewidth=1.5) #plot in meV



file_path = 'scripts/Data/phDispFullPath.txt'

lattice_constant = 1

# Initialize a list to hold the parsed data
dataOmega = []
# Open the file and read the data
with open(file_path, 'r') as file:
    file.readline()
    for line in file:
        line = line.strip('\n')
        parts = line.split()
        omegas = [float(w) for w in parts[1:]]

        dataOmega.append(omegas)


x = np.linspace(0, 1, len(dataOmega))

c = 'tab:orange'
conversion_fac = 4.136/ 33.356 

axs.plot(x, [w[0] * conversion_fac for w in dataOmega], color = c,linewidth=1.5)  # plot in meV
axs.plot(x, [w[1] * conversion_fac for w in dataOmega], color = c,linewidth=1.5)  # plot in meV
axs.plot(x, [w[2] * conversion_fac for w in dataOmega], color = c,linewidth=1.5)  # plot in meV

labels = [r'$\Gamma$', '$H$', '$N$', r'$\Gamma$', '$P$', '$H$']

axs.set_xticks([0, .2, .4, .6, .8, 1])
axs.set_xticklabels(labels)


for i in [0, .2, .4, .6, .8, 1]:
    axs.axvline(x=i, color='black', linestyle='solid', linewidth=0.5)



axs.set_xlim(-0.0,1.0)

axs.set_ylabel(r'$\mathrm{dispersion relation}$ $\mathrm{(meV)}$')

#plt.savefig('scripts/Figures/ph_disp_bccFe_new.pdf')
plt.show()
plt.clf()