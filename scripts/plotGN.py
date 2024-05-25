import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
from util import Power, Sqrt
import scienceplots



# Define a function to convert tuple strings to complex numbers
def parse_complex(tuple_str):
    tuple_str = tuple_str.replace(')', '')
    tuple_str = tuple_str.replace('(', '')
    tuple_str = tuple_str.split(',')
    return complex(float(tuple_str[0]), float(tuple_str[1]))


def retrieve_data_from_file(file_path, skip_first=True):

    rest = []

    with open(file_path, 'r') as file:
        if skip_first:
            file.readline()

        for line in file:
            line = line.strip('\n')
            parts = line.split(',')

            # Use regular expression to find all tuples of complex numbers
            matches = re.findall(r'\([-\d.e]+,[-\d.e]+\)', line)

            if len(matches) != 0:
                rest.append([parse_complex(z) for z in matches])
            else:
                rest.append([float(z) for z in parts])

    return rest



def retrieve_data_from_file_space(file_path, skip_first=True):

    rest = []

    with open(file_path, 'r') as file:
        if skip_first:
            file.readline()
        for line in file:
            line = line.strip('\n')
            parts = line.split(' ')

            # Use regular expression to find all tuples of complex numbers
            matches = re.findall(r'\([-\d.e]+,[-\d.e]+\)', line)

            if len(matches) != 0:
                rest.append([parse_complex(z) for z in matches])
            else:
                rest.append([float(z) for z in parts])

    return rest



eigenenergies_path = 'Outputs/GN_path/eigenenergies_GN.txt'
ph_dispersion_path = 'Outputs/GN_path/GN_phdisp.txt'

eigenenergies = retrieve_data_from_file(eigenenergies_path)

phdisp = retrieve_data_from_file_space(ph_dispersion_path)

magdisp = retrieve_data_from_file('Outputs/GN_path/mag.txt')


plt.style.use('science')


path = np.linspace(0, 1, len(eigenenergies))

fig, axs = plt.subplots(1, 2, figsize=(16/2.52, 5/2.52), gridspec_kw={'width_ratios': [1, 4]})

# Second subplot
inset = fig.add_subplot(axs[1])

# Inset for the second subplot
axin = inset.inset_axes([0.6, 0.1, 0.35, 0.76])

#ax2.indicate_inset_zoom(axin1, edgecolor="black" loc1=1)
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
mark_inset(axs[1], axin, loc1=2, loc2=3, lw=.4, alpha=0.5)


axin.plot(path, [ev[4] for ev in eigenenergies], color='tab:blue', lw = 1.5)
axin.plot(path, [ev[5] for ev in eigenenergies], color='tab:blue', lw = 1.5)
axin.plot(path, [ev[6] for ev in eigenenergies], color='tab:blue', lw = 1.5)
axin.plot(path, [ev[7] for ev in eigenenergies], color='tab:blue', lw = 1.5)

axin.set_xlim(0.048,0.105)
axin.set_ylim(1,3)

lw = 1.5
lw2 = 1
axs[0].plot(path, [ev[4] for ev in eigenenergies], color='tab:blue', lw = lw)
axs[0].plot(path, [ev[5] for ev in eigenenergies], color='tab:blue', lw = lw)
axs[0].plot(path, [ev[6] for ev in eigenenergies], color='tab:blue', lw = lw)
axs[0].plot(path, [ev[7] for ev in eigenenergies], color='tab:blue', lw = lw)

axs[0].set_ylabel('dispersion relation (meV)')
axs[0].set_xticks([0,1], labels=[r'$\Gamma$', 'N'])

axs[1].plot(path, [ev[4] for ev in eigenenergies], color='tab:blue', lw = lw)
axs[1].plot(path, [ev[5] for ev in eigenenergies], color='tab:blue', lw = lw)
axs[1].plot(path, [ev[6] for ev in eigenenergies], color='tab:blue', lw = lw)
axs[1].plot(path, [ev[7] for ev in eigenenergies], color='tab:blue', lw = lw)

phdisp.pop(0)
magdisp.pop(0)
axs[1].plot(path, [omega[3] for omega in phdisp], color='tab:orange', lw = lw2, linestyle='--', alpha=0.8)
axs[1].plot(path, [omega[4] for omega in phdisp], color='tab:orange', lw = lw2, linestyle='--', alpha=0.8)
axs[1].plot(path, [omega[5] for omega in phdisp], color='tab:orange', lw = lw2, linestyle='--', alpha=0.8)
axs[1].plot(path, [m[0].real for m in magdisp], color='tab:orange',   lw = lw2, linestyle='--', alpha=0.8)

axin.plot(path, [omega[3] for omega in phdisp], color='tab:orange', lw = lw2, linestyle='--', alpha=0.8)
axin.plot(path, [omega[4] for omega in phdisp], color='tab:orange', lw = lw2, linestyle='--', alpha=0.8)
axin.plot(path, [omega[5] for omega in phdisp], color='tab:orange', lw = lw2, linestyle='--', alpha=0.8)
axin.plot(path, [m[0].real for m in magdisp], color='tab:orange',   lw = lw2, linestyle='--', alpha=0.8)

axs[0].plot(path, [omega[3] for omega in phdisp], color='tab:orange', lw = lw2, linestyle='--', alpha=0.8)
axs[0].plot(path, [omega[4] for omega in phdisp], color='tab:orange', lw = lw2, linestyle='--', alpha=0.8)
axs[0].plot(path, [omega[5] for omega in phdisp], color='tab:orange', lw = lw2, linestyle='--', alpha=0.8)
axs[0].plot(path, [m[0].real for m in magdisp], color='tab:orange',   lw = lw2, linestyle='--', alpha=0.8)

axs[1].set_xlim(0,.3)
axs[1].set_ylim(0,5)

plt.savefig('scripts/Figures/GN_path.png', dpi=600)
plt.show()
        