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

filename8x8 = "Outputs/GHxyz/eigenenergies_GHz.txt"
filename_ph = "Outputs/GHxyz/GHz.txt"
filename_mag = "Outputs/GHxyz/mag.txt"

disp_ph = retrieve_data_from_file_space(filename_ph)
disp_mag = retrieve_data_from_file(filename_mag)


eigenenergies = retrieve_data_from_file(filename8x8)

eigenenergies = eigenenergies[int(len(eigenenergies)/2.0)-1:]

x_coupl = np.linspace(0, 1, len(eigenenergies))


plt.style.use('science')

a = .8
lw=1.7

# set up double plots
# Setup figure and grids
fig = plt.figure(figsize=(16/2.52, 5/2.52))  # Wider figure to accommodate both plots
gs = fig.add_gridspec(1, 2, width_ratios=[1, 4])  # Width ratio between first and second plot

# First subplot
ax1 = fig.add_subplot(gs[0])
for i in range(8):
    y_coupl = [np.abs(ev[i]) for ev in eigenenergies]
    ax1.plot(x_coupl, y_coupl, color='tab:blue', linewidth=lw, alpha=1)
ax1.set_xticks([0, 1])
ax1.set_xticklabels([r'$\Gamma$', r'$H$'])
ax1.set_ylabel('dispersion relation (meV)')

# Second subplot
ax2 = fig.add_subplot(gs[1])
for i in range(8):
    y_coupl = [np.abs(ev[i]) for ev in eigenenergies]
    ax2.plot(x_coupl, y_coupl, color='tab:blue', linewidth=2, alpha=1)





# Inset for the second subplot
axin1 = ax2.inset_axes([0.6, 0.1, 0.4, 0.76])
for i in range(8):
    y_coupl = [np.abs(ev[i]) for ev in eigenenergies]
    axin1.plot(x_coupl, y_coupl, color='tab:blue', linewidth=2, alpha=0.9)
axin1.set_ylim(1.75, 3.25)
axin1.set_xlim(0.04, 0.075)
#ax2.indicate_inset_zoom(axin1, edgecolor="black" loc1=1)
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
mark_inset(ax2, axin1, loc1=2, loc2=3, lw=.4, alpha=0.8)


#ax2.set_ylabel('dispersion relation (meV)')


ax2.plot(np.linspace(-1, 1, len(disp_ph)), [e_ph[3]+0.05 for e_ph in disp_ph], color='tab:orange', linewidth=lw, linestyle='dashed', alpha=a)
ax2.plot(np.linspace(-1, 1, len(disp_ph)), [e_ph[4]+0.03 for e_ph in disp_ph], color='tab:orange', linewidth=lw, linestyle='dashed', alpha=a)
ax2.plot(np.linspace(-1, 1, len(disp_ph)), [e_ph[5]+0.07 for e_ph in disp_ph], color='tab:orange', linewidth=lw, linestyle='dashed', alpha=a)
ax2.plot(np.linspace(-1, 1, len(disp_mag)), [e_mag[0].real for e_mag in disp_mag], color='tab:orange', linewidth=lw, linestyle='dashed', alpha=a)

ax1.plot(np.linspace(-1, 1, len(disp_ph)), [e_ph[3] for e_ph in disp_ph], color='tab:orange', linewidth=lw, linestyle='dashed', alpha=a)
ax1.plot(np.linspace(-1, 1, len(disp_ph)), [e_ph[4] for e_ph in disp_ph], color='tab:orange', linewidth=lw, linestyle='dashed', alpha=a)
ax1.plot(np.linspace(-1, 1, len(disp_ph)), [e_ph[5] for e_ph in disp_ph], color='tab:orange', linewidth=lw, linestyle='dashed', alpha=a)
ax1.plot(np.linspace(-1, 1, len(disp_mag)), [e_mag[0].real for e_mag in disp_mag], color='tab:orange', linewidth=lw, linestyle='dashed', alpha=a)


axin1.plot(np.linspace(-1, 1, len(disp_ph)), [e_ph[3]+0.03 for e_ph in disp_ph], color='tab:orange', linewidth=lw, linestyle='dashed', alpha=a)
#axin1.plot(np.linspace(-1, 1, len(disp_ph)), [e_ph[4]+0.03 for e_ph in disp_ph], color='tab:orange', linewidth=lw, linestyle='dashed', alpha=a)
#axin1.plot(np.linspace(-1, 1, len(disp_ph)), [e_ph[5]+1.07 for e_ph in disp_ph], color='tab:orange', linewidth=lw, linestyle='dashed', alpha=a)
axin1.plot(np.linspace(-1, 1, len(disp_mag)), [e_mag[0].real for e_mag in disp_mag], color='tab:orange', linewidth=lw, linestyle='dashed', alpha=a)


ax2.set_xlim(0,0.2)
ax2.set_ylim(-.3,10)
ax1.set_xlim(0,1)

axin1.set_yticks([2,3])

axin1.set_yticks([2,3])

plt.tight_layout()
plt.savefig('scripts/Figures/GH_path.png', dpi=600)
plt.show()







filename8x8 = "Outputs/GHxyz/eigenenergies_GHz.txt"
filename_ph = "Outputs/GHxyz/GHz.txt"
filename_mag = "Outputs/GHxyz/mag.txt"
eigenenergies = retrieve_data_from_file(filename8x8, skip_first=False)

disp_ph = retrieve_data_from_file_space(filename_ph)
disp_mag = retrieve_data_from_file(filename_mag)

disp_ph.pop(0)

print(len(disp_ph), len(disp_mag), len(eigenenergies))

x = np.linspace(-1, 1, len(disp_ph))
y = []

'''
for i in range(len(disp_ph)):
    minimum = 10

    for ev in eigenenergies[i]:
        for omega in disp_ph[i]:
            if abs(ev-omega) < minimum:
                minimum = abs(ev-omega)
        abs_to_mag = np.abs(disp_mag[i][0]-ev)
        if abs_to_mag < minimum:
            minimum = abs_to_mag

    y.append(minimum)
'''




y_phdisp_original = [omega[4] for omega in disp_ph]
y_trans_unpertubed_after_diag = [ev[5] for ev in eigenenergies]

difference = [ np.abs(y_phdisp_original[i]-y_trans_unpertubed_after_diag[i])  for i in range(len(y_phdisp_original))]

fig = plt.figure(figsize=(6/2.52, 5/2.52))



plt.plot(x, difference)
plt.xlim(0,1)
plt.ylabel(r'energy difference (meV)')
plt.xticks([0,1], labels=[r'$\Gamma$', r'$H$'])
plt.tight_layout()
plt.savefig('scripts/Figures/diag_diff_GH_transverse_no_anticross.png', dpi=600)
plt.show()


        