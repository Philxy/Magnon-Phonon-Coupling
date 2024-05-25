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


def retrieve_data_from_file(file_path):

    rest = []

    with open(file_path, 'r') as file:
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

filename8x8 = "Outputs/GHx/eigenenergies_GHx.txt"

eigenenergies = retrieve_data_from_file(filename8x8)


x_coupl = np.linspace(-1, 1, len(eigenenergies))

# set up double plots
# Setup figure and grids
fig = plt.figure(figsize=(12/2.52, 6/2.52))  # Wider figure to accommodate both plots
gs = fig.add_gridspec(1, 2, width_ratios=[1, 3.5])  # Width ratio between first and second plot

# First subplot
ax1 = fig.add_subplot(gs[0])
for i in range(8):
    y_coupl = [np.abs(ev[i]) for ev in eigenenergies]
    ax1.plot(x_coupl, y_coupl, color='tab:blue')
ax1.set_xticks([0, 1])
ax1.set_xticklabels([r'$\Gamma$', r'$H$'])
ax1.set_ylabel('dispersion relation (meV)')

# Second subplot
ax2 = fig.add_subplot(gs[1])
for i in range(8):
    y_coupl = [np.abs(ev[i]) for ev in eigenenergies]
    ax2.plot(x_coupl[:int(len(x_coupl)*0.12)], y_coupl[:int(len(x_coupl)*0.12)], color='tab:blue')

ax2.set_xticks([0, 0.05,0.125])


# Inset for the second subplot
axin1 = ax2.inset_axes([0.15, 0.42, 0.4, 0.5])
for i in range(8):
    y_coupl = [np.abs(ev[i]) for ev in eigenenergies]
    axin1.plot(x_coupl[:int(len(x_coupl)*0.12)], y_coupl[:int(len(x_coupl)*0.12)], color='tab:blue', linewidth=1.2)
axin1.set_ylim(2, 3)
axin1.set_xlim(0.05, 0.065)
axin1.set_xticks([])
ax2.indicate_inset_zoom(axin1, edgecolor="black")
ax2.set_ylabel('dispersion relation (meV)')

plt.tight_layout()
plt.savefig('scripts/Figures/diag_double_GHx.png', dpi=600)
plt.show()




dmi_like_file = 'Outputs/GHx/dmiLike_GHx.txt'

Dxx, Dxy, Dxz, Dyx, Dyy, Dyz = [], [], [], [], [], []

with open(dmi_like_file, 'r') as file:  
    file.readline()
    for line in file:
        line = line.strip('\n')
        parts = line.split(' ')

        data = [parse_complex(z) for z in parts if z != '']
        Dxx.append(data[0])
        Dxy.append(data[1])
        Dxz.append(data[2])
        Dyx.append(data[3])
        Dyy.append(data[4])
        Dyz.append(data[5])


plt.plot(x_coupl, [d.imag for d in Dxx], label='Dxx')
plt.plot(x_coupl, [d.imag for d in Dxy], label='Dxy')
plt.plot(x_coupl, [d.imag for d in Dxz], label='Dxz')


plt.legend()
plt.show()

plt.plot(x_coupl, [d.imag for d in Dyx], label='Dyx')
plt.plot(x_coupl, [d.imag for d in Dyy], label='Dyy')
plt.plot(x_coupl, [d.imag for d in Dyz], label='Dyz')
plt.legend()
plt.show()
exit()

# Trying to display the disp rel with correct color bar!

cmap = plt.get_cmap('rainbow')
norm = plt.Normalize(vmin=min*0.99, vmax=max*1.00)







for j in it:
    art_der_teilchen = []

    for i in it:
        y_coupl = [np.abs(ev[i]) for ev in eigenenergies]

        # art_der_teilchen = [ np.abs(line[i+j*8]) for line in S_matrix]

        art_der_teilchen = [(np.absolute(S_inv[i][j].imag))
                            ** 2 for S_inv in S_inv_matrices]
        # art_der_teilchen = [np.abs(S[j][i]) for S in S_matrices]

        data_series = pd.Series(art_der_teilchen, index=x_coupl)

        # Apply a simple moving average with a window of 5
        smooth_data = data_series.rolling(window=1).mean()

        plt.scatter(x_coupl, y_coupl, s=2, c=smooth_data,
                    cmap=cmap, norm=norm)  # RdYlBu

    # plt.savefig("scripts/Figures/Cmaps/" + str(j) + ".png")

    if j == 7:

        plt.ylabel('dispersion relation (meV)')
        #plt.xticks(ticks=[0, 1], labels=[r'$\Gamma$', r'$H$'])
        plt.tight_layout()
        plt.show()

    plt.clf()

exit()
x_min, x_max = 0.04, 0.08

for j in it:
    for i in it:
        y_coupl = [np.abs(ev[i]) for ev in eigenenergies]
        art_der_teilchen = [np.absolute(S_inv[i][j])
                            for S_inv in S_inv_matrices]
        data_series = pd.Series(art_der_teilchen, index=x_coupl)

        # Apply smoothing only within the specified x-axis segment
        smooth_data = data_series.copy()  # Start with a copy of the original series
        segment_mask = (data_series.index >= x_min) & (
            data_series.index <= x_max)
        smooth_data[segment_mask] = data_series[segment_mask].rolling(
            window=10, min_periods=1).mean()

        plt.scatter(x_coupl, y_coupl, s=4, c=smooth_data, cmap=cmap, norm=norm)

    if j == 7:
        plt.show()

    plt.clf()


exit()

for i in range(8):
    for j in range(8):
        y_coupl = [np.abs(ev[7]) for ev in eigenenergies]
        art_der_teilchen = [np.abs(line[i]) for line in S_matrix[1]]

        plt.scatter(x_coupl, y_coupl, s=0.8,
                    c=art_der_teilchen, cmap="rainbow")
        plt.savefig("scripts/Figures/Cmaps/" + str(i) + "_" + str(j) + ".png")
        plt.clf()


# plt.scatter(x_ph,y_ph1, s=0.5, color='red')
# plt.scatter(x_ph,y_ph2, s=0.5, color='red')
# plt.scatter(x_ph,y_ph3, s=0.5, color='red')
# plt.scatter(x_mag,y_mag, s=0.5, color='blue')

plt.show()
plt.clf()

# 11,14
for i in range(12, 13):

    art_der_teilchen1 = [np.abs(line[i*8+i]) for line in S_matrix[1]]
    x_art = np.linspace(0, 1, len(art_der_teilchen1))
    plt.scatter(x_art, art_der_teilchen1)

plt.show()




