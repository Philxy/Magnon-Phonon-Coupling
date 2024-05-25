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


eigenenergies_GHx = retrieve_data_from_file('Outputs/GHxyz/eigenenergies_GHx.txt')
eigenenergies_GHy = retrieve_data_from_file('Outputs/GHxyz/eigenenergies_GHy.txt')
eigenenergies_GHz = retrieve_data_from_file('Outputs/GHxyz/eigenenergies_GHz.txt')

horizontal_axis = np.linspace(-1, 1, len(eigenenergies_GHx))

dmi_like_file_GHx = 'Outputs/GHxyz/dmiLike_GHx.txt'
dmi_like_file_GHy = 'Outputs/GHxyz/dmiLike_GHy.txt'
dmi_like_file_GHz = 'Outputs/GHxyz/dmiLike_GHz.txt'

GHx_Dxx, GHx_Dxy, GHx_Dxz, GHx_Dyx, GHx_Dyy, GHx_Dyz = [], [], [], [], [], []
GHy_Dxx, GHy_Dxy, GHy_Dxz, GHy_Dyx, GHy_Dyy, GHy_Dyz = [], [], [], [], [], []
GHz_Dxx, GHz_Dxy, GHz_Dxz, GHz_Dyx, GHz_Dyy, GHz_Dyz = [], [], [], [], [], []

with open(dmi_like_file_GHx, 'r') as file:  
    file.readline()
    for line in file:
        line = line.strip('\n')
        parts = line.split(' ')

        data = [parse_complex(z) for z in parts if z != '']
        GHx_Dxx.append(data[0])
        GHx_Dxy.append(data[1])
        GHx_Dxz.append(data[2])
        GHx_Dyx.append(data[3])
        GHx_Dyy.append(data[4])
        GHx_Dyz.append(data[5])


with open(dmi_like_file_GHy, 'r') as file:
    file.readline()
    for line in file:
        line = line.strip('\n')
        parts = line.split(' ')

        data = [parse_complex(z) for z in parts if z != '']
        GHy_Dxx.append(data[0])
        GHy_Dxy.append(data[1])
        GHy_Dxz.append(data[2])
        GHy_Dyx.append(data[3])
        GHy_Dyy.append(data[4])
        GHy_Dyz.append(data[5])


with open(dmi_like_file_GHz, 'r') as file:
    file.readline()
    for line in file:
        line = line.strip('\n')
        parts = line.split(' ')

        data = [parse_complex(z) for z in parts if z != '']
        GHz_Dxx.append(data[0])
        GHz_Dxy.append(data[1])
        GHz_Dxz.append(data[2])
        GHz_Dyx.append(data[3])
        GHz_Dyy.append(data[4])
        GHz_Dyz.append(data[5])

import scienceplots
plt.style.use('science')
fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(15/2.52, 14/2.52))

#axs[0][0].text( -.1,8, r'$\Gamma-H_x$')
axs[0][0].plot(horizontal_axis, [d.imag for d in GHx_Dxx],color='r',linewidth=2, linestyle='dotted')
axs[0][0].plot(horizontal_axis, [d.imag for d in GHx_Dxy],color='b',linewidth=2, linestyle='dashed')
axs[0][0].plot(horizontal_axis, [d.imag for d in GHx_Dxz],color='g',linewidth=2, linestyle='dashdot')
#axs[0][1].text( -.1,8, r'$\Gamma-H_x$')
axs[0][1].plot(horizontal_axis, [d.imag for d in GHx_Dyx],color='r',linewidth=2, linestyle='dotted')
axs[0][1].plot(horizontal_axis, [d.imag for d in GHx_Dyy],color='b',linewidth=2, linestyle='dashed')
axs[0][1].plot(horizontal_axis, [d.imag for d in GHx_Dyz],color='g',linewidth=2, linestyle='dashdot')

#axs[1][0].text( -.1,8, r'$\Gamma-H_y$')
axs[1][0].plot(horizontal_axis, [d.imag for d in GHy_Dxx],color='r',linewidth=2, linestyle='dotted')
axs[1][0].plot(horizontal_axis, [d.imag for d in GHy_Dxy],color='b',linewidth=2, linestyle='dashed')
axs[1][0].plot(horizontal_axis, [d.imag for d in GHy_Dxz],color='g',linewidth=2, linestyle='dashdot')
#axs[1][1].text( -.1,8, r'$\Gamma-H_y$')
axs[1][1].plot(horizontal_axis, [d.imag for d in GHy_Dyx],color='r',linewidth=2, linestyle='dotted')
axs[1][1].plot(horizontal_axis, [d.imag for d in GHy_Dyy],color='b',linewidth=2, linestyle='dashed')
axs[1][1].plot(horizontal_axis, [d.imag for d in GHy_Dyz],color='g',linewidth=2, linestyle='dashdot')

#axs[2][0].text( -.1,8, r'$\Gamma-H_z$')
axs[2][0].plot(horizontal_axis, [d.imag for d in GHz_Dxx],color='r',linewidth=2, linestyle='dotted')
axs[2][0].plot(horizontal_axis, [d.imag for d in GHz_Dxy],color='b',linewidth=2, linestyle='dashed')
axs[2][0].plot(horizontal_axis, [d.imag for d in GHz_Dxz],color='g',linewidth=2, linestyle='dashdot')
#axs[2][1].text( -.1,8, r'$\Gamma-H_z$')
axs[2][1].plot(horizontal_axis, [d.imag for d in GHz_Dyx],color='r',linewidth=2, linestyle='dotted')
axs[2][1].plot(horizontal_axis, [d.imag for d in GHz_Dyy],color='b',linewidth=2, linestyle='dashed')
axs[2][1].plot(horizontal_axis, [d.imag for d in GHz_Dyz],color='g',linewidth=2, linestyle='dashdot')

axs[0][1].set_ylim(-10,10)
axs[1][1].set_ylim(-10,10)
axs[2][1].set_ylim(-10,10)
axs[0][0].set_ylim(-10,10)
axs[1][0].set_ylim(-10,10)
axs[2][0].set_ylim(-10,10)

axs[0][0].plot([],[], label='$x$', color='r',linewidth=2, linestyle='dotted')
axs[0][0].plot([],[], label='$y$', color='b',linewidth=2, linestyle='dashed')
axs[0][0].plot([],[], label='$z$', color='g',linewidth=2, linestyle='dashdot')

axs[0][0].legend(title=r'displacement axis $\mu$', ncols=3, loc='upper center', fancybox=True, shadow=True)

axs[0][0].set_xticks([-1,0,1], labels=[r'-$H_x$', r'$\Gamma$', r'$H_x$'])
axs[1][0].set_xticks([-1,0,1], labels=[r'-$H_y$', r'$\Gamma$', r'$H_y$'])
axs[2][0].set_xticks([-1,0,1], labels=[r'-$H_z$', r'$\Gamma$', r'$H_z$'])
axs[0][1].set_xticks([-1,0,1], labels=[r'-$H_x$', r'$\Gamma$', r'$H_x$'])
axs[1][1].set_xticks([-1,0,1], labels=[r'-$H_y$', r'$\Gamma$', r'$H_y$'])
axs[2][1].set_xticks([-1,0,1], labels=[r'-$H_z$', r'$\Gamma$', r'$H_z$'])

axs[0][0].set_ylabel(r'$\mathrm{Im}\mathcal{D}^{x\mu}_\mathbf{k}$')
axs[0][1].set_ylabel(r'$\mathrm{Im}\mathcal{D}^{y\mu}_\mathbf{k}$')
axs[1][0].set_ylabel(r'$\mathrm{Im}\mathcal{D}^{x\mu}_\mathbf{k}$')
axs[1][1].set_ylabel(r'$\mathrm{Im}\mathcal{D}^{y\mu}_\mathbf{k}$')
axs[2][0].set_ylabel(r'$\mathrm{Im}\mathcal{D}^{x\mu}_\mathbf{k}$')
axs[2][1].set_ylabel(r'$\mathrm{Im}\mathcal{D}^{y\mu}_\mathbf{k}$')


plt.tight_layout()
plt.savefig('scripts/Figures/dmiLike_GHxyz.png', dpi=600)
plt.show()


fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10/2.52, 9/2.52), sharex=False, sharey=False,gridspec_kw={'width_ratios': [1, 2]})


axs[0][0].plot(horizontal_axis, [np.abs(A) for A in eigenenergies_GHx], color='tab:blue')
axs[1][0].plot(horizontal_axis, [np.abs(B) for B in eigenenergies_GHy], color='tab:orange')
axs[0][0].set_xlim(0,1)
axs[1][0].set_xlim(0,1)

axs[0][1].plot(horizontal_axis, [np.abs(A) for A in eigenenergies_GHx], color='tab:blue')
axs[1][1].plot(horizontal_axis, [np.abs(B) for B in eigenenergies_GHy], color='tab:orange')
axs[0][1].set_xlim(0,.075)
axs[1][1].set_xlim(0,.075)
axs[0][1].set_ylim(0,4)
axs[1][1].set_ylim(0,4)


axs[0][0].set_xticks([0,1], labels=[r'$\Gamma$', '$H_x$'])
axs[1][0].set_xticks([0,1], labels=[r'$\Gamma$', '$H_y$'])

axs[0][0].set_ylabel(r'dispersion relation (meV)')
axs[1][0].set_ylabel(r'dispersion relation (meV)')

plt.tight_layout()  
plt.savefig('scripts/Figures/eigenenergies_GHxy.png', dpi=600)
plt.show()