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


file_path = 'Outputs/numbersJsymmTest.txt'

lattice_constant = 1

# Initialize a list to hold the parsed data
datak = []
dataJx = []
dataJy = []
dataJz = []
# Open the file and read the data
with open(file_path, 'r') as file:
    file.readline()
    for line in file:
        line = line.strip('\n')
        parts = line.split(',')

        kVector = [float(part) for part in parts[:3]]

        # Use regular expression to find all tuples
        matches = re.findall(r'\([-\d.e]+,[-\d.e]+\)', line)
        Jx, Jy, Jz = [parse_complex(z) for z in matches][:6]

        datak.append(kVector)
        dataJx.append(Jx)
        dataJy.append(Jy)
        dataJz.append(Jz)


x, y, z = zip(*datak)


import scienceplots

plt.style.use('science')


fig, axs = plt.subplots(ncols=5,nrows=1,figsize=(12/2.25,3/2.52))

plt.subplots_adjust(wspace=0.06, hspace=0.15)

xsplt = split_into_n_parts(dataJx,5)
print(xsplt)
lw = 2
al = 0.8
linex, =axs[0].plot(np.linspace(0,1,len(xsplt[0])), [x.imag for x in xsplt[0]], label=r'$x$',linewidth=lw, alpha=al, linestyle='dotted')
cx = linex.get_color()

axs[1].plot(np.linspace(0,1,len(xsplt[1])), [x.imag for x in xsplt[1]], label=r'$x$', color=cx,linewidth=lw, linestyle='dotted', alpha=al)
axs[2].plot(np.linspace(0,1,len(xsplt[2])), [x.imag for x in xsplt[2]], label=r'$x$', color=cx,linewidth=lw, linestyle='dotted', alpha=al)
axs[3].plot(np.linspace(0,1,len(xsplt[3])), [x.imag for x in xsplt[3]], label=r'$x$', color=cx,linewidth=lw, linestyle='dotted', alpha=al)
axs[4].plot(np.linspace(0,1,len(xsplt[4])), [x.imag for x in xsplt[4]], label=r'$x$', color=cx,linewidth=lw, linestyle='dotted', alpha=al)

ysplt = split_into_n_parts(dataJy,5)

liney, =axs[0].plot(np.linspace(0,1,len(ysplt[0])), [x.imag for x in ysplt[0]], label=r'$y$',linewidth=lw, alpha=al)
cy = liney.get_color()
axs[1].plot(np.linspace(0,1,len(ysplt[1])), [x.imag for x in ysplt[1]], label=r'$y$', color=cy,linewidth=lw, alpha=al)
axs[2].plot(np.linspace(0,1,len(ysplt[2])), [x.imag for x in ysplt[2]], label=r'$y$', color=cy,linewidth=lw, alpha=al)
axs[3].plot(np.linspace(0,1,len(ysplt[3])), [x.imag for x in ysplt[3]], label=r'$y$', color=cy,linewidth=lw, alpha=al)
axs[4].plot(np.linspace(0,1,len(ysplt[4])), [x.imag for x in ysplt[4]], label=r'$y$', color=cy,linewidth=lw, alpha=al)

zsplt = split_into_n_parts(dataJz,5)

linez, =axs[0].plot(np.linspace(0,1,len(zsplt[0])), [x.imag for x in zsplt[0]], label=r'$z$',linewidth=lw, alpha=al, linestyle='dashed')
cz = linez.get_color()
axs[1].plot(np.linspace(0,1,len(zsplt[1])), [x.imag for x in zsplt[1]], label=r'$z$', color=cz,linewidth=lw,  alpha=al, linestyle='dashed')
axs[2].plot(np.linspace(0,1,len(zsplt[2])), [x.imag for x in zsplt[2]], label=r'$z$', color=cz,linewidth=lw,  alpha=al, linestyle='dashed')
axs[3].plot(np.linspace(0,1,len(zsplt[3])), [x.imag for x in zsplt[3]], label=r'$z$', color=cz,linewidth=lw,  alpha=al, linestyle='dashed')
axs[4].plot(np.linspace(0,1,len(zsplt[4])), [x.imag for x in zsplt[4]], label=r'$z$', color=cz,linewidth=lw,  alpha=al, linestyle='dashed')


#axs[0][0].plot(np.linspace(0,1,len(splst[0][0])), , label=r'$J^{x\Gamma}$', color='tab:blue')

axs[0].set_xticks([0,1], [r'$\Gamma$', '$H$'])
axs[1].set_xticks([0,1], ['$H$', '$N$'])
axs[2].set_xticks([0,1], ['$N$',r'$\Gamma$'])
axs[3].set_xticks([0,1], [r'$\Gamma$', '$P$'])
axs[4].set_xticks([0,1], ['$P$', '$H$'])

axs[0].set_xticks([0,1], labels=['',''])
axs[1].set_xticks([0,1], labels=['',''])
axs[2].set_xticks([0,1], labels=['',''])
axs[3].set_xticks([0,1], labels=['',''])
axs[4].set_xticks([0,1], labels=['',''])


axs[0].set_yticks([0,.1,.2])
axs[1].set_yticks([0,.1,.2],labels=['','',''])
axs[2].set_yticks([0,.1,.2],labels=['','',''])
axs[3].set_yticks([0,.1,.2],labels=['','',''])
axs[4].set_yticks([0,.1,.2],labels=['','',''])

axs[0].set_xticks([0,1], [r'$\Gamma$', '$H$'])
axs[1].set_xticks([0,1], ['$H$', '$N$'])
axs[2].set_xticks([0,1], ['$N$',r'$\Gamma$'])
axs[3].set_xticks([0,1], [r'$\Gamma$', '$P$'])
axs[4].set_xticks([0,1], ['$P$', '$H$'])


axs[0].set_ylim(-0.09,0.25)
axs[1].set_ylim(-0.09,0.25)
axs[2].set_ylim(-0.09,0.25)
axs[3].set_ylim(-0.09,0.25)
axs[4].set_ylim(-0.09,0.25)





axs[0].set_ylabel(r'$\mathrm{Im}J^{\mathrm{S-off},\mu}_\mathbf{k}$ ($\frac{\mathrm{meV}}{\mathrm{a.u.}}$)')
#axs[0].legend(title=r'Axis $\mu$')
plt.savefig('scripts/Figures/ImJSoff_bccFe.pdf')
plt.show()
