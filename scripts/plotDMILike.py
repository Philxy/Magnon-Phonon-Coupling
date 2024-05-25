import matplotlib.pyplot as plt
import numpy as np
import regex as re

Dxx = []
Dxy = []
Dxz = []

Dyx = []
Dyy = []
Dyz = []


def parse_complex(tuple_str):
    tuple_str = tuple_str.replace(')', '')
    tuple_str = tuple_str.replace('(', '')
    tuple_str = tuple_str.split(',')
    return complex(float(tuple_str[0]), float(tuple_str[1]))


with open('Outputs/full_path/dmiLike_newCD.txt', 'r') as file:
    file.readline()
    for line in file:
        line = line.strip('\n')
        parts = line.split(' ')
        
        # Use regular expression to find all tuples
        matches = re.findall(r'\([-\d.e]+,[-\d.e]+\)', line)

        if len(matches) == 0:
            break

        Dxx.append(parse_complex(matches[0]))
        Dxy.append(parse_complex(matches[1]))
        Dxz.append(parse_complex(matches[2]))

        Dyx.append(parse_complex(matches[3]))
        Dyy.append(parse_complex(matches[4]))
        Dyz.append(parse_complex(matches[5]))
        

import scienceplots 
fig, axs = plt.subplots(ncols=1,nrows=2,figsize=(16/2.52,8/2.52))
plt.style.use("seaborn-v0_8-bright")

x = np.arange(0, len(Dxx))
axs[0].plot(x, [Dxx[i].imag for i in range(len(Dxx))], label='Dxx')
axs[0].plot(x, [Dxy[i].imag for i in range(len(Dxy))], label='Dxy')
axs[0].plot(x, [Dxz[i].imag for i in range(len(Dxz))], label='Dxz')
axs[1].plot(x, [Dyx[i].imag for i in range(len(Dxx))], label='Dyx')
axs[1].plot(x, [Dyy[i].imag for i in range(len(Dxy))], label='Dyy')
axs[1].plot(x, [Dyz[i].imag for i in range(len(Dxz))], label='Dyz')

axs[0].legend()
axs[1].legend()
plt.tight_layout()
plt.show()