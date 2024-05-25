import matplotlib.pyplot as plt
import numpy as np
import scienceplots

plt.style.use("science")

def parse_complex(tuple_str):
    tuple_str = tuple_str.replace(')', '')
    tuple_str = tuple_str.replace('(', '')
    tuple_str = tuple_str.split(',')
    return complex(float(tuple_str[0]), float(tuple_str[1]))

As,Bs = [],[]

axis = 'x'

with open('Outputs/GHxyz/AB_GH' + axis + '.txt', 'r') as file:
    file.readline()
    for line in file:
        line = line.strip('\n')
        parts = line.split(' ')
        
        A,B = parse_complex(parts[0]), parse_complex(parts[1])
        As.append(A)
        Bs.append(B)



x = np.linspace(-1,1,len(As))

fig, axs = plt.subplots(ncols=1,nrows=1,figsize=(8/2.52,4/2.52))


plt.plot(x, [np.abs(A) for A in As], label=r'$|A^{+}_{\boldsymbol{k}}|$')
plt.plot(x, [np.abs(B) for B in Bs], label=r'$|A^{-}_{\boldsymbol{k}}|$')

plt.xticks([-1,0,1], labels=[r'$-H_'  + axis + '$', r'$\Gamma$', r'$H_' + axis +  '$'])

plt.legend()

plt.ylabel(r'$|A^{+}_{\boldsymbol{k}}|$ (meV)')
plt.tight_layout()
plt.savefig('scripts/Figures/AB_GH' + axis + '.png', dpi=600)
plt.show()

