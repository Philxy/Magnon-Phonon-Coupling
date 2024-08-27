import matplotlib.pyplot as plt 
import numpy as np
import regex as re
import scienceplots

def parse_complex(tuple_str):
    tuple_str = tuple_str.replace(')', '')
    tuple_str = tuple_str.replace('(', '')
    tuple_str = tuple_str.split(',')
    return complex(float(tuple_str[0]), float(tuple_str[1]))


As, Bs = [], []

filepath = 'scripts/Data/set_pol_vec_input_output/AB_set_pol_vec_GP.txt'
filepath_colpa = "Outputs/colpa_GP/AB.txt"

with open(filepath_colpa, 'r') as file:
    file.readline()
    for line in file:
        line = line.strip('\n')
        parts = line.split(' ')
        
        # Use regular expression to find all tuples
        matches = re.findall(r'\([-\d.e]+,[-\d.e]+\)', line)

        if len(matches) == 0:
            break

        AB = [parse_complex(matches[0]), parse_complex(matches[1])]

        As.append(np.max([np.abs(AB[0]), np.abs(AB[1])]))
        Bs.append(np.min([np.abs(AB[0]), np.abs(AB[1])]))


x = np.linspace(0, 1, len(As))

plt.rc('text', usetex=True)


plt.style.use('science')

fig, axs = plt.subplots(ncols=1,nrows=1,figsize=(8/2.52,4/2.52))

axs.plot(x, [np.abs(A) for A in As], label=r"$|A^+_{\boldsymbol{k}}|$", linewidth=2)
axs.plot(x, [np.abs(B) for B in Bs], label=r"$|A^-_{\boldsymbol{k}}|$", linewidth=2)

plt.ylabel(r'$|A^{\pm}_{\boldsymbol{k}}|$ (meV)')
plt.legend()
plt.xticks([0,1], labels=[r'$\Gamma$', r'$P$'])
plt.tight_layout()
plt.savefig('scripts/Figures/AB_GP_colpa.png', dpi=600)
plt.show()