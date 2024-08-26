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
filepath = "Outputs/full_path/AB_set_pol_vec.txt" # white data
filepath = "Outputs/test_oldMehtod/AB.txt" #colpa data

with open(filepath, 'r') as file:
    file.readline()
    for line in file:
        line = line.strip('\n')
        parts = line.split(' ')
        
        # Use regular expression to find all tuples
        matches = re.findall(r'\([-\d.e]+,[-\d.e]+\)', line)

        if len(matches) == 0:
            break

        AB = [parse_complex(matches[0]), parse_complex(matches[1])]

        As.append(AB[0])
        Bs.append(AB[1])


x = np.linspace(0, 1, len(As))

plt.style.use('science')

fig, axs = plt.subplots(ncols=1,nrows=3,figsize=(16/2.52,8/2.52))

axs[0].plot(x, [A.real for A in As], label="real")
axs[0].plot(x, [A.imag for A in As], label="imag")
axs[1].plot(x, [B.real for B in Bs], label="real")
axs[1].plot(x, [B.imag for B in Bs], label="imag")
axs[2].plot(x, [np.abs(B.imag)/np.abs(A.imag) for A, B in zip(As, Bs)], label=r'$\frac{|A|}{|B|}$')
axs[2].set_ylim(0, 5)


axs[0].set_ylabel(r'$A$')
axs[1].set_ylabel(r'$B$')
axs[0].legend()
axs[1].legend()
plt.tight_layout()
#plt.savefig('scripts/Figures/AB.png', dpi=600)
plt.show()



        