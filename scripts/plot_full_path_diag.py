import scienceplots
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
from util import Power, Sqrt


# Enable LaTeX text rendering globally
# plt.rcParams['text.usetex'] = True
# plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amsfonts}'


def wP(AA, BB, CC, DD):
    return Sqrt(Power(AA, 2) + Power(BB, 2) + Sqrt(Power(Power(AA, 2) - Power(BB, 2), 2) + 16*AA*BB*CC*DD))/Sqrt(2)


def wM(AA, BB, CC, DD):
    return Sqrt(Power(AA, 2) + Power(BB, 2) - Sqrt(Power(Power(AA, 2) - Power(BB, 2), 2) + 16*AA*BB*CC*DD))/Sqrt(2)


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

            # Use regular expression to find all tuples of complex numbers
            matches = re.findall(r'\([-\d.e]+,[-\d.e]+\)', line)

            if len(matches) != 0:
                rest.append([parse_complex(z) for z in matches])
            else:
                parts = line.split(',')
                if parts != ['']:
                    rest.append([float(z) for z in parts])

    return rest


filename8x8 = "Outputs/full_path/eigenenergies_set_pol_vec.txt"

eigenenergies = retrieve_data_from_file(filename8x8)


CDs = retrieve_data_from_file("Outputs/full_path/CD_set_pol_vec.txt")
S_matrix = retrieve_data_from_file("Outputs/full_path/eigenVectors_set_pol_vec.txt")


Cs = [C[:3] for C in CDs]
Ds = [D[3:] for D in CDs]

x_Cs = np.linspace(0, 1, len(Cs))
x_Ds = np.linspace(0, 1, len(Ds))


x_coupl = np.linspace(0, 1, len(eigenenergies))


S_inv_matrices = [np.linalg.inv(
    np.array(S_matrix[i]).reshape(8, 8)) for i in range(len(S_matrix))]
S_matrices = [np.array(S_matrix[i]).reshape(8, 8)
              for i in range(len(S_matrix))]

plt.style.use('science')

plt.figure(figsize=(10/2.52, 6/2.52))



for i in range(8):
    y_coupl = [np.abs(ev[i]) for ev in eigenenergies]
    #art_der_teilchen = [ np.abs(line[i+8*7]) for line in S_matrix] # work kind of
    #art_der_teilchen = [ np.abs(line[i]) for line in S_matrix] # work kind of

    art_der_teilchen = [ S[i][0].real for S in S_inv_matrices]

    plt.scatter(x_coupl, y_coupl, c=art_der_teilchen, cmap='viridis', s=1)
    #plt.savefig('scripts/Figures/colors/' + str(ii) + "_" + str(jj)+ '.png', dpi= 200)
    #plt.clf()

plt.ylabel(r'\text{dispersion relation (meV)}')
# plt.savefig('scripts/Figures/asdf.png', dpi= 600)
plt.show()
plt.clf()

fig, axs = plt.subplots(2, 1, figsize=(16/2.52, 16/2.52))

for branch in range(0,3):
    C_values_imag = [(Cs[i][branch].imag) for i in range(len(Cs))]
    C_values_real = [(Cs[i][branch].real) for i in range(len(Cs))]
    axs[0].plot(x_Cs, C_values_imag, label=f'C{branch+1}')
    axs[1].plot(x_Cs, C_values_real, label=f'C{branch+1}')
    D_values = [(Ds[i][branch].imag) for i in range(len(Ds))]
    #plt.plot(x_Ds, D_values, label=f'D{branch+1}')
    #C_times_D_values = [np.abs(C_values[i]*D_values[i]) for i in range(len(C_values))]
    #C_plus_D_values = [np.abs(C_values[i].real) + np.abs(D_values[i].real) for i in range(len(C_values))]
    #plt.scatter(x_Cs, C_plus_D_values, label=f'C{branch+1}D{branch+1}', s=0.5)
    #plt.scatter(x_Cs, C_times_D_values, label=f'C{branch+1}D{branch+1}', s=0.5)

axs[0].set_ylabel(r'\text{Im $C_{\boldsymbol{k}\lambda}$ (meV)}')
axs[1].set_ylabel(r'\text{Re $C_{\boldsymbol{k}\lambda}$ (meV)}')
plt.tight_layout()
plt.savefig('scripts/Figures/coupling_param.png', dpi= 600)
plt.show()


