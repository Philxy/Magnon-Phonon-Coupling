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
            parts = line.split(',')

            # Use regular expression to find all tuples of complex numbers
            matches = re.findall(r'\([-\d.e]+,[-\d.e]+\)', line)

            if len(matches) != 0:
                rest.append([parse_complex(z) for z in matches])
            else:
                rest.append([float(z) for z in parts])

    return rest


filename8x8 = "Outputs/gH_1000/Eigenenergies.txt"
filename_ph = "Outputs/numbersPh_20x20x20.txt"
filename_mag = "Outputs/numbersJIso_20x20.txt"

eigenenergies = retrieve_data_from_file(filename8x8)

# eigenenergies_ph = retrieve_data_from_file(filename_ph)[1]
# eigenenergies_mag = retrieve_data_from_file(filename_mag)[1]


CDs = retrieve_data_from_file("Outputs/gH_1000/CD.txt")
S_matrix = retrieve_data_from_file("Outputs/gH_1000/EVec.txt")


Cs = [C[:3] for C in CDs]
Ds = [D[3:] for D in CDs]


C1_real, C2_real, C3_real = [np.abs(c[0].real) for c in Cs], [np.abs(
    c[1].real) for c in Cs], [np.abs(c[2].real) for c in Cs]
C1_imag, C2_imag, C3_imag = [np.abs(c[0].imag) for c in Cs], [np.abs(
    c[1].imag) for c in Cs], [np.abs(c[2].imag) for c in Cs]

D1_real, D2_real, D3_real = [np.abs(d[0].real) for d in Ds], [np.abs(
    d[1].real) for d in Ds], [np.abs(d[2].real) for d in Ds]
D1_imag, D2_imag, D3_imag = [np.abs(d[0].imag) for d in Ds], [np.abs(
    d[1].imag) for d in Ds], [np.abs(d[2].imag) for d in Ds]

x_Cs = np.linspace(0, 1, len(Cs))
x_Ds = np.linspace(0, 1, len(Ds))


# plt.plot(x_Cs,C1_real)
# plt.plot(x_Cs,C2_real)
# plt.plot(x_Cs,C3_real)
# plt.plot(x_Cs,C1_imag)
# plt.plot(x_Cs,C2_imag)
# plt.plot(x_Cs,C3_imag)

# Plot C*D
# plt.plot(x_Cs, [(c[0]*d[0]).real for c,d in zip(Cs,Ds)])
# plt.plot(x_Cs, [(c[1]*d[1]).real for c,d in zip(Cs,Ds)])
# plt.plot(x_Cs, [(c[2]*d[2]).real for c,d in zip(Cs,Ds)])

# plt.show()
# plt.clf()
#
# exit()


# x_ph = np.linspace(0,1,len(eigenenergies_ph))
# y_ph1 = [ ev[0] for ev in eigenenergies_ph]
# y_ph2 = [ ev[1] for ev in eigenenergies_ph]
# y_ph3 = [ ev[2] for ev in eigenenergies_ph]
#
# x_mag =np.linspace(0,1,len(eigenenergies_mag))
# y_mag = [ ev[0] for ev in eigenenergies_mag]

x_coupl = np.linspace(0, 1, len(eigenenergies))


"""
for branch in [0,1,2]:

    wpY = [wP(ph[branch],mag[0],C[branch],D[branch]).real  for ph,mag,C,D in zip(eigenenergies_ph,eigenenergies_mag,Cs,Ds)]
    wmY = [wM(ph[branch],mag[0],C[branch],D[branch]).real  for ph,mag,C,D in zip(eigenenergies_ph,eigenenergies_mag,Cs,Ds)]

    #plt.scatter(x_coupl, wpY, s=1, color='black')
    #plt.scatter(x_coupl, wmY, s=1, color='black')
"""

S_inv_matrices = [np.linalg.inv(
    np.array(S_matrix[i]).reshape(8, 8)) for i in range(len(S_matrix))]
S_matrices = [np.array(S_matrix[i]).reshape(8, 8)
              for i in range(len(S_matrix))]


for i in range(8):
    y_coupl = [np.abs(ev[i]) for ev in eigenenergies]
    # art_der_teilchen = [ np.abs(line[i+8*7]) for line in S_matrix] # work kind of
    # art_der_teilchen = [ np.abs(line[i+i*8]) for line in S_matrix] # work kind of

    # art_der_teilchen = [ np.abs(line[i+8*i]) for line in S_matrix]

    plt.scatter(x_coupl, y_coupl, s=1.4)

plt.ylabel(r'\text{dispersion relation (meV)}')
# plt.savefig('scripts/Figures/8x8diag.pdf', dpi= 1000)
# plt.show()
plt.clf()
it = range(8)


plt.style.use('science')

plt.figure(figsize=(10/2.52, 6/2.52))

min, max = +np.inf, -np.inf

for i in it:
    for j in it:
        art_der_teilchen = [(np.absolute(S_inv[i][j].imag))
                            ** 2 for S_inv in S_inv_matrices]
        curr_max = np.max(art_der_teilchen)
        curr_min = np.min(art_der_teilchen)
        if curr_max > max:
            max = curr_max
        if curr_min < min:
            min = curr_min

        if j == 7 and (i != 6 and i != 7):
            n = 1
            plt.scatter(x_coupl[::n], art_der_teilchen[::n],
                        c=art_der_teilchen[::n], cmap='rainbow', s=1.4)

cb = plt.colorbar()  # orientation='horizontal', location='top', pad=0.003


cb.set_ticks([0, 1])
cb.set_ticklabels([r'$\mathrm{phonon-like}$', r'$\mathrm{magnon-like}$'])

plt.xlim(0.015, 0.095)
plt.ylabel('transformation matrix elements')
plt.tight_layout()
plt.savefig(
    'scripts/Figures/8x8diag_transformation_matrix_elements_1000.png', dpi=600)
plt.show()
plt.clf()

print(min, max)

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
