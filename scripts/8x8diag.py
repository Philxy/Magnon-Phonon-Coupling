import matplotlib.pyplot as plt
import numpy as np
import re
from util import Power, Sqrt



# Enable LaTeX text rendering globally
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amsfonts}'


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


filename8x8 = "Outputs/4x4x4Eigenenergies.txt"
filename_ph = "Outputs/numbersPh_20x20x20.txt"
filename_mag = "Outputs/numbersJIso_20x20.txt"

eigenenergies = retrieve_data_from_file(filename8x8)

#eigenenergies_ph = retrieve_data_from_file(filename_ph)[1]
#eigenenergies_mag = retrieve_data_from_file(filename_mag)[1]

'''
CDs = retrieve_data_from_file("Outputs/4x4x4CD.txt")
S_matrix = retrieve_data_from_file("Outputs/4x4x4EVec.txt")

Cs = [C[:3] for C in CDs[1]]
Ds = [D[3:] for D in CDs[1]]


C1_real,C2_real,C3_real = [np.abs(c[0].real) for c in Cs], [np.abs(c[1].real) for c in Cs], [np.abs(c[2].real) for c in Cs]
C1_imag,C2_imag,C3_imag = [np.abs(c[0].imag) for c in Cs], [np.abs(c[1].imag) for c in Cs], [np.abs(c[2].imag) for c in Cs]

D1_real,D2_real,D3_real = [np.abs(d[0].real) for d in Ds], [np.abs(d[1].real) for d in Ds], [np.abs(d[2].real) for d in Ds]
D1_imag,D2_imag,D3_imag = [np.abs(d[0].imag) for d in Ds], [np.abs(d[1].imag) for d in Ds], [np.abs(d[2].imag) for d in Ds]

x_Cs = np.linspace(0,1,len(Cs))
x_Ds = np.linspace(0,1,len(Ds))


#plt.plot(x_Cs,C1_real)
#plt.plot(x_Cs,C2_real)
#plt.plot(x_Cs,C3_real)
#plt.plot(x_Cs,C1_imag)
#plt.plot(x_Cs,C2_imag)
#plt.plot(x_Cs,C3_imag)

# Plot C*D
#plt.plot(x_Cs, [(c[0]*d[0]).real for c,d in zip(Cs,Ds)])
#plt.plot(x_Cs, [(c[1]*d[1]).real for c,d in zip(Cs,Ds)])
#plt.plot(x_Cs, [(c[2]*d[2]).real for c,d in zip(Cs,Ds)])

#plt.show()
#plt.clf()
#
#exit()
'''

#x_ph = np.linspace(0,1,len(eigenenergies_ph))
#y_ph1 = [ ev[0] for ev in eigenenergies_ph]
#y_ph2 = [ ev[1] for ev in eigenenergies_ph]
#y_ph3 = [ ev[2] for ev in eigenenergies_ph]
#
#x_mag =np.linspace(0,1,len(eigenenergies_mag))
#y_mag = [ ev[0] for ev in eigenenergies_mag]

x_coupl = np.linspace(0,1,len(eigenenergies))


"""
for branch in [0,1,2]:

    wpY = [wP(ph[branch],mag[0],C[branch],D[branch]).real  for ph,mag,C,D in zip(eigenenergies_ph,eigenenergies_mag,Cs,Ds)]
    wmY = [wM(ph[branch],mag[0],C[branch],D[branch]).real  for ph,mag,C,D in zip(eigenenergies_ph,eigenenergies_mag,Cs,Ds)]

    #plt.scatter(x_coupl, wpY, s=1, color='black')
    #plt.scatter(x_coupl, wmY, s=1, color='black')
"""

for en in eigenenergies:
    print(en)

for i in range(8):
    y_coupl =  [ np.abs(ev[i]) for ev in eigenenergies]
    #art_der_teilchen = [ np.abs(line[i+8*7]) for line in S_matrix[1]] # work kind of
    #art_der_teilchen = [ np.abs(line[i+8*i]) for line in S_matrix[1]]

    plt.scatter(x_coupl,y_coupl, s=1) #

plt.ylabel(r'\text{dispersion relation (meV)}')
plt.savefig('scripts/Figures/8x8diag.pdf', dpi= 1000)
plt.show()
plt.clf()



exit()

for i in range(8):
    for j in range(8):
        y_coupl =  [ np.abs(ev[7]) for ev in eigenenergies]
        art_der_teilchen = [ np.abs(line[i]) for line in S_matrix[1]]

        plt.scatter(x_coupl,y_coupl, s=0.8, c=art_der_teilchen, cmap="rainbow")
        plt.savefig("scripts/Figures/Cmaps/" + str(i) + "_" + str(j) + ".png")
        plt.clf()


#plt.scatter(x_ph,y_ph1, s=0.5, color='red')
#plt.scatter(x_ph,y_ph2, s=0.5, color='red')
#plt.scatter(x_ph,y_ph3, s=0.5, color='red')
#plt.scatter(x_mag,y_mag, s=0.5, color='blue')

plt.show()
plt.clf()

#11,14
for i in range(12,13):


    art_der_teilchen1 = [ np.abs(line[i*8+i]) for line in S_matrix[1]]
    x_art = np.linspace(0,1,len(art_der_teilchen1))
    plt.scatter(x_art, art_der_teilchen1)

plt.show()