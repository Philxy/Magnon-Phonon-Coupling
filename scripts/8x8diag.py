import matplotlib.pyplot as plt
import numpy as np
import re

# Define a function to convert tuple strings to complex numbers
def parse_complex(tuple_str):
    tuple_str = tuple_str.replace(')', '')
    tuple_str = tuple_str.replace('(', '')
    tuple_str = tuple_str.split(',')
    return complex(float(tuple_str[0]), float(tuple_str[1]))


def retrieve_data_from_file(file_path):

    rest = []
    kVectors = []

    with open(file_path, 'r') as file:
        file.readline()
        for line in file:
            line = line.strip('\n')
            parts = line.split(',')

            kVector = [float(part) for part in parts[:3]]
            kVectors.append(kVector)

            # Use regular expression to find all tuples of complex numbers
            matches = re.findall(r'\([-\d.e]+,[-\d.e]+\)', line)

            if len(matches) != 0:
                rest.append([parse_complex(z) for z in matches])
            else:
                rest.append([float(z) for z in parts[3:]])

    return [kVectors, rest]


filename8x8 = "Outputs/8x8Eigenenergies.txt"

filename_ph = "Outputs/numbersPhTest.txt"
filename_mag = "Outputs/numbersJIso_test.txt"

eigenenergies = retrieve_data_from_file(filename8x8)[1]
eigenenergies_ph = retrieve_data_from_file(filename_ph)[1]
eigenenergies_mag = retrieve_data_from_file(filename_mag)[1]

x_ph = np.linspace(0,1,len(eigenenergies_ph))
y_ph1 = [ ev[0] for ev in eigenenergies_ph]
y_ph2 = [ ev[1] for ev in eigenenergies_ph]
y_ph3 = [ ev[2] for ev in eigenenergies_ph]

x_mag =np.linspace(0,1,len(eigenenergies_mag))
y_mag = [ ev[0] for ev in eigenenergies_mag]

x_coupl = np.linspace(0,1,len(eigenenergies))

for i in range(8):
    y_coupl =  [ np.abs(ev[i]) for ev in eigenenergies]
    plt.scatter(x_coupl,y_coupl, color='green', s=0.5)

plt.scatter(x_ph,y_ph1, s=0.5, color='red')
plt.scatter(x_ph,y_ph2, s=0.5, color='red')
plt.scatter(x_ph,y_ph3, s=0.5, color='red')
plt.scatter(x_mag,y_mag, s=0.5, color='blue')



plt.show()

