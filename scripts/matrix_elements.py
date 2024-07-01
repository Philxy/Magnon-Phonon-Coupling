from tabulate import tabulate
import scienceplots
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
from util import Power, Sqrt


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


    return rest


filename8x8 = "Outputs/full_path/eigenenergies_set_pol_vec.txt"
filename_ph = "Outputs/numbersPh_20x20x20.txt"
filename_mag = "Outputs/numbersJIso_20x20.txt"

#eigenenergies = retrieve_data_from_file(filename8x8)


S_matrix = retrieve_data_from_file("Outputs/GHxyz/eigenVectors_GHz.txt")

for i in range(8):
    for j in range(8):

        print(S_matrix[100][i+j*8])
    print("\n")


S_inv_matrices = [np.linalg.inv(
    np.array(S_matrix[i]).reshape(8, 8, order='F')) for i in range(len(S_matrix))]
S_matrices = [np.array(S_matrix[i]).reshape(8, 8, order='F')
              for i in range(len(S_matrix))]

wilket = int(700)

# Use f-strings for formatting with 3 decimal places
inv_matrix_as_str = []
for i in range(8):
    row = []
    for j in range(8):
        real = S_matrices[wilket][i][j].real
        imag = S_matrices[wilket][i][j].imag
        row.append(f"{real:.3f}+{imag:.3f}j")
        abs_val = real**2 + imag**2
        if abs_val > 0.01:
            print(f"Element ({i},{j}) has a value of {abs_val}")
    inv_matrix_as_str.append(row)

print("S_inv_matrices:")
print(tabulate(inv_matrix_as_str))

# Apply the same for S_inv_matrices
matrix_as_str = []
for i in range(8):
    row = []
    for j in range(8):
        real = S_inv_matrices[wilket][i][j].real
        imag = S_inv_matrices[wilket][i][j].imag
        row.append(f"{real:.3f}+{imag:.3f}j")
    matrix_as_str.append(row)

print("S_matrices:")
print(tabulate(matrix_as_str))



#test = np.array([1,2,3,4])
#test = test.reshape(2,2, order='F')
#print(test)
