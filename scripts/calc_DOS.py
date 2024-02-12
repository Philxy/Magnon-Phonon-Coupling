import matplotlib.pyplot as plt
import numpy as np
import re


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


file_path = 'Outputs/numbersBZSampling.txt'

lattice_constant = 1

# Initialize a list to hold the parsed data
datak = []
dataJ = []
# Open the file and read the data
with open(file_path, 'r') as file:
    file.readline()
    for line in file:
        line = line.strip('\n')
        parts = line.split(',')

        kVector = [float(part) for part in parts[:3]]

        # Use regular expression to find all tuples
        matches = re.findall(r'\([-\d.e]+,[-\d.e]+\)', line)
        J = [parse_complex(z) for z in matches]

        datak.append(kVector)
        dataJ.append(J[0])


energy_values = [Jk.real for Jk in dataJ]

num_bins = 100
energy_range = (np.min(energy_values)*0.9, np.max(energy_values)*1.01)  # Energy range


# Calculate the histogram
dos, bin_edges = np.histogram(energy_values, bins=num_bins, range=energy_range, density=True)


# Calculate bin centers from edges
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

bin_centers = bin_edges[1:]

# Plotting the DOS
plt.figure(figsize=(16/2.52, 8/2.52))
plt.plot(bin_centers, dos, drawstyle='steps-mid')
plt.xlabel('Energy')
plt.ylabel('DOS')
plt.tight_layout()
plt.savefig('scripts/Figures/DOS_magnons_bccFe.pdf')
plt.show()