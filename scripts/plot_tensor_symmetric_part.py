import matplotlib.pyplot as plt
import numpy as np
import re

def parse_complex(tuple_str):
    tuple_str = tuple_str.replace(')', '')
    tuple_str = tuple_str.replace('(', '')
    tuple_str = tuple_str.split(',')
    return complex(float(tuple_str[0]), float(tuple_str[1]))


def retrieve_data_from_file(file_path, skip_first=True):

    rest = []

    with open(file_path, 'r') as file:
        if skip_first:
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


data = retrieve_data_from_file('Outputs/full_path_new/Jsymm.txt')



fig, axs = plt.subplots(1, 5, figsize=(15, 5))


for axis in range(0,3):
    path1 = [d[axis].imag for d in data[:50]]
    path2 = [d[axis].imag for d in data[50:100]]
    path3 = [d[axis].imag for d in data[100:150]]
    path4 = [d[axis].imag for d in data[150:200]]
    path5 = [d[axis].imag for d in data[200:250]]

    axs[0].plot(path1, label='Path 1')
    axs[1].plot(path2, label='Path 2')
    axs[2].plot(path3, label='Path 3')
    axs[3].plot(path4, label='Path 4')
    axs[4].plot(path5, label='Path 5')

for i in range(5):
    axs[i].set_ylim(-.3,.3)

plt.show()