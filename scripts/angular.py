import numpy as np
import matplotlib.pyplot as plt


import re
import cmath

def parse_complex_numbers(line):
    # Pattern to match tuples of floats
    pattern = r'\((-?\d+\.?\d*e?-?\d*),(-?\d+\.?\d*e?-?\d*)\)'
    
    # Find all matches in the line
    matches = re.findall(pattern, line)
    
    # Convert matches to complex numbers
    complex_numbers = [complex(float(match[0]), float(match[1])) for match in matches]
    
    return complex_numbers

def read_angular_momentum_components(filename):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            parsed_line = parse_complex_numbers(line)
            data.append(parsed_line)

    return data 

# Example usage
filename = 'Outputs/gH_1000/AngularMomentum.txt'
angular_momentum_data = read_angular_momentum_components(filename)

x = np.arange(0, len(angular_momentum_data))

for j in range(8):
    
    plt.plot(x,[ (angular_momentum_data[i][j]).real for i in range(len(angular_momentum_data))], color='black')

plt.show()
plt.clf()




#plt.scatter(x, [np.sum(L) for L in angular_momentum_data])
#plt.show()