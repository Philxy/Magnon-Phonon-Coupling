import numpy as np
import matplotlib.pyplot as plt


import re
import scienceplots

plt.style.use('science')


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
        file.readline()
        for line in file:
            parsed_line = parse_complex_numbers(line)
            data.append(parsed_line)

    return data 

# Example usage
filename = 'Outputs/full_path/angularMomentum_set_pol_vec.txt'
angular_momentum_data = read_angular_momentum_components(filename)

x = np.arange(0, len(angular_momentum_data))

dict = {0:r'$\alpha_{\boldsymbol{k}1}^{(\dagger)}$',1:r'$\alpha_{\boldsymbol{k}2}^{(\dagger)}$',2:r'$\alpha_{\boldsymbol{k}3}^{(\dagger)}$', 3: r'$\beta{\boldsymbol{k}}^{(\dagger)}$'}

for j in range(4):
    plt.plot(x,[ 2* (angular_momentum_data[i][j]).real for i in range(len(angular_momentum_data))],linewidth=2)

plt.ylabel(r'angular momentum $\langle  L^{z} \rangle $ ($\hbar$)')
plt.yticks([-1,0,1])
plt.tight_layout()

plt.show()
plt.clf()




#plt.scatter(x, [np.sum(L) for L in angular_momentum_data])
#plt.show()