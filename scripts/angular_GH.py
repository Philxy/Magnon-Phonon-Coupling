import numpy as np
import matplotlib.pyplot as plt


import re
import scienceplots

plt.style.use('science')


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
filename = 'scripts/Data/set_pol_vec_input_output/angularMomentum_set_pol_vec.txt'
filename = 'Outputs/test_GH/ang_mom.txt'
angular_momentum_data = read_angular_momentum_components(filename)


plt.figsize=(12/2.52, 2/2.52)

x = np.arange(0, len(angular_momentum_data))

dict = {0:r'$\alpha_{\boldsymbol{k}1}^{(\dagger)}$',1:r'$\alpha_{\boldsymbol{k}2}^{(\dagger)}$',2:r'$\alpha_{\boldsymbol{k}3}^{(\dagger)}$', 3: r'$\beta{\boldsymbol{k}}^{(\dagger)}$'}

plt.plot(x,[ 2* (angular_momentum_data[i][0]).real for i in range(len(angular_momentum_data))],linewidth=2, alpha=0.8, linestyle='--')
plt.plot(x,[ 2* (angular_momentum_data[i][1]).real for i in range(len(angular_momentum_data))],linewidth=2, alpha=0.8)
plt.plot(x,[ 2* (angular_momentum_data[i][2]).real for i in range(len(angular_momentum_data))],linewidth=2, alpha=0.8)
plt.plot(x,[ 2* (angular_momentum_data[i][3]).real for i in range(len(angular_momentum_data))],linewidth=2, alpha=0.8)

plt.xticks([0,499], labels=[r'$\Gamma$',r'$H$'])
plt.ylabel(r'angular momentum ($\hbar$)')
plt.yticks([-1,0,1])
plt.xlim(0,499)
plt.tight_layout()
#plt.savefig('scripts/Figures/angular_momentum_GH_colpa.png', dpi=600)
plt.show()
plt.clf()




#plt.scatter(x, [np.sum(L) for L in angular_momentum_data])
#plt.show()