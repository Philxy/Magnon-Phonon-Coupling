import numpy as np
import matplotlib.pyplot as plt


filename = 'Outputs/time_evolut_energies.txt'

with open(filename, 'r') as file:
    data = file.readlines()


total_energy = []
magnetic_energy = []
phonon_energy = []
time = []

for line in data:
    values = line.split()
    # Process the values as needed
    t = float(values[0])
    magE = float(values[1])
    phE = float(values[2])
    totE = float(values[3])

    total_energy.append(totE)
    magnetic_energy.append(magE)
    phonon_energy.append(phE)
    time.append(t)


plt.plot(time, total_energy, label='Total Energy')
plt.plot(time, magnetic_energy, label='Magnetic Energy')
plt.plot(time, phonon_energy, label='Phonon Energy')
plt.legend()
plt.show()
