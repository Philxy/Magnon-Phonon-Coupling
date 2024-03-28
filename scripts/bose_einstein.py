import matplotlib.pyplot as plt
import numpy as np
import scipy as sc


x = np.linspace(-10, 10, 100)


def bose_einstein(x, E, mu):
    return 1 / (np.exp((x-mu) / E) - 1)


def gaussian(x):
    return 1/np.sqrt(np.pi) * np.exp(-(x)**2 / 0.1)


#for E in np.logspace(0.1, 10, 3):
#    plt.plot(x, bose_einstein(x, E, 0))


plt.plot(x, gaussian(x))

#plt.yscale('log')
plt.show()
