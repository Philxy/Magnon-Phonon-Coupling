from util import Power, Sqrt, gaussian
import cmath
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm


a = 10E-10

# phonon disp relation


def omega(k, v):
    return np.cos(k*a/2.0)


# magnon disp relation
def OMEGA(k):
    return Sqrt(np.sin(k*a/2.0)**2)


def AA(k):
    return omega(k, 0)


def BB(k):
    return OMEGA(k)


def CC(k):
    # return 0
    return 0.1*gaussian(k*a/2.0, 0, 0.5)


def DD(k):
    return 0.1
    return 0.1*gaussian(k*a/2.0, 0, 0.5)


def wP(k):
    return Sqrt(Power(AA(k), 2) + Power(BB(k), 2) + Sqrt(Power(Power(AA(k), 2) - Power(BB(k), 2), 2) + 16*AA(k)*BB(k)*CC(k)*DD(k)))/Sqrt(2)


def wM(k):
    return Sqrt(Power(AA(k), 2) + Power(BB(k), 2) - Sqrt(Power(Power(AA(k), 2) - Power(BB(k), 2), 2) + 16*AA(k)*BB(k)*CC(k)*DD(k)))/Sqrt(2)


k_values = np.linspace(-np.pi/a, np.pi/a, 1000)

energy_all = [wP(k) + wM(k) for k in k_values]
energy_m = [wM(k) for k in k_values]
energy_p = [wP(k) for k in k_values]

# plt.plot(k_values,energy_all, label='all')
plt.plot(k_values, energy_m, label='$\omega_{-}$', linewidth=3)
plt.plot(k_values, energy_p, label='$\omega_{+}$', linewidth=3)
plt.plot(k_values, [AA(k) for k in k_values], linestyle='dashed')
plt.plot(k_values, [BB(k) for k in k_values], linestyle='dashed')

# plt.plot(k_values, [CC(k) for k in k_values])

plt.legend()
plt.show()
