import numpy as np
import matplotlib.pyplot as plt


def deltaDistrApprox(x, a=0.1):
    return 1.0 / (a * np.sqrt(np.pi)) * np.exp(-(x / a) * (x / a))


x = np.linspace(-3, 3, 1000)

for a in [0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0]:
    y = deltaDistrApprox(x, a=a)

    plt.plot(x, y, label=f'a={a}')
plt.legend()
plt.show()
