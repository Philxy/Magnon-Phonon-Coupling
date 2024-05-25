import numpy as np
import matplotlib.pyplot as plt
import scienceplots


e1 = np.array([1, 0, 0])
e2 = np.array([0, 1, 0])
e3 = np.array([0, 0, 1])


phiX = np.pi/2
phiY = -np.pi/2
phiZ = 0
epX = 1
epY = 1
time = np.linspace(0, 2*np.pi, 100)
x = [epX * (np.exp(1j*(t+phiX)) + np.exp(-1j*(t+phiX))) for t in time]
y = [epY * (np.exp(1j*(t+phiY)) + np.exp(-1j*(t+phiY))) for t in time]


def u(t, phiX=phiX, phiY=phiY, phiZ=phiZ, e1=e1, e2=e2, e3=e3, A1=1, A2=1):
    pol_vecs = [e1, e2, e3]
    result = np.zeros(3)
    for e in pol_vecs:
        result += np.real((A1 * np.exp(1j*(t+phiX)) +
                          A2 * np.exp(-1j*(t+phiY))) * e)

    return result



def u_simple(t, delta_phi):
    e = 1/np.sqrt(2)*np.array([np.exp(delta_phi *1.0j), 1, 0])
    return e*np.exp(1j*(t))

plt.style.use('science')

fig, ax = plt.subplots(figsize=(7/2.52, 8.3/2.52))


dphi = -1.5726



ut = [u_simple(t, dphi) for t in time]
utx = [u[0] for u in ut]
uty = [u[1] for u in ut]
ax.plot(np.real(utx), np.real(uty), label=r'$\Gamma-P=(2\pi/a,2\pi/a,2\pi/a)$', linewidth=3)


ut = [u_simple(t, np.pi/2) for t in time]
utx = [u[0] for u in ut]
uty = [u[1] for u in ut]
ax.plot(np.real(utx), np.real(uty), label=r'$\Gamma-H=(0,0,2\pi/a)$', linestyle='--', linewidth=3)

ax.set_xticks([-1/np.sqrt(2), 0, 1/np.sqrt(2)], labels=[r'$-\frac{1}{\sqrt{2}}$', '0', r'$\frac{1}{\sqrt{2}}$'])
ax.set_yticks([-1/np.sqrt(2), 0, 1/np.sqrt(2)], labels=[r'$-\frac{1}{\sqrt{2}}$', '0', r'$\frac{1}{\sqrt{2}}$'])

ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.3)) 
ax.set_xlabel(r'displacement $u_x$')
ax.set_ylabel(r'displacement $u_y$')
plt.tight_layout()
plt.show()


fig, ax = plt.subplots(figsize=(6, 6))


phiX = 0
phiY = 0
phiZ = 0
ut = [u(t, phiX=phiX, phiY=phiY) for t in time]
utx = [u[0] for u in ut]
uty = [u[1] for u in ut]

ax.plot(np.real(utx), np.real(uty))


A1 = 1/np.sqrt(2)*(1+1.0j)
A2 = 1/np.sqrt(2)*(1-1.0j)
phiX = 0
phiY = 0
phiZ = 0
ut = [u(t, phiX=phiX, phiY=phiY, A1=A1, A2=A2) for t in time]
utx = [u[0] for u in ut]
uty = [u[1] for u in ut]

ax.plot(np.real(utx), np.real(uty))


ax.set_xlabel('x')
ax.set_ylabel('y')
plt.show()
