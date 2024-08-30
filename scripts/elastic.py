import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

import scienceplots
plt.style.use('science')

# Define a function to convert tuple strings to complex numbers
def parse_complex(tuple_str):
    tuple_str = tuple_str.replace(')', '')
    tuple_str = tuple_str.replace('(', '')
    tuple_str = tuple_str.split(',')
    return complex(float(tuple_str[0]), float(tuple_str[1]))

magnon_file = 'Outputs/GHxyz/mag.txt'
phonon_file = 'Outputs/GHxyz/GHz_formatted.txt'

ks = []
Js = []
omegas = []

with open(magnon_file, 'r') as file:
    file.readline()
    for line in file:
        line = line.strip('\n')
        parts = line.split(',')

        k = [float(parts[2])]
        J = float(parts[3].replace('(', '').replace(')', ''))
        Js.append(J)
        ks.append(k)


with open(phonon_file, 'r') as file:
    for line in file:
        line = line.strip('\n')
        parts = line.split(' ')

        omega_t = float(parts[3])
        omegas.append(omega_t)




Js = np.array(Js[498:])
ks = np.array(ks[498:])
omegas = np.array(omegas[499:])

def func(w):
    w_spin,w_t,sigma = 1,1,1
    return (w-w_spin)*(w*w-w_t*w_t)-sigma**2 *.5*w_t


def func_k_minus(w_spin, w_t, sigma):
    return lambda w: (w-w_spin)*(w*w-w_t*w_t)-sigma**2*.5*w_t


def func_k_plus(w_spin, w_t, sigma):
    return lambda w: (w+w_spin)*(w*w-w_t*w_t)+sigma**2*.5*w_t

def wa(w_spin, w_t, sigma):
    return (w_t+w_spin)/2+ .5*((w_t-w_spin)**2 + sigma**2)**.5


def wb(w_spin, w_t, sigma):
    return (w_t+w_spin)/2- .5*((w_t-w_spin)**2 + sigma**2)**.5


# Iron parameters

hbar = 6.6E-34/(2*np.pi)
meV_in_J = 1.602*1E-22

b2 = 7.62 * 1E6 # in J /m^3  from https://www.sciencedirect.com/science/article/pii/S0039602800001047
M = 1.75 * 1E6 # A/m from https://link.springer.com/article/10.1140/epjp/s13360-020-00294-y
rho = 7.874 * 1E-3 * 1E6 # g/cm^3 to kg/m^3
v_t = 2969 # * 9 # cheating factor
D =  155.3 * 1.602*1E-22 * (1E-10)**2 / hbar  # * 50 # cheating factor in meV A^2 from https://link.springer.com/article/10.1140/epjp/s13360-020-00294-y
gamma = 8.681 * 1E6 # 1/(sT) 
lattice_constant = 2.856 * 1E-10
H0 = 0* 1000.0/(4*np.pi) * .27 * 1E2 # 2 kOe to A/m

plot_res_min = []
plot_res_plu = []

min_dist_of_wk_and_ws = 1E8
min_sigma = 1E8
min_dist_kvec = 1E8

k_vectors = 2*np.pi/lattice_constant * np.linspace(0,.1,50)

k_vectors = np.append(k_vectors,np.linspace(1.2534E9, 1.2537E9, 1000))



for kvec in k_vectors:
    wk = kvec*v_t
    ws =  D * kvec**2 + gamma * H0
    sigma = b2 * np.sqrt(2*gamma*kvec/(rho*v_t*M)) *1

    dist = np.abs(wk-ws)
    if dist < min_dist_of_wk_and_ws and kvec > np.max(k_vectors[:10]):

        min_dist_of_wk_and_ws = dist
        min_sigma = sigma
        min_dist_kvec = kvec

    roots_plus = fsolve(func_k_minus(w_spin=ws,w_t=wk,sigma=sigma), [wk,wk,ws])
    roots_minus = fsolve(func_k_plus(w_spin=ws,w_t=wk,sigma=sigma), [wk,wk,ws])

    #plt.scatter(kvec,wa(ws,wk,0),label='w_a',color='black')
    #plt.scatter(kvec,wb(ws,wk,0),label='w_b',color='black')

    for root in roots_plus:
        #plt.scatter(kvec,root,label='w_plus',color='red')
        plot_res_plu.append(root)
        continue


    for root in roots_minus:
        #plt.scatter(kvec,root,label='w_minus',color='blue')
        plot_res_min.append(root)

        continue

print(min_sigma*hbar/meV_in_J)

#plt.axvline(x=min_dist_kvec,color='black',linestyle='--')

fig = plt.figure(figsize=(8/2.52, 6/2.52))

plt.plot(ks/(lattice_constant),meV_in_J*Js/ hbar, color='red', lw=3, alpha=0.5)
plt.plot(ks/(lattice_constant),meV_in_J*omegas/ hbar, color='blue', lw=3, alpha=0.5)

for i in range(3):
    plt.scatter(k_vectors, plot_res_min[i::3], color='black', s=1)
    plt.scatter(k_vectors, plot_res_plu[i::3], color='black', s=1)
    continue



plt.ylim(0,np.max(k_vectors)*v_t)
plt.xlim(0,np.max(k_vectors))
#plt.plot(k[:,2],omegas[:,0],label='1')
#plt.plot(k[:,2],omegas[:,1],label='2')
#plt.plot(k[:,2],omegas[:,2],label='3')
#plt.plot(k[:,2],J,label='J')

# axis in meV
yticks = np.linspace(0,np.max(plot_res_min),5)
yticklabels = [f'{t*hbar/meV_in_J:.2f}' for t in yticks]
plt.yticks(yticks, yticklabels)


xticks = np.linspace(0,np.max(k_vectors),5)
xticklabels = [f'{t*lattice_constant/(2*np.pi):.3f}' for t in xticks]

plt.xticks(xticks, xticklabels)
plt.xlabel(r'$k$ ($2\pi /a$)')
plt.ylabel("dispersion relation (meV)")
plt.tight_layout()
plt.savefig("scripts/Figures/magnetoelastic1.png", dpi=500)
plt.show()


fig = plt.figure(figsize=(8/2.52, 6/2.52))


for i in range(3):
    plt.scatter(k_vectors, plot_res_min[i::3], color='black', s=.1)
    plt.scatter(k_vectors, plot_res_plu[i::3], color='black', s=.1)
    continue



plt.xlim(1.2534E9, 1.2537E9)
plt.ylim(3.721E12,3.7225E12)


# axis in meV
yticks = np.linspace(3.721E12,3.7225E12,5)
yticklabels = [f'{t*hbar/meV_in_J:.4f}' for t in yticks]
plt.yticks(yticks, yticklabels)

xticks = np.linspace(1.2534E9, 1.2537E9,5)
xticklabels = [f'{t*lattice_constant/(2*np.pi):.3f}' for t in xticks]




plt.xticks(xticks, xticklabels)
plt.xlabel(r'$k$ ($2\pi /a$)')
plt.ylabel("dispersion relation (meV)")
plt.tight_layout()
plt.savefig("scripts/Figures/magnetoelastic2.png", dpi=500)

plt.show()