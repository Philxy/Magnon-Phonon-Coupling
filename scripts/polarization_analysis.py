import matplotlib.pyplot as plt
import numpy as np
import scienceplots

plt.style.use('science')

pol_vec_file_path = "Outputs/GP_path/polVec.txt"


def parse_complex(tuple_str):
    tuple_str = tuple_str.replace(')', '')
    tuple_str = tuple_str.replace('(', '')
    tuple_str = tuple_str.split(',')
    return complex(float(tuple_str[0]), float(tuple_str[1]))


file_content_list_parsed = []


with open(pol_vec_file_path, "r") as file:
    file.readline()
    for line in file:
        line = line.strip('\n')
        line = line.split(' ')

        line_parsed = [parse_complex(z) for z in line if z != '']

        file_content_list_parsed.append(line_parsed)



pol_vec_for_exemplary_pol_vec = file_content_list_parsed[190]


pol_vectors = []

for branch in range(8):
    pvec = np.zeros(3, dtype=complex)
    for i in range(3):
        pvec[i] = pol_vec_for_exemplary_pol_vec[branch*3+i]

    if np.abs(np.linalg.norm(pvec)- 1)  < 0.001:
        pol_vectors.append(pvec)



print(pol_vectors)

jones_vector1,jones_vector2,jones_vector3 = pol_vectors
time = np.linspace(0, 2*np.pi, 1000)

u1x = [(jones_vector1[0] * np.exp(1.0j*t)).real for t in time]
u1y = [(jones_vector1[1] * np.exp(1.0j*t)).real for t in time]
u1z = [(jones_vector1[2] * np.exp(1.0j*t)).real for t in time]

u2x = [(jones_vector2[0] * np.exp(1.0j*t)).real for t in time]
u2y = [(jones_vector2[1] * np.exp(1.0j*t)).real for t in time]
u2z = [(jones_vector2[2] * np.exp(1.0j*t)).real for t in time]

u3x = [(jones_vector3[0] * np.exp(1.0j*t)).real for t in time]
u3y = [(jones_vector3[1] * np.exp(1.0j*t)).real for t in time]
u3z = [(jones_vector3[2] * np.exp(1.0j*t)).real for t in time]

import cmath

polar_forms = [cmath.polar(np.sum(p)) for p in pol_vectors]

for i,p1 in enumerate(polar_forms):
    for j,p2 in enumerate(polar_forms):
        print(i,j,np.abs(p1[1]-p2[1])/np.pi)

print("polar forms:", polar_forms)

plt.plot(u1x, u1y, )
plt.plot(u2x, u2y, linestyle='solid',linewidth=2)
plt.plot(u3x, u3y, linestyle='dashed',linewidth=2)


plt.xlim(-1.2/np.sqrt(2), 1.2/np.sqrt(2))
plt.ylim(-1.2/np.sqrt(2), 1.2/np.sqrt(2))

ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.draw()

plt.figsize=(8/2.52, 8/2.52)
plt.xlabel(r'$x$ $\left(\sqrt{\frac{\hbar}{2m\omega_{\boldsymbol{k}\lambda}}}\right)$')
plt.ylabel(r'$y$ $\left(\sqrt{\frac{\hbar}{2m\omega_{\boldsymbol{k}\lambda}}}\right)$')
plt.yticks([-1/np.sqrt(2), 0, 1/np.sqrt(2)], labels=[r'$-\frac{1}{\sqrt{2}}$', '0', r'$\frac{1}{\sqrt{2}}$'])
plt.xticks([-1/np.sqrt(2), 0, 1/np.sqrt(2)], labels=[r'$-\frac{1}{\sqrt{2}}$', '0', r'$\frac{1}{\sqrt{2}}$'])


a = 0.5
plt.axvline(-1/np.sqrt(2), color='black', linewidth=0.5, alpha=a)
plt.axvline(1/np.sqrt(2), color='black', linewidth=0.5, alpha=a)
plt.axhline(-1/np.sqrt(2), color='black', linewidth=0.5, alpha=a)
plt.axhline(1/np.sqrt(2), color='black', linewidth=0.5, alpha=a)


plt.tight_layout()
plt.savefig('scripts/Figures/polVec_GP.png', dpi=600)
plt.show()


# phase shift

pol_vectors = []


for i in range(len(file_content_list_parsed)):

    pol_vec_for_exemplary_pol_vec = file_content_list_parsed[i]

    for branch in range(8):
        pvec = np.zeros(3, dtype=complex)
        for i in range(3):
            pvec[i] = pol_vec_for_exemplary_pol_vec[branch*3+i]

        if np.abs(np.linalg.norm(pvec)- 1)  < 0.001:
            pol_vectors.append(pvec)


phase_shifts = []


for pol_vec in pol_vectors:

    components_in_polar_form = [cmath.polar(component) for component in pol_vec]

    print("----------------")
    print("Phase shift 1 2:", np.abs(components_in_polar_form[0][1] - components_in_polar_form[1][1])/np.pi)
    print("Phase shift 2 3:", np.abs(components_in_polar_form[1][1] - components_in_polar_form[2][1])/np.pi)
    print("Phase shift: 1 3:", np.abs(components_in_polar_form[2][1] - components_in_polar_form[0][1])/np.pi)

    phase_shifts.append(np.abs(components_in_polar_form[0][1] - components_in_polar_form[1][1])/np.pi)
    #phase_shifts.append(np.abs(components_in_polar_form[1][1] - components_in_polar_form[2][1])/np.pi)
    #phase_shifts.append(np.abs(components_in_polar_form[2][1] - components_in_polar_form[0][1])/np.pi)

plt.scatter(range(len(phase_shifts)), phase_shifts)
plt.show()

'''
AB_file_path = "Outputs/GH_path/AB.txt"

data_AB = []

with open(AB_file_path, "r") as file:

    file.readline()
    for line in file:

        line = line.strip('\n')
        line = line.split(' ')
        line_parsed = [parse_complex(a) for a in line]

        data_AB.append(line_parsed)

print(np.abs(data_AB[400][0])/np.abs(data_AB[400][1]))
print(np.abs(data_AB[100][0])/np.abs(data_AB[100][1]))
print(np.abs(data_AB[800][1])/np.abs(data_AB[800][0]))


for i, ab_pair in enumerate(data_AB):
    
    plt.scatter(i, np.abs(ab_pair[1]), color='black')
    plt.scatter(i, np.abs(ab_pair[0]), color='red')



plt.show()
'''