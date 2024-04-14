import numpy as np
import matplotlib.pyplot as plt
import regex as re
from util import Sqrt, Power
import scienceplots


def omega_p(A,B,C,D):
    return np.sqrt(np.power(A,2) + np.power(B,2) + np.sqrt(np.power(np.power(A,2) - np.power(B,2),2) + 16*A*B*(C*D).real))/np.sqrt(2)


def omega_m(A,B,C,D):
    return np.sqrt(np.power(A,2) + np.power(B,2) - np.sqrt(np.power(np.power(A,2) - np.power(B,2),2) + 16*A*B*(C*D).real))/np.sqrt(2)


def wP(A,B,C,D):
    return omega_p(A,B,C,D)


def wM(A,B,C,D):
    return omega_m(A,B,C,D)

def List(*x):
    return [*x]

def Smat(AA,BB,CC,DD):


    kk = Sqrt(Power(Power(AA,2) - Power(BB,2),2) + 16*AA*BB*CC*DD)
    wPlus = wP(AA,BB,CC,DD)
    wMinus = wM(AA,BB,CC,DD)

    return List(List((kk + (AA + BB)*(-AA - BB + 2*wMinus))/(4.*AA*CC),(-kk + (AA + BB)*(-AA - BB + 2*wPlus))/(4.*AA*CC),(kk - (AA + BB)*(AA + BB + 2*wMinus))/(4.*AA*CC),
    -0.25*(kk + (AA + BB)*(AA + BB + 2*wPlus))/(AA*CC)),List((-8*AA*CC*DD + Power(AA,2)*(2*BB - 2*wMinus) + (-Power(BB,2) + kk)*(2*BB - 2*wMinus))/(8.*AA*Power(CC,2)),
    (-8*AA*CC*DD + Power(AA,2)*(2*BB - 2*wPlus) + (Power(BB,2) + kk)*(-2*BB + 2*wPlus))/(8.*AA*Power(CC,2)),
    (-8*AA*CC*DD + Power(AA,2)*(2*BB + 2*wMinus) + (-Power(BB,2) + kk)*(2*BB + 2*wMinus))/(8.*AA*Power(CC,2)),
    (-8*AA*CC*DD + Power(AA,2)*(2*BB + 2*wPlus) - (Power(BB,2) + kk)*(2*BB + 2*wPlus))/(8.*AA*Power(CC,2))),
   List((-kk + (AA - BB)*(AA - BB + 2*wMinus))/(4.*AA*CC),(kk + (AA - BB)*(AA - BB + 2*wPlus))/(4.*AA*CC),(-kk + (AA - BB)*(AA - BB - 2*wMinus))/(4.*AA*CC),
    (kk + (AA - BB)*(AA - BB - 2*wPlus))/(4.*AA*CC)),List(1,1,1,1))


def Sinvmat(AA,BB,CC,DD):
    kk = Sqrt(Power(Power(AA,2) - Power(BB,2),2) + 16*AA*BB*CC*DD)
    wPlus = wP(AA,BB,CC,DD)
    wMinus = wM(AA,BB,CC,DD)

    return List(List(-0.25*(CC*(Power(AA,2) + Power(BB,2) - kk + 2*BB*wMinus - 2*AA*(BB + wMinus)))/(kk*wMinus),-((AA*Power(CC,2))/(kk*wMinus)),
    (CC*(kk - (AA + BB)*(AA + BB + 2*wMinus)))/(4.*kk*wMinus),(-4*AA*CC*DD + Power(AA,2)*(BB + wMinus) - (Power(BB,2) - kk)*(BB + wMinus))/(4.*kk*wMinus)),
   List((CC*(Power(AA,2) + Power(BB,2) + kk + 2*BB*wPlus - 2*AA*(BB + wPlus)))/(4.*kk*wPlus),(AA*Power(CC,2))/(kk*wPlus),
    (CC*(kk + (AA + BB)*(AA + BB + 2*wPlus)))/(4.*kk*wPlus),(4*AA*CC*DD - Power(AA,2)*(BB + wPlus) + (Power(BB,2) + kk)*(BB + wPlus))/(4.*kk*wPlus)),
   List((CC*(-kk + (AA - BB)*(AA - BB + 2*wMinus)))/(4.*kk*wMinus),(AA*Power(CC,2))/(kk*wMinus),(CC*(-kk + (AA + BB)*(AA + BB - 2*wMinus)))/(4.*kk*wMinus),
    (4*AA*CC*DD + (Power(BB,2) - kk)*(BB - wMinus) + Power(AA,2)*(-BB + wMinus))/(4.*kk*wMinus)),
   List(-0.25*(CC*(kk + (AA - BB)*(AA - BB + 2*wPlus)))/(kk*wPlus),-((AA*Power(CC,2))/(kk*wPlus)),(CC*(-kk - (AA + BB)*(AA + BB - 2*wPlus)))/(4.*kk*wPlus),
    (-4*AA*CC*DD + Power(AA,2)*(BB - wPlus) - (Power(BB,2) + kk)*(BB - wPlus))/(4.*kk*wPlus)))



# Define a function to convert tuple strings to complex numbers
def parse_complex(tuple_str):
    tuple_str = tuple_str.replace(')', '')
    tuple_str = tuple_str.replace('(', '')
    tuple_str = tuple_str.split(',')
    return complex(float(tuple_str[0]), float(tuple_str[1]))


def retrieve_data_from_file(file_path, separator=','):

    rest = []

    with open(file_path, 'r') as file:
        file.readline()
        for line in file:
            line = line.strip('\n')
            parts = line.split(separator)

            # Use regular expression to find all tuples of complex numbers
            matches = re.findall(r'\([-\d.e]+,[-\d.e]+\)', line)

            if len(matches) != 0:
                rest.append([parse_complex(z) for z in matches])
            else:
                rest.append([float(z) for z in parts])

    return rest


filename_eigen = "scripts/Data/analytical/Eigenenergies.txt"
filename_ph = "scripts/Data/analytical/path.txt"
filename_mag = "scripts/Data/analytical/magDisp.txt"

eigenenergies = retrieve_data_from_file(filename_eigen)


CDs = retrieve_data_from_file("scripts/Data/analytical/CD.txt")
phDisp = retrieve_data_from_file(filename_ph, separator=' ')
magDisp = retrieve_data_from_file(filename_mag)
eigenenergies = retrieve_data_from_file(filename_eigen)
S_matrix = retrieve_data_from_file("scripts/Data/analytical/EVec.txt")

mag, ph1, ph2, ph3, eig, C1, C2, C3, D1, D2, D3 = [], [], [], [], [], [], [], [], [], [], []


for i in range(len(eigenenergies)):
    ph1.append(phDisp[i][3])
    ph2.append(phDisp[i][4])
    ph3.append(phDisp[i][5])

    C1.append(CDs[i][0])
    C2.append(CDs[i][1])
    C3.append(CDs[i][2])

    D1.append(CDs[i][3])
    D2.append(CDs[i][4])
    D3.append(CDs[i][5])

    mag.append(magDisp[i][0].real)


modes_ph = [ph1, ph2, ph3]
modes_C = [C1, C2, C3]
modes_D = [D1, D2, D3] 

x = np.arange(0, len(mag), 1)



#  calc the transition matrix elements 
y_S_mat = [[],[],[]]
y_S_inv_mat = [[],[],[]]

for idx in range(len(mag)):
    for branch in range(3):
        A = mag[idx]
        B = modes_ph[branch][idx]
        C = modes_C[branch][idx]
        D = modes_D[branch][idx]

        s_mat = Smat(A,B,C,D)    
        s_inv_mat = Sinvmat(A,B,C,D)

        y_S_mat[branch].append(s_mat)
        y_S_inv_mat[branch].append(s_inv_mat)


def mid_of_two(list1, list2):
    return [np.max([list1[i], list2[i]]) - (np.abs(list1[i] - list2[i])/2.0) for i in range(len(list1))]

y_trans1_mode1 = [y_S_inv_mat[0][i][1][3] for i in range(len(mag))]
y_trans1_mode2 = [y_S_inv_mat[1][i][1][3] for i in range(len(mag))]

y_trans2_mode1 = [y_S_inv_mat[0][i][0][3] for i in range(len(mag))]
y_trans2_mode2 = [y_S_inv_mat[1][i][0][3] for i in range(len(mag))]

y_S_inv_long1 = [y_S_inv_mat[2][i][1][3] for i in range(len(mag))]
y_S_inv_long2 = [y_S_inv_mat[2][i][0][3] for i in range(len(mag))]

y_S_inv1_mat = mid_of_two(y_trans1_mode1, y_trans1_mode2)
y_S_inv2_mat = mid_of_two(y_trans2_mode1, y_trans2_mode2)


# plotting the transition matrix elements
plt.style.use('science')
fig, ax = plt.figure(figsize=(16/2.52, 7/2.52)), plt.gca()


from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
# Define colormap and normalizer
cmap = plt.get_cmap('rainbow_r')
norm = Normalize(vmin=0, vmax=1)
mappable = ScalarMappable(norm=norm, cmap=cmap)
# Define colormap and normalizer
cmap = plt.get_cmap('rainbow_r')
norm = Normalize(vmin=0, vmax=1)
mappable = ScalarMappable(norm=norm, cmap=cmap)

# Plotting your data (make sure x, y_S_inv1_mat, y_S_inv2_mat are defined)
ax.plot(x[10:], y_S_inv1_mat[10:], color = 'black')
ax.plot(x[10:], y_S_inv2_mat[10:], color = 'black')

# Set tick labels with proper LaTeX formatting
ax.set_xticks([10, 1000])
ax.set_xticklabels([r'$\Gamma$', 'H'])

# Add colorbar
cbar = fig.colorbar(mappable, ax=ax, orientation='horizontal', location='top', pad=0.06, anchor=(0.5, 0.7))

cbar.set_ticks([0, 1])
cbar.set_ticklabels(['magnon-like', 'phonon-like'])
plt.tight_layout()

ax.set_ylabel(r'transformation matrix elements')

plt.savefig('scripts/Figures/analytical/transition_matrix_elements_ana.png', dpi=600)
plt.show()

plt.clf()


# plotting the dispersion relation

fig, axs = plt.subplots(ncols=2,nrows=1,figsize=(16/2.52,32/2.52), sharey=True)
plt.style.use('science')


# transversal modes

y_minus_mode1 = [omega_m(modes_ph[0][i], mag[i], modes_C[0][i], modes_D[0][i]) for i in range(len(mag))]
y_plus_mode1 = [omega_p(modes_ph[0][i], mag[i], modes_C[0][i], modes_D[0][i])  for i in range(len(mag))]

y_minus_mode2 = [omega_m(modes_ph[1][i], mag[i], modes_C[1][i], modes_D[1][i]) for i in range(len(mag))]
y_plus_mode2 = [omega_p(modes_ph[1][i], mag[i], modes_C[1][i], modes_D[1][i])  for i in range(len(mag))]

#y_transversal_plus = [np.max([y_plus_mode1[i], y_plus_mode2[i]]) - (np.abs(y_plus_mode1[i] - y_plus_mode2[i])/2.0) for i in range(len(mag))]
#y_transversal_minus = [np.max([y_minus_mode1[i], y_minus_mode2[i]]) - (np.abs(y_minus_mode1[i] - y_minus_mode2[i])/2.0) for i in range(len(mag))]


y_transversal_plus = [complex(mid_of_two(y_plus_mode1, y_plus_mode2)[i]).real for i in range(len(mag))]
y_transversal_minus = [complex(mid_of_two(y_minus_mode1, y_minus_mode2)[i]).real for i in range(len(mag))]

y_longitudinal_plus = [omega_p(modes_ph[2][i], mag[i], modes_C[2][i], modes_D[2][i]) for i in range(len(mag))]
y_longitudinal_minus = [omega_m(modes_ph[2][i], mag[i], modes_C[2][i], modes_D[2][i]) for i in range(len(mag))]



cmap = plt.get_cmap('rainbow_r')
cmap.norm = Normalize(vmin=0, vmax=1)


# Create inset axis
ax1_inset = axs[0].inset_axes([0.1, -1, 0.85, 0.85])
ax2_inset = axs[1].inset_axes([0.1, -1, 0.85, 0.85])



for i in range(len(mag)-1):

    c1 = np.abs(complex(y_S_inv1_mat[i+1]).real)
    c2 = np.abs(complex(y_S_inv2_mat[i+1]).real)
    c3 = np.abs(complex(y_S_inv_long1[i+1]).real)
    c4 = np.abs(complex(y_S_inv_long2[i+1]).real)

    # transversal 
    axs[0].plot([x[i], x[i+1]], [y_transversal_plus[i], y_transversal_plus[i+1]], color=cmap(c1), linewidth=1.8)
    axs[0].plot([x[i], x[i+1]], [y_transversal_minus[i], y_transversal_minus[i+1]], color=cmap(c2), linewidth=1.8)

    ax1_inset.plot([x[i], x[i+1]], [y_transversal_plus[i], y_transversal_plus[i+1]], color=cmap(c1), linewidth=3)
    ax1_inset.plot([x[i], x[i+1]], [y_transversal_minus[i], y_transversal_minus[i+1]], color=cmap(c2), linewidth=3)

    # longitudinal
    axs[1].plot([x[i], x[i+1]], [y_longitudinal_plus[i], y_longitudinal_plus[i+1]], color=cmap(c3), linewidth=1.8)
    axs[1].plot([x[i], x[i+1]], [y_longitudinal_minus[i], y_longitudinal_minus[i+1]], color=cmap(c4), linewidth=1.8)

    ax2_inset.plot([x[i], x[i+1]], [y_longitudinal_plus[i], y_longitudinal_plus[i+1]], color=cmap(c3), linewidth=3)
    ax2_inset.plot([x[i], x[i+1]], [y_longitudinal_minus[i], y_longitudinal_minus[i+1]], color=cmap(c4), linewidth=3)


axs[0].set_title('transversal modes')
axs[1].set_title('longitudinal mode')

ax1_inset.set_xlim(40,80)
ax2_inset.set_xlim(40,150)
ax1_inset.set_ylim(1.75,3.75)
ax2_inset.set_ylim(2,12)

axs[0].indicate_inset_zoom(ax1_inset, edgecolor="black")
axs[1].indicate_inset_zoom(ax2_inset, edgecolor="black")

axs[0].set_ylabel(r'dispersion relation (meV)')

axs[0].set_xticks([0,1000])
axs[1].set_xticks([0,1000])
axs[0].set_xticklabels(['$\Gamma$', 'H'])
axs[1].set_xticklabels(['$\Gamma$', 'H'])

ax1_inset.set_xticks([])
ax2_inset.set_xticks([])

plt.tight_layout()
plt.savefig('scripts/Figures/analytical/dispersion_relation_ana.png', dpi=600)
plt.show()
plt.clf()



