import matplotlib.pyplot as plt
import numpy as np
import re
from util import Power, Sqrt


def wP(AA, BB, CC, DD):
    return Sqrt(Power(AA, 2) + Power(BB, 2) + Sqrt(Power(Power(AA, 2) - Power(BB, 2), 2) + 16*AA*BB*CC*DD))/Sqrt(2)


def wM(AA, BB, CC, DD):
    return Sqrt(Power(AA, 2) + Power(BB, 2) - Sqrt(Power(Power(AA, 2) - Power(BB, 2), 2) + 16*AA*BB*CC*DD))/Sqrt(2)


# Define a function to convert tuple strings to complex numbers
def parse_complex(tuple_str):
    tuple_str = tuple_str.replace(')', '')
    tuple_str = tuple_str.replace('(', '')
    tuple_str = tuple_str.split(',')
    return complex(float(tuple_str[0]), float(tuple_str[1]))


def retrieve_data_from_file(file_path):

    rest = []
    kVectors = []

    with open(file_path, 'r') as file:
        file.readline()
        for line in file:
            line = line.strip('\n')
            parts = line.split(',')

            kVector = [float(part) for part in parts[:3]]
            kVectors.append(kVector)

            # Use regular expression to find all tuples of complex numbers
            matches = re.findall(r'\([-\d.e]+,[-\d.e]+\)', line)

            if len(matches) != 0:
                rest.append([parse_complex(z) for z in matches])
            else:
                rest.append([float(z) for z in parts[3:]])

    return [kVectors, rest]


def get_file_length(file_path):

    line_counter = 0
    with open(file_path, 'r') as file:
        file.readline()
        for line in file:
            line_counter += 1

    return line_counter


def main():

    # open relevant files:

    # magnon dispersion relation
    file_magnon = 'scripts/Data/numbersJIso.txt'
    # phonon dispersion relation
    file_phonon = 'scripts/Data/numbersPh.txt'
    # DMI-like coupling parameters
    file_DMIlike = 'scripts/Data/numbersD.txt'


    file_DMIlike = "Outputs/numbersD.txt"
    file_phonon = "Outputs/numbersPh.txt"
    file_magnon = "Outputs/numbersJIso.txt"

    data_ph = retrieve_data_from_file(file_phonon)
    data_mag = retrieve_data_from_file(file_magnon)
    data_D = retrieve_data_from_file(file_DMIlike)

    num_kVec = get_file_length(file_magnon)

    # contains the k vectors
    kVectors = np.zeros(shape=(num_kVec, 3), dtype=float)
    # contains the phonon energy dispersion relation
    ph_energy = np.zeros(shape=(num_kVec, 3), dtype=float)
    # contains the magnon energy dispersion relation
    mag_energy = np.zeros(shape=(num_kVec), dtype=float)
    # polarization vectors where the columns represent the pol vectors: e1 = [pol_vec[0][0], pol_vec[1][0], pol_vec[2][0]] and so on
    pol_vec = np.zeros(shape=(num_kVec, 3, 3), dtype=float)
    # ph-mg coupling parameters
    DMI_like_coupl = np.zeros(shape=(num_kVec, 2, 3), dtype=complex)
    DMI_like_coupl_neg_k = np.zeros(shape=(num_kVec, 2, 3), dtype=complex) # for the k values with the opposite sign as above


    assert len(data_D[0]) == len(data_ph[0]) and len(data_mag[0]) == len(
        data_ph[0]),  "The data should have the same number of k Vectors!"

    # Assign the arrays with the data given in the files:
    for idx in range(num_kVec):

        kVectors[0], kVectors[1], kVectors[2] = data_ph[0][idx]
        ph_energy[idx][0], ph_energy[idx][1], ph_energy[idx][2] = data_ph[1][idx][:3]
        mag_energy[idx] = data_mag[1][idx][0]

        pol_vec[idx][0][0] = data_ph[1][idx][3]
        pol_vec[idx][1][0] = data_ph[1][idx][4]
        pol_vec[idx][2][0] = data_ph[1][idx][5]
        pol_vec[idx][0][1] = data_ph[1][idx][6]
        pol_vec[idx][1][1] = data_ph[1][idx][7]
        pol_vec[idx][2][1] = data_ph[1][idx][8]
        pol_vec[idx][0][2] = data_ph[1][idx][9]
        pol_vec[idx][1][2] = data_ph[1][idx][10]
        pol_vec[idx][2][2] = data_ph[1][idx][11]

        DMI_like_coupl[idx][0][0] = data_D[1][idx][0]
        DMI_like_coupl[idx][0][1] = data_D[1][idx][1]
        DMI_like_coupl[idx][0][2] = data_D[1][idx][2]
        DMI_like_coupl[idx][1][0] = data_D[1][idx][3]
        DMI_like_coupl[idx][1][1] = data_D[1][idx][4]
        DMI_like_coupl[idx][1][2] = data_D[1][idx][5]

        DMI_like_coupl_neg_k[idx][0][0] = data_D[1][idx][6]
        DMI_like_coupl_neg_k[idx][0][1] = data_D[1][idx][7]
        DMI_like_coupl_neg_k[idx][0][2] = data_D[1][idx][8]
        DMI_like_coupl_neg_k[idx][1][0] = data_D[1][idx][9]
        DMI_like_coupl_neg_k[idx][1][1] = data_D[1][idx][10]
        DMI_like_coupl_neg_k[idx][1][2] = data_D[1][idx][11]


    imaginary_unit = complex(0, 1)
    
    S = 4.8
    atomic_mass = 55.845  # in Dalton
    d_aniso = 6.97 * 1E-3  # anisotropy energy in meV

    # assuming the units are given as follows:
    # - magnon and phonon eigenenergies in milli eV
    # - DMI-like coupling constants in milli eV / (a.u.)

    x_values = []
    A_coefficients, B_coefficients, C_coefficients, D_coefficients = [], [], [], []

    conjugate_check = []

    D_minus_mu_coefficients = []

    for idx in range(num_kVec):

        # correct the magnon dispersion relation so that it includes anisotropy and the max number of spin excitations S
        mag_energy[idx] = S * (mag_energy[idx] + 2*d_aniso)

        # skip all points in k space where we have negative or null phonon energies:
        if ph_energy[idx][0] <= 0.0 or ph_energy[idx][1] <= 0.0 or ph_energy[idx][2] <= 0.0:
            continue

        # scale phonon energy
        #ph_energy[idx][0] *= .45
        #ph_energy[idx][1] *= .45
        #ph_energy[idx][2] *= .45

        C_coeff_all_branches = [0, 0, 0]
        D_coeff_all_branches = [0, 0, 0]

        D_minus_mu_coeff = [0,0,0]

        conjugate_difference = []

        for branch in range(3):
            C_coeff = 0
            D_coeff = 0

            for axis in range(3):

                D_mk_plus_mu = DMI_like_coupl_neg_k[idx][0][axis] + imaginary_unit * DMI_like_coupl_neg_k[idx][1][axis]
                D_k_minus_mu = DMI_like_coupl[idx][0][axis] - imaginary_unit * DMI_like_coupl[idx][1][axis]

                C_coeff +=  2*imaginary_unit/np.sqrt(2*S) * (pol_vec[idx][axis][branch]) * 3.8636 * np.sqrt(1.0 / (2 * atomic_mass * ph_energy[idx][branch])) * D_k_minus_mu
                D_coeff += -2*imaginary_unit/np.sqrt(2*S) * (pol_vec[idx][axis][branch]) * 3.8636 * np.sqrt(1.0 / (2 * atomic_mass * ph_energy[idx][branch])) * D_mk_plus_mu

                D_minus_mu_coeff[axis] = D_k_minus_mu

            C_coeff_all_branches[branch] = C_coeff
            D_coeff_all_branches[branch] = D_coeff

            conjugate_difference.append( np.conjugate(C_coeff) - D_coeff)

        conjugate_check.append(conjugate_difference  )

        D_minus_mu_coefficients.append(D_minus_mu_coeff)

        x_values.append(idx)
        A_coefficients.append([ph_energy[idx][0], ph_energy[idx][1], ph_energy[idx][2]])
        B_coefficients.append(mag_energy[idx])
        C_coefficients.append(C_coeff_all_branches)
        D_coefficients.append(D_coeff_all_branches)

    fig, axs = plt.subplots(ncols=1, nrows=3)


    for branch in range(0, 1):
        wp = [wP(A_coefficients[idx][branch], B_coefficients[idx], C_coefficients[idx]
                 [branch], D_coefficients[idx][branch]).real for idx in range(len(x_values))]
        wm = [wM(A_coefficients[idx][branch], B_coefficients[idx], C_coefficients[idx]
                 [branch], D_coefficients[idx][branch]).real for idx in range(len(x_values))]

        axs[0].scatter(x_values, [np.abs(C_coefficients[idx][branch].imag) for idx in range(len(x_values))], s=.5, label='imag')
        axs[0].scatter(x_values, [np.abs(C_coefficients[idx][branch].real) for idx in range(len(x_values))], s=.5, label='real')
        
        # multiply the C and D coefficients: as they are complex conjugates, the product should be positive and real
        #axs[1].scatter(x_values, [(D_coefficients[idx][branch]*C_coefficients[idx][branch]).imag for idx in range(len(x_values))], s=.5, label='imag')
        #axs[1].scatter(x_values, [(D_coefficients[idx][branch]*C_coefficients[idx][branch]).real for idx in range(len(x_values))], s=.5, label='real')

        # Check weather the coefficients are complex conjugates
        #axs[1].scatter(x_values, [x[branch].real for x in conjugate_check], label='real')
        #axs[1].scatter(x_values, [x[branch].imag for x in conjugate_check],label='imag')


        axs[2].scatter(x_values, wm, s=.5, label='wM', c='b' )
        axs[2].scatter(x_values, wp, s=.5, label='wM', c='b' )

    for axis in range(3):
        axs[1].plot(x_values, [D_minus_mu_coefficients[idx][axis].real for idx in range(len(x_values))], label=str(axis+1), linestyle='--')
        axs[1].plot(x_values, [D_minus_mu_coefficients[idx][axis].imag for idx in range(len(x_values))], label=str(axis+1))

    axs[0].legend()
    axs[1].legend()

    axs[0].set_ylabel(r"$C_\mathbf{k}$")
    axs[1].set_ylabel(r"$\mathcal{D}^{-,\mu}_\mathbf{k}$")

    plt.show()

    exit()

    print(DMI_like_coupl)

    """
    # Transform the units:
    mRy_in_meV = 13.605693
    mag_data_conversion_factor = mRy_in_meV  # from mRy to meV
    ph_data_conversion_factor = 1.907  # from mRy/(a.u.)^2 to meV
    D_data_conversion_factor = 1

    data_mag[1] = [[E[0] * mag_data_conversion_factor]  for E in data_mag[1]]

    for idx in range(len(data_ph[1])):
        data_ph[1][idx][:3]  = [E*ph_data_conversion_factor for E in data_ph[1][idx][:3]]
    """

    A_coefficients, B_coefficients, C_coefficients, D_coefficients = [], [], [], []

    fig, axs = plt.subplots(ncols=1, nrows=3)

    for idx in range(len(data_ph[0])):

        kVector = data_ph[0][idx]

        # phonon eigen energies and pol vectors
        ph_eigen_energies = data_ph[1][idx][:3]
        ph_pol_vec1, ph_pol_vec2, ph_pol_vec3 = data_ph[1][idx][3:
                                                                6], data_ph[1][idx][6:9], data_ph[1][idx][9:12]
        ph_pol_vectors = [ph_pol_vec1, ph_pol_vec2, ph_pol_vec3]

        if ph_eigen_energies[0] == 0.0 or ph_eigen_energies[1] == 0.0 or ph_eigen_energies[2] == 0.0:
            continue

        # magnon eigen energies
        mag_eigen_energy = S * (data_mag[1][idx][0] + 2 * d_aniso)

        # DMI-like coupling paramters
        Dxx, Dxy, Dxz, Dyx, Dyy, Dyz = data_D[1][idx][:6]
        D_matrix = [[Dxx, Dxy, Dxz], [Dyx, Dyy, Dyz]]

        # ph_eigen_energies = [E/5 for E in ph_eigen_energies] # scale the phonon eigen energies

        C_coeff_all_branches_curr_k_vec = []
        D_coeff_all_branches_curr_k_vec = []

        for branch in range(3):

            C_coeff = 0
            D_coeff = 0

            for axis in range(3):
                C_coeff += 2 * imaginary_unit / np.sqrt(2*S) * (ph_pol_vectors[branch][axis]) * 3.8636 * np.sqrt(1.0/(
                    2*abs(ph_eigen_energies[branch]*atomic_mass))) * (D_matrix[0][axis] - imaginary_unit * D_matrix[1][axis])
                D_coeff += - 2 * imaginary_unit / np.sqrt(2*S) * (ph_pol_vectors[branch][axis]) * 3.8636 * np.sqrt(1.0/(
                    2*abs(ph_eigen_energies[branch]*atomic_mass))) * (D_matrix[0][axis] + imaginary_unit * D_matrix[1][axis])

            C_coeff_all_branches_curr_k_vec.append(C_coeff)
            D_coeff_all_branches_curr_k_vec.append(D_coeff)

        A_coefficients.append(ph_eigen_energies)
        B_coefficients.append(mag_eigen_energy)
        C_coefficients.append(C_coeff_all_branches_curr_k_vec)
        D_coefficients.append(D_coeff_all_branches_curr_k_vec)
        kVectors.append(kVector)

    x_values = np.linspace(0, 1, len(kVectors))

    # for branch in range(3):
    #    plt.scatter(x_values, [E[branch] for E in A_coefficients], s=1, color='blue')

    # plt.scatter(x_values, B_coefficients, s=1)

    for branch in range(1):

        wp = [wP(A_coefficients[idx][branch], B_coefficients[idx], C_coefficients[idx]
                 [branch], D_coefficients[idx][branch]).real for idx in range(len(x_values))]
        wm = [wM(A_coefficients[idx][branch], B_coefficients[idx], C_coefficients[idx]
                 [branch], D_coefficients[idx][branch]).real for idx in range(len(x_values))]

        # phonon_or_magnon_like = [b/a for a,b in zip(wp,wm)]
        # plt.scatter(x_values, wp, s=1, label='wP', c=phonon_or_magnon_like, cmap='plasma')
        # plt.scatter(x_values, wm, s=1, label='wM', c=phonon_or_magnon_like, cmap='plasma')

        axs[0].scatter(x_values, [(D_coefficients[idx][branch].real)
                       for idx in range(len(x_values))], s=.5)
        axs[1].scatter(x_values, [(D_coefficients[idx][branch].imag)
                       for idx in range(len(x_values))], s=.5)

        axs[2].scatter(x_values, wm, s=1, label='wM', c='red')
        axs[2].scatter(x_values, wp, s=1, label='wM', c='blue')

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
