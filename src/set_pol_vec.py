import numpy as np
import matplotlib.pyplot as plt

# Goal: Given a path in the bcc Brillouin zone: G-H-N-G-P-H


num_kpoints_per_segment = 1000
num_paths = 5

file_name = 'scripts/Data/set_pol_vec_input_output/ph_data.txt'


omega = np.zeros((num_paths * num_kpoints_per_segment, 3))
pol_vec = np.zeros((num_paths * num_kpoints_per_segment, 3, 3))
kvec = np.zeros((num_paths * num_kpoints_per_segment, 3))


with open(file_name, 'r') as file:
    file.readline()


    for i, line in enumerate(file):

        if i >= num_kpoints_per_segment*num_paths:
            break

        line = line.strip('\n')
        line = line.split(' ')

        kvec[i] = np.array([float(line[0]), float(line[1]), float(line[2])])
        omega[i] = np.array([float(line[3]), float(line[4]), float(line[5])])
        pol_vec_1 = np.array([float(line[6]), float(line[7]), float(line[8])])
        pol_vec_2 = np.array([float(line[9]), float(line[10]), float(line[11])])
        pol_vec_3 = np.array([float(line[12]), float(line[13]), float(line[14])])


        pol_vec[i] = np.array([pol_vec_1, pol_vec_2, pol_vec_3])



 
e3_GH = np.array([ 0.0, 1.0, 0.0])
e2_GH = np.array([1, 0.0, 0])
e1_GH = np.array([0.0, 0.0, 1])


e3_HN = np.array([-0.707107, -0.707107, 0.0])
e2_HN = np.array([0.0, 0.0, 1.0 ])
e1_HN = np.array([-0.707107, 0.707107, 0.0 ])


e3_NG = np.array([-0.707107, -0.707107, 0.0])
e2_NG = np.array([0.0 ,0.0, 1.0 ])
e1_NG = np.array([-0.707107 ,0.707107, 0.0])


e3_GP = np.array([ 0.57735, 0.57735, 0.57735]) # longitidunal
e2_GP = np.array([ 0.816495, -0.40669, -0.409804])
e1_GP = np.array([0.001798, -0.708004, 0.706206 ])

e3_PH =   np.array([-0.829181, -0.558981, 0.0])
e2_PH =   np.array([-0.558981, 0.829181, 0.0 ])
e1_PH =   np.array([0.0, 0.0, 1.0 ])


print(np.linalg.norm(e1_GH))
print(np.linalg.norm(e2_GH))
print(np.linalg.norm(e3_GH))
print(np.linalg.norm(e1_GP))
print(np.linalg.norm(e2_GP))
print(np.linalg.norm(e3_GP))
print(np.linalg.norm(e1_PH))
print(np.linalg.norm(e2_PH))
print(np.linalg.norm(e3_PH))


omegas_output = np.zeros((num_paths * num_kpoints_per_segment, 3))
pol_vec_output = np.zeros((num_paths * num_kpoints_per_segment, 3, 3))
kvec_output = np.zeros((num_paths * num_kpoints_per_segment, 3))

dict_olt_idx_to_new_idx = {0: 2, 1: 0, 2: 1}


for i in range(num_paths * num_kpoints_per_segment):
    
    if i < num_kpoints_per_segment:
        # G-H

        #omega_trans = np.min([omega[i][0], omega[i][1], omega[i][2]])
        #omega_long = np.max([omega[i][0], omega[i][1], omega[i][2]])
        omegas_output[i] = omega[i]

        kvec_output[i][0] = kvec[i][dict_olt_idx_to_new_idx[0]]
        kvec_output[i][1] = kvec[i][dict_olt_idx_to_new_idx[1]] 
        kvec_output[i][2] = kvec[i][dict_olt_idx_to_new_idx[2]]


        pol_vec_output[i][0][0] = e1_GH[dict_olt_idx_to_new_idx[0]]
        pol_vec_output[i][1][0] = e1_GH[dict_olt_idx_to_new_idx[1]]
        pol_vec_output[i][2][0] = e1_GH[dict_olt_idx_to_new_idx[2]]
        pol_vec_output[i][0][1] = e2_GH[dict_olt_idx_to_new_idx[0]]
        pol_vec_output[i][1][1] = e2_GH[dict_olt_idx_to_new_idx[1]]
        pol_vec_output[i][2][1] = e2_GH[dict_olt_idx_to_new_idx[2]]
        pol_vec_output[i][0][2] = e3_GH[dict_olt_idx_to_new_idx[0]]
        pol_vec_output[i][1][2] = e3_GH[dict_olt_idx_to_new_idx[1]]
        pol_vec_output[i][2][2] = e3_GH[dict_olt_idx_to_new_idx[2]]

    elif i < 2*num_kpoints_per_segment:
        # H-N
        
        #omegas_output[i][0] = omega[i][dict_olt_idx_to_new_idx[0]]
        #omegas_output[i][1] = omega[i][dict_olt_idx_to_new_idx[1]]
        #omegas_output[i][2] = omega[i][dict_olt_idx_to_new_idx[2]]

        omegas_output[i] = omega[i]


        kvec_output[i][0] = kvec[i][dict_olt_idx_to_new_idx[0]]
        kvec_output[i][1] = kvec[i][dict_olt_idx_to_new_idx[1]] 
        kvec_output[i][2] = kvec[i][dict_olt_idx_to_new_idx[2]]

        pol_vec_output[i][0][0] = e1_HN[dict_olt_idx_to_new_idx[0]]
        pol_vec_output[i][1][0] = e1_HN[dict_olt_idx_to_new_idx[1]]
        pol_vec_output[i][2][0] = e1_HN[dict_olt_idx_to_new_idx[2]]
        pol_vec_output[i][0][1] = e2_HN[dict_olt_idx_to_new_idx[0]]
        pol_vec_output[i][1][1] = e2_HN[dict_olt_idx_to_new_idx[1]]
        pol_vec_output[i][2][1] = e2_HN[dict_olt_idx_to_new_idx[2]]
        pol_vec_output[i][0][2] = e3_HN[dict_olt_idx_to_new_idx[0]]
        pol_vec_output[i][1][2] = e3_HN[dict_olt_idx_to_new_idx[1]]
        pol_vec_output[i][2][2] = e3_HN[dict_olt_idx_to_new_idx[2]]



    elif i < 3*num_kpoints_per_segment:
        # N-G
        #omegas_output[i][0] = omega[i][dict_olt_idx_to_new_idx[0]]
        #omegas_output[i][1] = omega[i][dict_olt_idx_to_new_idx[1]]
        #omegas_output[i][2] = omega[i][dict_olt_idx_to_new_idx[2]]
        omegas_output[i] = omega[i]


        kvec_output[i][0] = kvec[i][dict_olt_idx_to_new_idx[0]]
        kvec_output[i][1] = kvec[i][dict_olt_idx_to_new_idx[1]] 
        kvec_output[i][2] = kvec[i][dict_olt_idx_to_new_idx[2]]

        pol_vec_output[i][0][0] = e1_NG[dict_olt_idx_to_new_idx[0]]
        pol_vec_output[i][1][0] = e1_NG[dict_olt_idx_to_new_idx[1]]
        pol_vec_output[i][2][0] = e1_NG[dict_olt_idx_to_new_idx[2]]
        pol_vec_output[i][0][1] = e2_NG[dict_olt_idx_to_new_idx[0]]
        pol_vec_output[i][1][1] = e2_NG[dict_olt_idx_to_new_idx[1]]
        pol_vec_output[i][2][1] = e2_NG[dict_olt_idx_to_new_idx[2]]
        pol_vec_output[i][0][2] = e3_NG[dict_olt_idx_to_new_idx[0]]
        pol_vec_output[i][1][2] = e3_NG[dict_olt_idx_to_new_idx[1]]
        pol_vec_output[i][2][2] = e3_NG[dict_olt_idx_to_new_idx[2]]


    elif i < 4*num_kpoints_per_segment:
        # G-P
        #omega_trans = np.min([omega[i][0], omega[i][1], omega[i][2]])
        #omega_long = np.max([omega[i][0], omega[i][1], omega[i][2]])
        #omegas_output[i] = np.array([omega_trans,omega_trans,omega_long])

        omegas_output[i] = omega[i]

        assert(np.abs(omega[i][2]) >= np.abs(omega[i][0])) , f'{omega[i][2]} {omega[i][0]}'
        assert(np.abs(omega[i][2]) >= np.abs(omega[i][1])) , f'{omega[i][2]} {omega[i][1]}'

        kvec_output[i][0] = kvec[i][dict_olt_idx_to_new_idx[0]]
        kvec_output[i][1] = kvec[i][dict_olt_idx_to_new_idx[1]] 
        kvec_output[i][2] = kvec[i][dict_olt_idx_to_new_idx[2]]

        pol_vec_output[i][0][0] = e1_GP[dict_olt_idx_to_new_idx[0]]
        pol_vec_output[i][1][0] = e1_GP[dict_olt_idx_to_new_idx[1]]
        pol_vec_output[i][2][0] = e1_GP[dict_olt_idx_to_new_idx[2]]
        pol_vec_output[i][0][1] = e2_GP[dict_olt_idx_to_new_idx[0]]
        pol_vec_output[i][1][1] = e2_GP[dict_olt_idx_to_new_idx[1]]
        pol_vec_output[i][2][1] = e2_GP[dict_olt_idx_to_new_idx[2]]
        pol_vec_output[i][0][2] = e3_GP[dict_olt_idx_to_new_idx[0]]
        pol_vec_output[i][1][2] = e3_GP[dict_olt_idx_to_new_idx[1]]
        pol_vec_output[i][2][2] = e3_GP[dict_olt_idx_to_new_idx[2]]

    else:
        # P-H
        #omega_trans = np.max([omega[i][0], omega[i][1], omega[i][2]])
        #omega_long = np.min([omega[i][0], omega[i][1], omega[i][2]])
        #omegas_output[i] = np.array([omega_trans,omega_trans,omega_long])

        omegas_output[i] = omega[i]


        kvec_output[i][0] = kvec[i][dict_olt_idx_to_new_idx[0]]
        kvec_output[i][1] = kvec[i][dict_olt_idx_to_new_idx[1]] 
        kvec_output[i][2] = kvec[i][dict_olt_idx_to_new_idx[2]]

        pol_vec_output[i][0][0] = e1_PH[dict_olt_idx_to_new_idx[0]]
        pol_vec_output[i][1][0] = e1_PH[dict_olt_idx_to_new_idx[1]]
        pol_vec_output[i][2][0] = e1_PH[dict_olt_idx_to_new_idx[2]]
        pol_vec_output[i][0][1] = e2_PH[dict_olt_idx_to_new_idx[0]]
        pol_vec_output[i][1][1] = e2_PH[dict_olt_idx_to_new_idx[1]]
        pol_vec_output[i][2][1] = e2_PH[dict_olt_idx_to_new_idx[2]]
        pol_vec_output[i][0][2] = e3_PH[dict_olt_idx_to_new_idx[0]]
        pol_vec_output[i][1][2] = e3_PH[dict_olt_idx_to_new_idx[1]]
        pol_vec_output[i][2][2] = e3_PH[dict_olt_idx_to_new_idx[2]]


out_file = 'scripts/Data/set_pol_vec_input_output/set_pol_vec_100.txt'

# check normalization
for i in range(num_paths * num_kpoints_per_segment):
    a = 2
    #print(f' {pol_vec_output[i][0][0]} {pol_vec_output[i][1][0]} {pol_vec_output[i][2][0]} | {pol_vec_output[i][0][1]} {pol_vec_output[i][1][1]} {pol_vec_output[i][2][1]} | {pol_vec_output[i][0][2]} {pol_vec_output[i][1][2]} {pol_vec_output[i][2][2]}\n')
    norm1 = np.linalg.norm(np.array([pol_vec_output[i][0][0], pol_vec_output[i][1][0], pol_vec_output[i][2][0]]))
    norm2 = np.linalg.norm(np.array([pol_vec_output[i][0][1], pol_vec_output[i][1][1], pol_vec_output[i][2][1]]))
    norm3 = np.linalg.norm(np.array([pol_vec_output[i][0][2], pol_vec_output[i][1][2], pol_vec_output[i][2][2]]))
    #print(f'{norm1} {norm2} {norm3}')


with open(out_file, 'w') as file:
    file.write("qx qy qz freq1 freq2 freq3 v1x v1y v1z v2x v2y v2z v3x v3y v3z\n")
    for i in range(num_paths * num_kpoints_per_segment):
        if i % 100 == 0:
            file.write(f'{kvec_output[i][0]} {kvec_output[i][1]} {kvec_output[i][2]} {omegas_output[i][0]} {omegas_output[i][1]} {omegas_output[i][2]} {pol_vec_output[i][0][0]} {pol_vec_output[i][1][0]} {pol_vec_output[i][2][0]} {pol_vec_output[i][0][1]} {pol_vec_output[i][1][1]} {pol_vec_output[i][2][1]} {pol_vec_output[i][0][2]} {pol_vec_output[i][1][2]} {pol_vec_output[i][2][2]}\n')

