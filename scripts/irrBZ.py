import numpy as np
from collections import Counter
import matplotlib.pyplot as plt


def read_symmetry_matrices(file_path):
    matrices = []
    with open(file_path) as f:
        for line in f:
            matrix = np.array(line.split(), dtype=int)
            matrix = matrix.reshape((3, 3))
            matrices.append(matrix)
    return matrices


def read_irr_BZ_vectors(file_path):
    vectors = []
    with open(file_path) as f:
        for line in f:
            vector = np.array(line.split(), dtype=float)
            vectors.append(vector)

    return vectors


def read_irr_BZ_vectors_mar(file_path):
    vectors = []
    with open(file_path) as f:
        for line in f:
            vector = np.array(line.split()[1:4], dtype=float)
            vectors.append(vector)

    return vectors


# basis of quantum espresso output
a1 = [2*np.pi, 0, 0]
a2 = [0, 2*np.pi, 0]
a3 = [0, 0, 2*np.pi]

# basis of reciprocal lattice vectors
b1 = 2*np.pi * np.cross(a2, a3) / np.dot(a1, np.cross(a2, a3))
b2 = 2*np.pi * np.cross(a3, a1) / np.dot(a1, np.cross(a2, a3))
b3 = 2*np.pi * np.cross(a1, a2) / np.dot(a1, np.cross(a2, a3))

print(b1, b2, b3)

def findRep(v, irr_BZ_vectors, symm_matrices):
    
    m  = 2

    for n in range(len(symm_matrices)):

        v_trans = np.dot(symm_matrices[n], v)

        for m1 in range(-m, m+1):
            for m2 in range(-m, m+1):
                for m3 in range(-m, m+1):
                
                    for i in range(len(irr_BZ_vectors)):
                
                        v_trans_shifted = np.zeros(3)
                        v_trans_shifted[0] = v_trans[0] + m1 * b1[0] + m2 * b2[0] + m3 * b3[0]
                        v_trans_shifted[1] = v_trans[1] + m1 * b1[1] + m2 * b2[1] + m3 * b3[1] 
                        v_trans_shifted[2] = v_trans[2] + m1 * b1[2] + m2 * b2[2] + m3 * b3[2]

                        if np.linalg.norm(v_trans_shifted - irr_BZ_vectors[i]) < 0.001:
                            return i
    
    return -1



file_path = 'Parameters/symmetryMatricesAll.txt'
matrices = read_symmetry_matrices(file_path)


file_path_irr = "Parameters/8x8x8/irrPoints.txt"

irr_BZ_vectors = read_irr_BZ_vectors(file_path_irr)

irr_BZ_vectors_mar = read_irr_BZ_vectors_mar("Parameters/lifetimes_10_qgrid_10_kgrid_g_ud_g_uu_g_dd.dat")

# Print irreducible BZ vectors
for i, v in enumerate(irr_BZ_vectors_mar):
    print(v)

# Print symmetry matrices
for i, m in enumerate(matrices):
    print(i)
    print(m)



for i, irrv1 in enumerate( irr_BZ_vectors_mar):
    for j, irrv2 in enumerate(irr_BZ_vectors_mar):
        sum = irrv1 - irrv2
        rep = findRep(sum, irr_BZ_vectors_mar, matrices)
        if rep == -1:
            print("No rep found")
            print(irrv1, irrv2, sum)
        else: 
            print("Rep found")
            print(irrv1, irrv2, sum, irr_BZ_vectors_mar[rep])

exit()

applied_symmetries = []

for v in irr_BZ_vectors_mar:
    for n in range(len(matrices)):
        v_trans = np.array([v[0], v[1], v[2]])
        
        for i in range(3):
            v_trans[i] = v[0]*matrices[n][i][0] + v[1]*matrices[n][i][1] + v[2]*matrices[n][i][2]

        #v_trans_norm = v_trans / np.linalg.norm(v_trans) * np.linalg.norm(v)

        for i in range(3):
            if v_trans[i] > 1:
                v_trans[i] -= 1
            elif v_trans[i]< -1:
                v_trans[i] += 1

        applied_symmetries.append(v_trans)
        

        #print(v, v_trans, np.linalg.norm(v), np.linalg.norm(v_trans))


def remove_duplicates(applied_symmetries):
    unique_symmetries = []
    for symmetry in applied_symmetries:
        if not any(np.array_equal(symmetry, unique) for unique in unique_symmetries):
            unique_symmetries.append(symmetry)
   
    return unique_symmetries


unique_applied_symmetries = remove_duplicates(applied_symmetries)
print(len(unique_applied_symmetries))
print(len(applied_symmetries))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter([point[0] for point in applied_symmetries], [point[1] for point in applied_symmetries], [point[2] for point in applied_symmetries])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
plt.clf()

plt.scatter([point[0] for point in applied_symmetries], [point[1] for point in applied_symmetries])
plt.show()

print(len(applied_symmetries))

print(8**3)


