import numpy as np
import numpy as np


def read_symmetry_matrices(file_path):
    matrices = []
    with open(file_path) as f:
        for line in f:
            matrix = np.array(line.split(), dtype=int)
            matrix = matrix.reshape((3, 3))
            matrices.append(matrix)
    return matrices


file_path = 'Parameters/symmetrieMatrices.txt'
matrices = read_symmetry_matrices(file_path)


# basis of quantum espresso output
a1 = [2*np.pi, 0, 0]
a2 = [0, 2*np.pi, 0]
a3 = [0, 0, 2*np.pi]

# basis of reciprocal lattice vectors
b1 = 2*np.pi * np.cross(a2, a3) / np.dot(a1, np.cross(a2, a3))
b2 = 2*np.pi * np.cross(a3, a1) / np.dot(a1, np.cross(a2, a3))
b3 = 2*np.pi * np.cross(a1, a2) / np.dot(a1, np.cross(a2, a3))


old_basis = np.array([b1, b2, b3])
new_basis = np.array([a1, a2, a3])

A_inv = np.linalg.inv(new_basis)
transformation_matrix = np.dot(A_inv, old_basis)
transformation_matrix_inv = np.linalg.inv(transformation_matrix)


transformed_matrices = []


for matrix in matrices:
    matrix_in_new_basis = np.dot(
        transformation_matrix_inv, np.dot(matrix, transformation_matrix))
    transformed_matrices.append(matrix_in_new_basis)


with open('Parameters/symmetrieMatrices_transformed.txt', 'w') as f:
    for matrix in transformed_matrices:
        f.write(' '.join([str(x) for x in matrix.flatten()]) + '\n')
