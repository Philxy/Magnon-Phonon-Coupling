import re
import numpy as np
import os


def extract_q_vectors_and_matrices(file_path):
    with open(file_path, 'r') as file:
        content = file.read()

    # Pattern to match q vector and its matrix
    pattern = r"q = \(\s*([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)\s*\)\s*\n\s*1\s*1\s*\n([\s\S]+?)\n\n"
    matches = re.findall(pattern, content)

    results = []

    for match in matches:
        q_vector = match[:3]  # First three groups are the q vector components
        q_vector = [float(qi) for qi in q_vector]
        matrix_str = match[3]  # Fourth group is the matrix

        # Extract non-zero matrix elements
        matrix_lines = matrix_str.strip().split('\n')
        matrix = []
        for line in matrix_lines:
            row = [float(val)
                   for i, val in enumerate(line.split()) if i % 2 == 0]
            if row:  # Only add non-empty rows
                matrix.append(row)

        # Store the q vector and its corresponding matrix
        results.append((q_vector, matrix))

    return results


def print_vec_and_mat_as_one_line(vec, mat):
    result = ''
    for i in range(3):
        result += str(vec[i]) + ","

    for i in range(3):
        for j in range(3):
            result += str(mat[i][j]) + ","

    result = result[:-1] + '\n'

    return result


a = 1
b1 = [2 * np.pi / a, 2 * np.pi / a, 0]
b2 = [0, 2 * np.pi / a, 2 * np.pi / a]
b3 = [2 * np.pi / a, 0, 2 * np.pi / a]

a1 = [1.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0]
a2 = [-1.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0]
a3 = [-1.0 / 2.0, -1.0 / 2.0, 1.0 / 2.0]


b_1 = 2*np.pi * np.cross(a2, a3) / np.dot(a1, np.cross(a2, a3))
b_2 = 2*np.pi * np.cross(a3, a1) / np.dot(a1, np.cross(a2, a3))
b_3 = 2*np.pi * np.cross(a1, a2) / np.dot(a1, np.cross(a2, a3))

directory_path = 'scripts/Data/dyn20x20x20'

out_file = open('Parameters/dynMat_20x20x20.txt', 'w')

out_file.write('kx,ky,kz,Dxx,Dxy,Dxz,Dyx,Dyy,Dyz,Dzx,Dzy,Dzz\n')

multiplicity_file = open('Parameters/multiplicity_20x20x20.txt', 'w')


for filename in os.listdir(directory_path):
    file_path = os.path.join(directory_path, filename)
    if os.path.isfile(file_path):
        print(f"Processing {filename}...")
        q_vectors_and_matrices = extract_q_vectors_and_matrices(file_path)

        multiplicity = len(q_vectors_and_matrices)
        multiplicity_file.write(f"{filename},{multiplicity}\n")

        for q_vector, matrix in q_vectors_and_matrices:

            #u, v, w = q_vector
            #q_vector[0] = u * b_1[0] + v * b_2[0] + w * b_3[0]
            #q_vector[1] = u * b_1[1] + v * b_2[1] + w * b_3[1]
            #q_vector[2] = u * b_1[2] + v * b_2[2] + w * b_3[2]

            q_vector[0] *= 2*np.pi
            q_vector[1] *= 2*np.pi
            q_vector[2] *= 2*np.pi

            out_file.write(print_vec_and_mat_as_one_line(q_vector, matrix))


exit()

# Replace 'file_path' with the actual path to your file
file_path = 'scripts/Data/bccfe.dyn1'
q_vectors_and_matrices = extract_q_vectors_and_matrices(file_path)


for q_vector, matrix in q_vectors_and_matrices:
    print(f"q vector: {q_vector}")
    print("Matrix:")
    for row in matrix:
        print(row)
    print()
    print(print_vec_and_mat_as_one_line(q_vector, matrix))
