import numpy as np
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


def apply_symmetry_operation(operation, vector):
    """
    Apply a symmetry operation to a vector.
    
    Parameters:
        operation (str): A string representation of the operation (e.g., "-x,y,-z").
        vector (np.array): A numpy array representing the vector to which the operation is applied.
        
    Returns:
        np.array: The transformed vector.
    """
    # Map the operation symbols to the vector components.
    op_map = {'x': vector[0], 'y': vector[1], 'z': vector[2], '-x': -vector[0], '-y': -vector[1], '-z': -vector[2]}
    
    # Split the operation into its components.
    ops = operation.split(',')
    
    # Apply the operations.
    transformed_vector = np.array([op_map[ops[0]], op_map[ops[1]], op_map[ops[2]]])
    
    return transformed_vector

def read_operations_file(filename):
    """
    Read the operations from a file.
    
    Parameters:
        filename (str): The path to the file containing the operations.
        
    Returns:
        list: A list of operations (as strings).
    """
    operations = []
    with open(filename, 'r') as file:
        for line in file:
            operation = line.strip('\n')
            print(operation)
            operations.append(operation.strip('\n'))
    return operations


def findRep(v, irr_BZ_vectors, symm_operations):
    
    m  = 2

    for idx, operation in enumerate(symm_operations, start=1):

        transformed_vector = apply_symmetry_operation(operation, v)
        #print(operation, v, transformed_vector)
        #print(np.linalg.norm(v)-np.linalg.norm(transformed_vector))

        for m1 in range(-m, m+1):
            for m2 in range(-m, m+1):
                for m3 in range(-m, m+1):
                
                    for i in range(len(irr_BZ_vectors)):
                
                        v_trans_shifted = np.zeros(3)
                        v_trans_shifted[0] = transformed_vector[0] + m1 * b1[0] + m2 * b2[0] + m3 * b3[0]
                        v_trans_shifted[1] = transformed_vector[1] + m1 * b1[1] + m2 * b2[1] + m3 * b3[1] 
                        v_trans_shifted[2] = transformed_vector[2] + m1 * b1[2] + m2 * b2[2] + m3 * b3[2]

                        if np.linalg.norm(v_trans_shifted - irr_BZ_vectors[i]) < 0.001:
                            return i
    
    return -1


def main(filename='scripts/Data/tu_graz_symm_bcc.txt'):
    """
    Main function to apply all symmetry operations to a given vector.
    
    Parameters:
        vector (tuple): The vector to which the operations will be applied.
        filename (str): The path to the file containing the operations.
    """
    operations = read_operations_file(filename)
    
    #for idx, operation in enumerate(operations, start=1):
    #    transformed_vector = apply_symmetry_operation(operation, vector_np)
    #    print(f"Operation {idx}: {operation} -> Transformed Vector: {transformed_vector}")



    # Irr BZ stuff


    file_path_irr = "Parameters/8x8x8/irrPoints.txt"
    irr_BZ_vectors = read_irr_BZ_vectors(file_path_irr)
    irr_BZ_vectors_mar = read_irr_BZ_vectors_mar("Parameters/lifetimes_10_qgrid_10_kgrid_g_ud_g_uu_g_dd.dat")

    # Print irreducible BZ vectors
    for i, v in enumerate(irr_BZ_vectors_mar):
        print(v)


    x = []

    for i, irrv1 in enumerate( irr_BZ_vectors):
        for op in operations:
            transformed_vector = apply_symmetry_operation(op, irrv1)
            x.append(transformed_vector)
    

    # Convert the x list to a numpy array
    x = np.array(x)

    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the x, y, and z coordinates
    ax.scatter(x[:, 0], x[:, 1], x[:, 2])

    # Set labels for the x, y, and z axes
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Show the plot
    plt.show()

    fail_count = 0
    success_count = 0

    for i, irrv1 in enumerate( irr_BZ_vectors):
        for j, irrv2 in enumerate(irr_BZ_vectors):
            sum = irrv1 + irrv2
            rep = findRep(sum, irr_BZ_vectors, operations)
            if rep == -1:
                print("No rep found")
                print(irrv1, irrv2, sum)
                fail_count += 1
            else: 
                #print("Rep found")
                #print(irrv1, irrv2, sum, irr_BZ_vectors_mar[rep])
                success_count += 1

            print(f"Success: {success_count}, Fail: {fail_count}")


main()




