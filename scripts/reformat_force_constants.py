import numpy as np



def gen_nn(n):
        # lattice vectors for bcc
    a1 = [-0.5, 0.5, 0.5]
    a2 = [0.5, -0.5, 0.5]
    a3 = [0.5, 0.5, -0.5]

    # contains points of the lattice relative to the atom located at (0,0,0)
    lattice = []
    distances = []  # contains the distances of the atom from the origin in the order given in 'lattice'

    for i in range(-n, n+1, 1):
        for j in range(-n, n+1, 1):
            for k in range(-n, n+1, 1):
                if i == 0 and j == 0 and k == 0:  # skip origin
                    continue
                x = i*a1[0] + j*a2[0] + k*a3[0]
                y = i*a1[1] + j*a2[1] + k*a3[1]
                z = i*a1[2] + j*a2[2] + k*a3[2]
                dist = np.sqrt(x*x + y*y + z*z)
                distances.append(dist)
                lattice.append([x, y, z])


    # sort the list containing the atoms of the generated lattice based on the relative distance from the atom at the origin
    combined = sorted(zip(distances, lattice))
    a_sorted, b_sorted = zip(*combined)
    a_sorted = list(a_sorted)
    b_sorted = list(b_sorted)


    # dictionary given as: 'distance: [[lattice_point],[lattice_point],...]
    dict = {}

    for i in range(len(distances)):
        d = distances[i]
        dict[d] = []

    for i in range(len(distances)):
        d = distances[i]
        dict[d].append(lattice[i])

    sorted_dict = {k: dict[k] for k in sorted(dict)}

    return sorted_dict





file = open("Parameters/force_constants_bccFe_unformatted_1.txt"  , 'r')


# write to output file:
out_file = open('Parameters/force_constants_bccFe1.txt', 'w')
out_file.write('x,y,x,Phi_xx,Phi_xy,Phi_xz,Phi_yx,Phi_yy,Phi_yz,Phi_zx,Phi_zy,Phi_zz\n')



# Delete the first lines manually to bring the file in a format something like this:
# 6 6 6 
# 1 1 1 1
# 1 1 1 2.24353E-03
# ...
#



grid_line = file.readline() # expexting the first line to contain the grid information, grid is assumed to be quadratic

q_grid_size = int(grid_line.split()[0])


nearest_neighbors = gen_nn(4)


# lattice vectors for bcc
a1 = [-0.5, 0.5, 0.5]
a2 = [0.5, -0.5, 0.5]
a3 = [0.5, 0.5, -0.5]


Phi_tens = np.ndarray(shape=(3,3,q_grid_size**3,2)) # dimension (3,3,total_num_neigh_considered,2) where last dimension is couppling and distance

dist_dict = {}

Phi = np.ndarray(shape=(3,3))


pos = []

for alpha_beta_pairs in range(3*3):
    " There are 3*3 pairs of alpha and beta in toal"
    
    line = file.readline().split()
    alpha = int(line[0])
    beta = int(line[1])

    for site in range(q_grid_size**3):
        line = file.readline().split()
        i, j, k, phi = float(line[0]), float(line[1]), float(line[2]), float(line[3])

        #i = i-1
        #j = j-1
        #k = k-1

        #if i > q_grid_size/2:
        #    i -= q_grid_size
        #if j > q_grid_size/2:
        #    j -= q_grid_size
        #if k > q_grid_size/2:
        #    k -= q_grid_size

        #if i == 0:
        #    i = -q_grid_size/2
        #if j == 0:
        #    j = -q_grid_size/2
        #if k == 0:
        #    k = -q_grid_size/2

        Rx = i * a1[0] + j * a2[0] + k * a3[0]
        Ry = i * a1[1] + j * a2[1] + k * a3[1]
        Rz = i * a1[2] + j * a2[2] + k * a3[2]

        pos.append([Rx,Ry,Rz])

        distance = np.sqrt(Rx**2+Ry**2+Rz**2)

        Phi_tens[alpha-1][beta-1][site][0] = phi
        Phi_tens[alpha-1][beta-1][site][1] = distance

        #if distance not in dist_dict:
        #    dist_dict[distance] = np.ndarray(shape=(3,3))
        #    dist_dict[distance][alpha-1][beta-1] = phi
        #else:
        #    if dist_dict[distance][alpha-1][beta-1] != phi:
        #        print(dist_dict[distance][alpha-1][beta-1], phi)
        #        
        #    dist_dict[distance][alpha-1][beta-1] = phi

        #if distance == 1.4142135623730951 and alpha == 1 and beta == 1:
        ##    print(Phi_tens[alpha-1][beta-1][site][0])
        #    print(i,j,k)

        #    Phi_tens[alpha-1][beta-1][site][0] = 0
        #    Phi_tens[alpha-1][beta-1][site][1] = 0

#for dist in dist_dict:
#    phi = dist_dict[dist]
#
#    for position in nearest_neighbors[dist]:
#        
#        x,y,z = position
#        Phi_xx = phi[0][0]
#        Phi_xy = phi[0][1]
#        Phi_xz = phi[0][2]
#        Phi_yx = phi[1][0]
#        Phi_yy = phi[1][1]
#        Phi_yz = phi[1][2]
#        Phi_zx = phi[2][0]
#        Phi_zy = phi[2][1]
#        Phi_zz = phi[2][2]
#        line = ''
#        line += str(x) + "," + str(y) + str(",") + str(y) + ","
#        line += str(Phi_xx) + "," + str(Phi_xy) + "," + str(Phi_xz) + ","
#        line += str(Phi_yx) + "," + str(Phi_yy) + "," + str(Phi_yz) + ","
#        line += str(Phi_zx) + "," + str(Phi_zy) + "," + str(Phi_zz) + "\n"
#        out_file.write(line)
#
#exit()


for site in range(q_grid_size**3):
    x,y,z = pos[site]
    Phi_xx = Phi_tens[0][0][site][0]
    Phi_xy = Phi_tens[0][1][site][0]
    Phi_xz = Phi_tens[0][2][site][0]
    Phi_yx = Phi_tens[1][0][site][0]
    Phi_yy = Phi_tens[1][1][site][0]
    Phi_yz = Phi_tens[1][2][site][0]
    Phi_zx = Phi_tens[2][0][site][0]
    Phi_zy = Phi_tens[2][1][site][0]
    Phi_zz = Phi_tens[2][2][site][0]
    line = ''
    line += str(x) + "," + str(y) + str(",") + str(y) + ","
    line += str(Phi_xx) + "," + str(Phi_xy) + "," + str(Phi_xz) + ","
    line += str(Phi_yx) + "," + str(Phi_yy) + "," + str(Phi_yz) + ","
    line += str(Phi_zx) + "," + str(Phi_zy) + "," + str(Phi_zz) + "\n"
    out_file.write(line)


exit()

for site in range(q_grid_size**3):
    x,y,z = pos[site]
    Phi_xx = Phi_tens[0][0][site][0]
    Phi_xy = Phi_tens[0][1][site][0]
    Phi_xz = Phi_tens[0][2][site][0]
    Phi_yx = Phi_tens[1][0][site][0]
    Phi_yy = Phi_tens[1][1][site][0]
    Phi_yz = Phi_tens[1][2][site][0]
    Phi_zx = Phi_tens[2][0][site][0]
    Phi_zy = Phi_tens[2][1][site][0]
    Phi_zz = Phi_tens[2][2][site][0]
    line = ''
    line += str(x) + "," + str(y) + str(",") + str(y) + ","
    line += str(Phi_xx) + "," + str(Phi_xy) + "," + str(Phi_xz) + ","
    line += str(Phi_yx) + "," + str(Phi_yy) + "," + str(Phi_yz) + ","
    line += str(Phi_zx) + "," + str(Phi_zy) + "," + str(Phi_zz) + "\n"
    out_file.write(line)