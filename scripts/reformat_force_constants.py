import numpy as np


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


# lattice vectors for bcc
a1 = [-0.5, 0.5, 0.5]
a2 = [0.5, -0.5, 0.5]
a3 = [0.5, 0.5, -0.5]

Phi_tens = np.ndarray(shape=(3,3,q_grid_size**3)) # dimension (3,3,total_num_neigh_considered)
Coordinates = []

for alpha_beta_pairs in range(3*3):
    " There are 3*3 pairs of alpha and beta in toal"
    
    line = file.readline().split()
    alpha = int(line[0])
    beta = int(line[1])


    for site in range(q_grid_size**3):
        line = file.readline().split()
        i, j, k, phi = float(line[0]), float(line[1]), float(line[2]), float(line[3])

        i = i-1
        j = j-1
        k = k-1
        Rx = i * a1[0] + j * a2[0] + k * a3[0]
        Ry = i * a1[1] + j * a2[1] + k * a3[1]
        Rz = i * a1[2] + j * a2[2] + k * a3[2]

        Coordinates.append([Rx,Ry,Rz])

        Phi_tens[alpha-1][beta-1][site] = phi

        if alpha_beta_pairs > 0:
            continue
        print([Rx,Ry,Rz], np.sqrt(Rx**2+Ry**2+Rz**2))




for site in range(q_grid_size**3):
    x = Coordinates[site][0]
    y = Coordinates[site][1]
    z = Coordinates[site][2]
    Phi_xx = Phi_tens[0][0][site]
    Phi_xy = Phi_tens[0][1][site]
    Phi_xz = Phi_tens[0][2][site]
    Phi_yx = Phi_tens[1][0][site]
    Phi_yy = Phi_tens[1][1][site]
    Phi_yz = Phi_tens[1][2][site]
    Phi_zx = Phi_tens[2][0][site]
    Phi_zy = Phi_tens[2][1][site]
    Phi_zz = Phi_tens[2][2][site]
    line = ''
    line += str(x) + "," + str(y) + str(",") + str(y) + ","
    line += str(Phi_xx) + "," + str(Phi_xy) + "," + str(Phi_xz) + ","
    line += str(Phi_yx) + "," + str(Phi_yy) + "," + str(Phi_yz) + ","
    line += str(Phi_zx) + "," + str(Phi_zy) + "," + str(Phi_zz) + "\n"
    out_file.write(line)