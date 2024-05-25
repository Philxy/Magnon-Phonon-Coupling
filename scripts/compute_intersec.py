import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize

kVectors = []
energies = []
mag_disp = []
ph_disp = []



with open("Outputs/wholeBZ/20x20x20BZ_mag.txt", "r") as file:
    file.readline()
    for line in file:
        data = line.split(',')

        mag_disp.append(float(data[3].replace('\n','').replace('(','')))
        


with open("Outputs/wholeBZ/20x20x20_disp.txt", "r") as file:
    file.readline()
    for line in file:
        data = line.split(' ')

        ph_disp.append([float(d) for d in data[3:6]])



with open("Outputs/wholeBZ/20x20x20_disp.txt", "r") as file:
    file.readline()
    for line in file:
        data = line.split(' ')

        kVectors.append([float(k) for k in data[:3]])



def find_intersection_kVec(kVectors, mag_disp, ph_disp):

    threshold = 1
    kVecs_close_to_interception = []

    for i in range(len(mag_disp)):
        for branch in range(3):
            if abs(mag_disp[i] - ph_disp[i][branch]) < threshold:
                kVecs_close_to_interception.append(kVectors[i])

    return kVecs_close_to_interception


def sample_close_to_some_kVec(kVec, range, num_setps):

    kVectors = []

    for i in np.linspace(-range, range, num_setps):
        for j in np.linspace(-range, range, num_setps):
            for k in np.linspace(-range, range, num_setps):
                kVectors.append([kVec[0]+i, kVec[1]+j, kVec[2]+k])

    return kVectors


def get_kVecs_close_to_intercep(kVecs, threshold, range, num_steps, mag_disp, ph_disp):
    
        kVecs_intersec = find_intersection_kVec(kVecs, mag_disp, ph_disp)
    
        res = []

        for kVec in kVecs_intersec:
            res.extend(sample_close_to_some_kVec(kVec, range, num_steps))
    
        return res


def print_kVec_for_QE_to_file(kVectors, filepath):

    with open(filepath, "w") as file:
        for k in kVectors:
            file.write(f'{k[0]/(2*np.pi)} {k[1]/(2*np.pi)} {k[2]/(2*np.pi)} 1\n')

        

# range near intersecting points: hybridization range is roughly 0.04 * 2pi/a.
# 
kVecs_close_to_all_intercep = get_kVecs_close_to_intercep(kVectors, threshold=1, range=0.1*2*np.pi, num_steps=10, mag_disp=mag_disp, ph_disp=ph_disp)

print('num kVecs close to all intersections: ', len(kVecs_close_to_all_intercep))


print_kVec_for_QE_to_file(kVecs_close_to_all_intercep, "Outputs/wholeBZ/kVecs_close_to_all_intercep_acc.txt")

interception_kVec = find_intersection_kVec(kVectors, mag_disp, ph_disp)

print('num intersecting k Points: ' ,len(interception_kVec))

