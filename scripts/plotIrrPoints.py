import matplotlib.pyplot as plt
import numpy as np


filename = "Outputs/symmApplied.txt"

points = []

with open(filename, 'r') as file:
    file.readline()
    for line in file:
        line = line.strip('\n')

        data = line.split(" ")
        x = float(data[0])
        y = float(data[1])
        z = float(data[2])

        points.append([x, y, z])


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter([point[0] for point in points], [point[1] for point in points], [point[2] for point in points])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.savefig('scripts/Figures/irrPoints/irrPoints_3d.png')
plt.show()


#plt.xlim(-1,1)
#plt.ylim(-1,1)
plt.plot([point[0] for point in points], [point[1] for point in points], 'o')
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('scripts/Figures/irrPoints/irrPoints_xy_20.png')

plt.show()

#plt.xlim(-1,1)
#plt.ylim(-1,1)
plt.plot([point[1] for point in points], [point[2] for point in points], 'o')
plt.xlabel('y')
plt.ylabel('z')
plt.savefig('scripts/Figures/irrPoints/irrPoints_yz_20.png')
plt.show()

#plt.xlim(-1,1)
#plt.ylim(-1,1)
plt.plot([point[0] for point in points], [point[2] for point in points], 'o')
plt.xlabel('x')
plt.ylabel('z')
plt.savefig('scripts/Figures/irrPoints/irrPoints_xz_20.png')

plt.show()



