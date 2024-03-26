#include "../include/globals.h"

Vector3D b1(2 * 3.14159265358979323846, 2 * 3.14159265358979323846, 0);
Vector3D b2(0, 2 * 3.14159265358979323846, 2 * 3.14159265358979323846);
Vector3D b3(2 * 3.14159265358979323846, 0, 2 * 3.14159265358979323846);

Eigen::Vector3d a_1(1.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0);
Eigen::Vector3d a_2(-1.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0);
Eigen::Vector3d a_3(-1.0 / 2.0, -1.0 / 2.0, 1.0 / 2.0);

Eigen::Vector3d b_1 = 2 * pi * a_2.cross(a_3) / (a_3.dot(a_1.cross(a_2)));
Eigen::Vector3d b_2 = 2 * pi * a_3.cross(a_1) / (a_3.dot(a_1.cross(a_2)));
Eigen::Vector3d b_3 = 2 * pi * a_1.cross(a_2) / (a_3.dot(a_1.cross(a_2)));
