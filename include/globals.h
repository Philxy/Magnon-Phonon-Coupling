#pragma once
#include "path.h"

// This file contains global variables

extern Vector3D b1;
extern Vector3D b2;
extern Vector3D b3;

extern Eigen::Vector3d b_1;
extern Eigen::Vector3d b_2;
extern Eigen::Vector3d b_3;

constexpr size_t N = 3; // size of the sampling grid (N,N,N)

constexpr double atomicMass = 55.845; // in Dalton

constexpr double S = 1.115;
