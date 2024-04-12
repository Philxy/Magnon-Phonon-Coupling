#pragma once
#include "path.h"
#include "Eigen/Dense"

// This file contains global variables

extern Vector3D b1;
extern Vector3D b2;
extern Vector3D b3;

extern Eigen::Vector3d a_1;
extern Eigen::Vector3d a_2;
extern Eigen::Vector3d a_3;

extern Eigen::Vector3d b_1;
extern Eigen::Vector3d b_2;
extern Eigen::Vector3d b_3;

constexpr double pi = 3.14159265358979323846;

constexpr size_t N = 3; // size of the sampling grid (N,N,N)

constexpr double atomicMass = 55.845; // in Dalton

constexpr double S = 2 * 2.25 / 4;

constexpr int nFT = 135160; // number of Fourier coefficients per displacement axis