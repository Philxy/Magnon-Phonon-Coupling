#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "util.h"


Eigen::EigenSolver<Eigen::Matrix3d> getEigenValVec(Eigen::Matrix3d dynMat);