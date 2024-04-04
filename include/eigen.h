#pragma once
#include <iostream>
#include "EigenPCH.h"
#include "util.h"


Eigen::EigenSolver<Eigen::Matrix3d> getEigenValVec(Eigen::Matrix3d dynMat);