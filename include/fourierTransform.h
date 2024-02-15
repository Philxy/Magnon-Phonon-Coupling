#pragma once
#include <complex>
#include <vector>
#include "couplingParameters.h"
#include <cassert>
#include <Eigen/Dense>

// D_{ijk}^{\nu,\mu}

std::complex<double> FTD(double kx, double ky, double kz, std::vector<CouplingParameter> params, Axis nu, Axis mu);

std::complex<double> FTD(double kx1, double ky1, double kz1, double kx2, double ky2, double kz2, std::vector<CouplingParameter> params, Axis nu, Axis mu);

std::complex<double> FTJiso(double kx, double ky, double kz, std::vector<CouplingParameter> params);

Eigen::Matrix3d dynMat(double kx, double ky, double kz, std::vector<CouplingParameter> params);

Eigen::Matrix3d forceMatrix(double x, double y, double z, std::vector<CouplingParameter> params);
