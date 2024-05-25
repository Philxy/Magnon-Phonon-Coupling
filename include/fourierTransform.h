#pragma once
#include <complex>
#include <vector>
#include <cassert>
#include "EigenPCH.h"
#include <omp.h>
#include "util.h"
#include "couplingParameters.h"


struct DMILike
{
    double kx, ky, kz;
    std::complex<double> D[3][3]; // first index: nu, second index: mu
    DMILike(double kx, double ky, double kz, const std::vector<CouplingParameter> &parameters);
};



// Represents the DMI-like anisotropic phonon-magnon coupling parameter in the format D_k^{\nu\mu} where \nu = x,y and \mu=x,y,z
struct DMILikeCouplingParam
{
    double kx, ky, kz;
    std::complex<double> D[3][2];

    DMILikeCouplingParam(double kx, double ky, double kz, const std::vector<CouplingParameter> &parameters);
};

std::complex<double> FTD(double kx, double ky, double kz, const std::vector<CouplingParameter> &params, Axis nu, Axis mu);

std::complex<double> FTJiso(double kx, double ky, double kz, const std::vector<CouplingParameter> &params);

std::complex<double> FTJ(double kx, double ky, double kz, const std::vector<CouplingParameter> &params, Axis alpha, Axis beta, Axis mu);


Eigen::Matrix3d dynMat(double kx, double ky, double kz, const std::vector<CouplingParameter> &params);

Eigen::Matrix3d forceMatrix(double x, double y, double z, const std::vector<CouplingParameter> &params);

std::complex<double> J_kq(double kx, double ky, double kz, double qx, double qy, double qz, const std::vector<CouplingParameter> &params, Axis alpha, Axis beta, Axis mu);