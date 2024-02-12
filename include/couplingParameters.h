#pragma once
#include <iostream>
#include <complex>
#include <fftw3.h>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "couplingParameters.h"

// Kartesian coordinate axes
enum Axis
{
    X,
    Y,
    Z
};

struct CouplingParameter
{
    double x, y, z;
    double J_xx, J_xy, J_xz, J_yx, J_yy, J_yz, J_zx, J_zy, J_zz;
    double x_ij, y_ij, z_ij, x_jk, y_jk, z_jk, x_ki, y_ki, z_ki;
    double J[3][3];
    double J_iso;
    double force_constant;
    Axis displacementDirection;

    std::string printParameter();
};

std::vector<CouplingParameter> readCouplingParameters(const std::string &filename, Axis displacementDirection);
std::vector<CouplingParameter> readCouplingParametersIso(const std::string &filename);
