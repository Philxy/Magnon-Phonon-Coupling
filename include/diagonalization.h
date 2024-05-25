#pragma once
#include <iostream>
#include "EigenPCH.h"
#include "path.h"
#include "dispersion.h"
#include "util.h"
#include "globals.h"

struct Diagonalization
{
    PhononDispParam phDisp;
    MagnonDispParam magDisp;
    Vector3D k;
    std::vector<std::complex<double>> C, D;
    std::vector<CouplingParameter> couplingParameters;
    Eigen::MatrixXcd matrixHamiltonian;
    std::vector<std::complex<double>> angularMomentum;
    std::vector<double> eigenEnergies;
    Eigen::MatrixXcd eigenvectors;
    Eigen::MatrixXcd eigenvectors_inverse;
    std::vector<std::complex<double>> dmiLike; // Dxx, Dxy, Dxz, Dyx, Dyy, Dyz

    std::complex<double> A, B;

    Diagonalization(const std::vector<CouplingParameter> &couplingParameters, const PhononDispParam &phDispParam, const MagnonDispParam &magDispParam, const Vector3D &kVec);

    void calcCD();
    void calcCD_new();
    void calcMatrixHamiltonian();
    void diagonalize();
    void calcDMILike();
    void calcAB();
    void calcAngularMomentum();


};

std::vector<std::vector<double>> diagonalizeHamiltonian(const std::vector<Vector3D> &path, const std::vector<PhononDispParam> &phononDispersion, const std::vector<MagnonDispParam> &magnonDispersion, const std::vector<CouplingParameter> &parameters);

Eigen::MatrixXd getMatrix_g();

