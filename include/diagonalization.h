#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
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

    std::vector<double> eigenEnergies;
    Eigen::MatrixXcd eigenvectors;
    Eigen::MatrixXcd eigenvectors_inverse;

    Diagonalization(const std::vector<CouplingParameter> &couplingParameters, const PhononDispParam &phDispParam, const MagnonDispParam &magDispParam, const Vector3D &kVec);

    void calcCD();
    void calcMatrixHamiltonian();
    void diagonalize();
};

std::vector<std::vector<double>> diagonalizeHamiltonian(const std::vector<Vector3D> &path, const std::vector<PhononDispParam> &phononDispersion, const std::vector<MagnonDispParam> &magnonDispersion, const std::vector<CouplingParameter> &parameters);

Eigen::MatrixXd getMatrix_g();

