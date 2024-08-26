#pragma once
#include <iostream>
#include "EigenPCH.h"
#include "path.h"
#include "dispersion.h"
#include "util.h"
#include "globals.h"
#include <iomanip>
#include "Eigen/Cholesky"
#include <Eigen/LU> 

// Define the components of the SLC tensor to be included
enum Interaction
{
    DMI,
    ANISOTROPY,
    ALL,
};

enum ComputationInfo
{
    SUCCESS,
    FAILURE
};

struct Diagonalization
{
    PhononDispParam phDisp;
    MagnonDispParam magDisp;
    Vector3D k;
    std::vector<std::complex<double>> C, D; // coefficients of the hybridization term in the Hamiltonian
    std::vector<CouplingParameter> couplingParameters;
    Eigen::MatrixXcd matrixHamiltonian;                // Hamiltonian in matrix representation
    std::vector<std::complex<double>> angularMomentum; // Lx, Ly, Lz parallel to the k vector
    std::vector<Eigen::Vector3d> angularMomentumFromEigenvectors;
    std::vector<Eigen::Vector3cd> polVectors;

    std::vector<double> eigenEnergies; //
    Eigen::MatrixXcd eigenvectors;
    Eigen::MatrixXcd eigenvectors_inverse;
    std::vector<std::complex<double>> dmiLike; // Dxx, Dxy, Dxz, Dyx, Dyy, Dyz

    std::complex<double> A, B;

    Diagonalization(const std::vector<CouplingParameter> &couplingParameters, const PhononDispParam &phDispParam, const MagnonDispParam &magDispParam, const Vector3D &kVec);

    void roundMatrixHamiltonian();

    void calcCD();
    void calculateCD(Interaction interaction);
    void calcMatrixHamiltonian();
    void diagonalize();
    ComputationInfo diagonalizeColpa();
    void calcDMILike();
    void calcAB();
    void calcAngularMomentum();
    void calcAngularMomentumFromEigenvectors();
    void calcPolVectors();
};

std::vector<std::vector<double>> diagonalizeHamiltonian(const std::vector<Vector3D> &path, const std::vector<PhononDispParam> &phononDispersion, const std::vector<MagnonDispParam> &magnonDispersion, const std::vector<CouplingParameter> &parameters);

Eigen::MatrixXd getMatrix_g();
