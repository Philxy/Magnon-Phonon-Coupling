#pragma once
#include "fourierTransform.h"
#include "dispersion.h"
#include "EigenPCH.h"
#include "globals.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <functional>

struct BZVec
{
    Eigen::Vector3d kVec;    // k vector
    double n_ph[3], n_mag;   // phonon and magnon occupation number
    double multiplicity;     // multiplicity of the k vector
    PhononDispParam phDisp;  // energy and pol vectors of the phonon
    MagnonDispParam magDisp; // energy of the magnon mode

    BZVec(Eigen::Vector3d kVec, PhononDispParam phDisp, MagnonDispParam magDisp, double multiplicity)
    {
        this->kVec = kVec;
        this->phDisp = phDisp;
        this->magDisp = magDisp;
        this->multiplicity = multiplicity;
    };
};

struct SymmetrieGroup
{
    std::vector<Eigen::Matrix3d> symmetrieOperations;

    void init(std::string symmetryFile);
    std::vector<Eigen::Vector3d> applySymmetries(const Eigen::Vector3d &kVec);
};

struct IrreducibleBZ
{

    std::vector<Eigen::Vector3d> irreducibleBZVectors; // irreducible BZ vectors
    std::vector<double> multiplicities;                // multiplicities of the irreducible BZ vectors

    std::vector<MagnonDispParam> magDisp; // magnon dispersion
    std::vector<PhononDispParam> phDisp;  // phonon dispersion

    std::vector<double> magOccNumbers;               // magnon occupation numbers
    std::vector<std::array<double, 3>> phOccNumbers; // phonon occupation numbers

    // COEFFICIENTS

    // first order in mag variables
    std::vector<std::array<std::complex<double>, 3>> CGrid;
    std::vector<std::array<std::complex<double>, 3>> DGrid;
    std::vector<std::array<std::complex<double>, 3>> CGrid_negative_sign;
    std::vector<std::array<std::complex<double>, 3>> DGrid_negative_sign;
    // second order in mag variables
    std::vector<std::vector<std::array<std::complex<double>, 3>>> gammaPlusGrid;  // \Gamma^+(k,q)
    std::vector<std::vector<std::array<std::complex<double>, 3>>> gammaMinusGrid; // \Gamma^-(k,q)
    std::vector<std::vector<std::array<std::complex<double>, 3>>> gammaZGrid;     // \Gamma^z(k,q)

    std::vector<std::vector<std::array<std::complex<double>, 3>>> gammaPlusGrid_negativeSign;  // \Gamma^+(k,-q)
    std::vector<std::vector<std::array<std::complex<double>, 3>>> gammaMinusGrid_negativeSign; // \Gamma^-(k,-q)
    std::vector<std::vector<std::array<std::complex<double>, 3>>> gammaZGrid_negativeSign;     // \Gamma^z(k,-q)

    // precomputed representatives satisfying pseudo conservation laws: G = +- k +- k' +- q
    std::vector<std::vector<int>> k_prime_representatives_gammaZ_minus_q, k_prime_representatives_gammaM_minus_q, k_prime_representatives_gammaP_minus_q, k_prime_representatives_gammaZ_plus_q, k_prime_representatives_gammaM_plus_q, k_prime_representatives_gammaP_plus_q;

    SymmetrieGroup symmetryGroup; // symmetry group

    void init(std::string irreducibleBZFile); // initializes the irreducible BZ vectors
    void initMagnonDisp(std::string couplingParameterFile);                              // initializes the magnon dispersion
    void initPhononDisp(std::string dynamicMatricesFile, std::string nextNeighbourFile); // initializes the phonon dispersion by doing the diagonalization

    // initializes the magnon and phonon dispersion from files
    void initMagnonPhononDispFromFile(std::string filepathPh, std::string filepathMag, std::string magnonDispOutputPath);

    void initOccNumbers(double thermalEnergyPh, double thermalEnergyMag);                                                            // initializes the occupation numbers
    void initCoefficients(const std::vector<CouplingParameter> &parameters, int ftN); // initializes the coefficients
    void readMultiplicities(const std::string &filename);                             // retrieves the multiplicities from a file
    void init_k_prime();

    int findRepresentative(const Eigen::Vector3d &kVec); // returns the index of the representative in the irreducibleBZVectors vector
    Eigen::Vector3d findClosestVector(const Eigen::Vector3d &kVec);
    int findClosestVectorIndex(const Eigen::Vector3d &kVec);

    Eigen::Vector3d getG_Z(const Eigen::Vector3d &k, const Eigen::Vector3d &q);
    Eigen::Vector3d getG_gammaM(const Eigen::Vector3d &kVec, const Eigen::Vector3d &qVec);
    Eigen::Vector3d getG_gammaP(const Eigen::Vector3d &kVec, const Eigen::Vector3d &qVec);

    void saveCoefficientsAsSqrtAbs(std::string filename);
    void readCoefficients(std::string filename);

    void integrate();


    std::vector<std::function<Eigen::Vector3d(Eigen::Vector3d)>> operations;

    void initSymmOp(std::string filepath);
    
    std::vector<Eigen::Vector3d> applySymmOp(const Eigen::Vector3d &vec);



};

double distance(const Eigen::Vector3d &k1, const Eigen::Vector3d &k2);
