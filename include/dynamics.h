#pragma once
#include "fourierTransform.h"
#include "dispersion.h"
#include <Eigen/Dense>
#include "globals.h"
#include <fstream>
#include <sstream>

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
    std::vector<Eigen::Vector3d> applySymmetries(Eigen::Vector3d kVec);
};

struct IrreducibleBZ
{

    std::vector<Eigen::Vector3d> irreducibleBZVectors;

    std::vector<MagnonDispParam> magDisp;
    std::vector<PhononDispParam> phDisp;

    SymmetrieGroup symmetrieGroup;

    void init(std::string irreducibleBZFile);
    void initMagnonDisp(std::string couplingParameterFile);
    void initPhononDisp(std::string dynamicMatricesFile, std::string nextNeighbourFile);

    int findRepresentative(Eigen::Vector3d kVec); // returns the index of the representative in the irreducibleBZVectors vector
};

bool insideFirstBZ(Eigen::Vector3d kVec);
Eigen::Vector3d mapToFirstBZ(Eigen::Vector3d kVec);
double distance(Eigen::Vector3d k1, Eigen::Vector3d k2);