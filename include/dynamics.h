#pragma once
#include "sampling.h"
#include "fourierTransform.h"

struct Coefficients
{

    std::complex<double> getGammaPlus(BZVector k1, BZVector k2);  // term second order in magnon creation operator
    std::complex<double> getGammaMinus(BZVector k1, BZVector k2); // term second order in magnon annihilation operator
    std::complex<double> getGammaZ(BZVector k1, BZVector k2);     // term with creation and annihilation of magnon

    std::vector<Vector3D> reciprocalLatticeVectors;

    double S = 1.1;
    double ftN = 1; // welchen wert hat das???
    double atomic_mass = 55.56;

    SamplingGrid samplingGrid;

    void initReciprocalLatticeVectors();                         // must be initialized!
    void init(const std::vector<CouplingParameter> &parameters); // must be initialized!

    std::complex<double> gammaPlus[N][N][N][N][N][N][3];
    std::complex<double> gammaMinus[N][N][N][N][N][N][3];
    std::vector<std::complex<double>[N][N][N][N][N][N][3]> gammaZValues;
};


