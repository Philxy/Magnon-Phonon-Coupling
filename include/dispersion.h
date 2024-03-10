#pragma once
#include <iostream>
#include <complex>
#include <fftw3.h>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "couplingParameters.h"
#include "util.h"
#include "path.h"
#include "fourierTransform.h"
#include "util.h"

struct PhononDispParam
{
    double kx, ky, kz;
    double E[3];             // energy for the phonon modes
    double polVectors[3][3]; // polarization vectors - each column represents a pol vec
};

struct MagnonDispParam
{
    double kx, ky, kz;
    double energy;
};

std::vector<MagnonDispParam> getMagneticDispersion(std::string couplingParameterFile, std::string magnonDispOutputPath, std::vector<Vector3D> path, double S);
std::vector<MagnonDispParam> getMagneticDispersion(std::string couplingParameterFile, std::vector<Vector3D> path, double S);

std::vector<PhononDispParam> getPhononDispersion(std::string dynamicMatricesFile, std::string nextNeighbourFile, std::string phononDispOutputPath, std::vector<Vector3D> path, double atomic_mass);
std::vector<PhononDispParam> getPhononDispersion(std::string dynamicMatricesFile, std::string nextNeighbourFile, std::vector<Vector3D> path, double atomic_mass);
