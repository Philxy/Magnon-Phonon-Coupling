#include "../include/dynamics.h"

double BoseEinstein(double energy, double thermalEnergy)
{
    return 1.0 / (exp(energy / thermalEnergy) - 1);
}

// calculates the thermal energy from the Bose-Einstein distribution given the energy and the occupation number
double BoseEinsteinInverted(double energy, double occupationNumber)
{
    return energy / log(1.0 / occupationNumber + 1);
}

double distance(const Eigen::Vector3d &k1, const Eigen::Vector3d &k2)
{
    Eigen::Vector3d diff = k1 - k2;

    return sqrt(diff(0) * diff(0) + diff(1) * diff(1) + diff(2) * diff(2));
}

std::vector<Eigen::Vector3d> SymmetrieGroup::applySymmetries(const Eigen::Vector3d &kVec)
{

    std::vector<Eigen::Vector3d> kVecs;
    for (Eigen::Matrix3d symOp : symmetrieOperations)
    {
        kVecs.push_back(symOp * kVec);
    }
    return kVecs;
}

int IrreducibleBZ::findRepresentative(const Eigen::Vector3d &kVec)
{
    const double dist_threshold = 0.001;

    const int indices[9] = {0, 1, -1, 2, -2, 3, -3, 4, -4};

    for (auto transformedVec : applySymmOp(kVec))
    {
        for (int i : indices)
        {
            for (int j : indices)
            {
                for (int k : indices)
                {
                    Eigen::Vector3d shifted_transformed_kVec = transformedVec + i * Eigen::Vector3d(2 * pi, 0, 0) + j * Eigen::Vector3d(0, 2 * pi, 0) + k * Eigen::Vector3d(0, 0, 2 * pi);

                    for (int n = 0; n < irreducibleBZVectors.size(); n++)
                    {
                        double dist = distance(irreducibleBZVectors.at(n), shifted_transformed_kVec);
                        if (dist < dist_threshold)
                        {
                            return n;
                        }
                    }
                }
            }
        }
    }
    std::cout << "Could not find rep in the IBZ for: " << kVec(0) / (2 * pi) << "," << kVec(1) / (2 * pi) << "," << kVec(2) / (2 * pi) << std::endl;

    return -1;

    /*

    // try mapping to first BZ
    for (int i = -m; i <= m; i++)
    {
        for (int j = -m; j <= m; j++)
        {
            for (int k = -m; k <= m; k++)
            {
                // shift the transformed k vector by the reciprocal lattice vectors
                const Eigen::Vector3d shifted_kVec = kVec + i * Eigen::Vector3d(2 * pi, 0, 0) + j * Eigen::Vector3d(0, 2 * pi, 0) + k * Eigen::Vector3d(0, 0, 2 * pi); // i * b_1 + j * b_2 + k * b_3; //

                for (Eigen::Vector3d transformed_shifted_kVec : symmetryGroup.applySymmetries(shifted_kVec))
                {


                    for (int n = 0; n < irreducibleBZVectors.size(); n++)
                    {
                        double dist = distance(irreducibleBZVectors.at(n), transformed_shifted_kVec);
                        // double dist = distance(irreducibleBZVectors.at(n), transformed_shifted_kVec);
                        if (dist < dist_threshold)
                        {
                            // std::cout << "dist: " << dist << std::endl;
                            return n;
                        }
                        // if (dist < curr_min_dist)
                        //{
                        //     curr_min_dist = dist;
                        //     curr_min_idx = n;
                        //     curr_transformed_kVec = transformed_shifted_kVec;
                        // }
                    }
                }
            }
        }
    }


    */

    // if (curr_min_dist > dist_threshold)
    //{
    //     std::cout << "Failure: Did not find representative" << std::endl;
    //     std::cout << "dist = " << curr_min_dist << std::endl;
    //     std::cout << "vec = " << kVec(0) / (2 * pi) << "," << kVec(1) / (2 * pi) << "," << kVec(2) / (2 * pi) << std::endl;
    //     std::cout << "rep = " << irreducibleBZVectors.at(curr_min_idx)(0) / (2 * pi) << "," << irreducibleBZVectors.at(curr_min_idx)(1) / (2 * pi) << "," << irreducibleBZVectors.at(curr_min_idx)(2) / (2 * pi) << std::endl;
    //
    //     return -1;
    // }
    // return -1;

    // std::cout << "Success: Found representative" << std::endl;
    // std::cout << "dist = " << curr_min_dist << std::endl;
    // std::cout << "v = " << kVec(0) / (2 * pi) << "," << kVec(1) / (2 * pi) << "," << kVec(2) / (2 * pi) << std::endl;
    // std::cout << "r = " << irreducibleBZVectors.at(curr_min_idx)(0) / (2 * pi) << "," << irreducibleBZVectors.at(curr_min_idx)(1) / (2 * pi) << "," << irreducibleBZVectors.at(curr_min_idx)(2) / (2 * pi) << std::endl;

    std::cout << "Could not find rep in the IBZ for: " << kVec(0) / (2 * pi) << "," << kVec(1) / (2 * pi) << "," << kVec(2) / (2 * pi) << std::endl;
    return -1;
    // return curr_min_idx;

    // std::cout << "k = " << kVec(0) << "," << kVec(1) << "," << kVec(2) << ","
    //           << "Error: Could not find representative in the irreducible BZ" << std::endl;
}

void IrreducibleBZ::init(std::string irreducibleBZFile)
{
    const double PI = 3.14159265358979323846;

    std::ifstream file(irreducibleBZFile);
    if (!file.is_open())
    {
        throw std::runtime_error("Error: Could not open file " + irreducibleBZFile);
    }

    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream ss(line);
        double x, y, z;

        if (!(ss >> x >> y >> z))
        {
            throw std::runtime_error("Error: Could not parse line " + line);
        }

        Eigen::Vector3d kVec(2 * pi * x, 2 * pi * y, 2 * pi * z); //= x * b_1 + y * b_2 + z * b_3; // //=// //
        irreducibleBZVectors.push_back(kVec);
    }

    file.close();
}

void SymmetrieGroup::init(std::string symmetryFile)
{
    std::ifstream file(symmetryFile);
    if (!file.is_open())
    {
        throw std::runtime_error("Error: Could not open file " + symmetryFile);
    }

    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream ss(line);
        Eigen::Matrix3d matrix;

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                if (!(ss >> matrix(i, j)))
                {
                    throw std::runtime_error("Error: Could not parse line " + line);
                }
            }
        }

        symmetrieOperations.push_back(matrix);
    }

    //// add the negative of the symmetrie matrices as well
    // std::vector<Eigen::Matrix3d> symmetrieOperationsMinus;
    // for (Eigen::Matrix3d symOp : symmetrieOperations)
    //{
    //     symmetrieOperationsMinus.push_back(-1 * symOp); // INVERS ODER EINFACH MINUS???
    // }
    //
    // symmetrieOperations.insert(symmetrieOperations.end(), symmetrieOperationsMinus.begin(), symmetrieOperationsMinus.end());

    file.close();
}

void IrreducibleBZ::initMagnonDisp(std::string couplingParameterFile)
{

    std::vector<Vector3D> path;

    for (Eigen::Vector3d vec : irreducibleBZVectors)
    {
        path.push_back(Vector3D(vec(0), vec(1), vec(2)));
    }

    std::vector<MagnonDispParam> magDisp = getMagneticDispersion(couplingParameterFile, path);

    this->magDisp = magDisp;
}

void IrreducibleBZ::initPhononDisp(std::string dynamicMatricesFile, std::string nextNeighbourFile)
{

    std::vector<Vector3D> path;

    for (Eigen::Vector3d vec : irreducibleBZVectors)
    {
        path.push_back(Vector3D(vec(0), vec(1), vec(2)));
    }

    std::vector<PhononDispParam> phDisp = getPhononDispersion(dynamicMatricesFile, nextNeighbourFile, path);

    this->phDisp = phDisp;
}

void IrreducibleBZ::readMultiplicities(const std::string &filename)
{

    std::vector<double> curr_multiplicities;

    std::ifstream file(filename);

    if (!file.is_open())
    {
        throw std::runtime_error("Error: Could not open file " + filename);
    }

    std::string line;
    while (std::getline(file, line))
    {
        int value = std::stoi(line);

        curr_multiplicities.push_back(value);
    }

    file.close();

    multiplicities = curr_multiplicities;

    if (multiplicities.size() != irreducibleBZVectors.size())
    {
        throw std::runtime_error("Error: Multiplicities do not match the number of irreducible BZ vectors");
    }

    std::cout << "Multiplicities read successfully" << std::endl;
}

Eigen::Vector3d IrreducibleBZ::findClosestVector(const Eigen::Vector3d &kVec)
{
    if (irreducibleBZVectors.empty())
    {
        throw std::runtime_error("Error: No vectors in the irreducible BZ");
    }

    double minDist = distance(kVec, irreducibleBZVectors[0]);
    Eigen::Vector3d closestVec = irreducibleBZVectors[0];

    for (const auto &vec : irreducibleBZVectors)
    {
        double dist = distance(kVec, vec);
        if (dist < minDist)
        {
            minDist = dist;
            closestVec = vec;
        }
    }

    return closestVec;
}
int IrreducibleBZ::findClosestVectorIndex(const Eigen::Vector3d &kVec)
{
    if (irreducibleBZVectors.empty())
    {
        throw std::runtime_error("Error: No vectors in the irreducible BZ");
    }

    double minDist = distance(kVec, irreducibleBZVectors[0]);
    int closestIndex = 0;

    for (int i = 1; i < irreducibleBZVectors.size(); i++)
    {
        double dist = distance(kVec, irreducibleBZVectors.at(i));
        if (dist < minDist)
        {
            minDist = dist;
            closestIndex = i;
        }
    }

    return closestIndex;
}

void IrreducibleBZ::initOccNumbers(double thermalEnergyPh, double thermalEnergyMag)
{
    phOccNumbers.resize(irreducibleBZVectors.size());
    magOccNumbers.resize(irreducibleBZVectors.size());

    for (int i = 0; i < irreducibleBZVectors.size(); i++)
    {
        phOccNumbers.at(i)[0] = BoseEinstein(phDisp.at(i).E[0], thermalEnergyPh);
        phOccNumbers.at(i)[1] = BoseEinstein(phDisp.at(i).E[1], thermalEnergyPh);
        phOccNumbers.at(i)[2] = BoseEinstein(phDisp.at(i).E[2], thermalEnergyPh);
        magOccNumbers.at(i) = BoseEinstein(magDisp.at(i).energy, thermalEnergyMag);
    }
    std::cout << "Occupation numbers initialized successfully" << std::endl;
}

void IrreducibleBZ::initMagnonPhononDispFromFile(std::string filepathPh, std::string filepathMag, std::string magnonDispOutputPath)
{
    std::vector<PhononDispParam> phData = readPhononDispParams(filepathPh);
    std::vector<MagnonDispParam> magData = getMagnonDispFromPhononDisp(phData, filepathMag, magnonDispOutputPath);

    this->phDisp = phData;
    this->magDisp = magData;
    std::cout << "Successfully initialized magnon and phonon data from file" << std::endl;
}

/*
Calculates the coupling coefficients \Gamma_{\pm} and \Gamma_{Z} for combinations of the given irreducible BZ vectors
- The units of the parameters are expected to be in meV/a.u.
- The units of the phonon dispersion are expected to be in meV
- The resulting coefficients are in units of meV
- The results are stored in the class variables gammaPlusGrid, gammaMinusGrid, and gammaZGrid (and their negative sign counterparts)
- The negative sign represents the negative sign of the second q vector in the coefficients index
*/
void IrreducibleBZ::initCoefficients(const std::vector<CouplingParameter> &parameters, int ftN)
{
    // Resize coefficients
    gammaPlusGrid.resize(irreducibleBZVectors.size());
    gammaMinusGrid.resize(irreducibleBZVectors.size());
    gammaZGrid.resize(irreducibleBZVectors.size());

    gammaPlusGrid_negativeSign.resize(irreducibleBZVectors.size());
    gammaMinusGrid_negativeSign.resize(irreducibleBZVectors.size());
    gammaZGrid_negativeSign.resize(irreducibleBZVectors.size());

    CGrid.resize(irreducibleBZVectors.size());
    DGrid.resize(irreducibleBZVectors.size());
    CGrid_negative_sign.resize(irreducibleBZVectors.size());
    DGrid_negative_sign.resize(irreducibleBZVectors.size());

    for (int i = 0; i < irreducibleBZVectors.size(); i++)
    {
        gammaPlusGrid.at(i).resize(irreducibleBZVectors.size());
        gammaMinusGrid.at(i).resize(irreducibleBZVectors.size());
        gammaZGrid.at(i).resize(irreducibleBZVectors.size());

        gammaPlusGrid_negativeSign.at(i).resize(irreducibleBZVectors.size());
        gammaMinusGrid_negativeSign.at(i).resize(irreducibleBZVectors.size());
        gammaZGrid_negativeSign.at(i).resize(irreducibleBZVectors.size());
    }

    std::complex<double> i(0, 1); // imaginary unit

#pragma omp parallel for
    for (int idx = 0; idx < irreducibleBZVectors.size(); idx++)
    {
        // first order in magnon variables
        {
            Eigen::Vector3d qVec = irreducibleBZVectors.at(idx);

            std::complex<double> C[] = {0, 0, 0};
            std::complex<double> D[] = {0, 0, 0};
            std::complex<double> C_negativeSign[] = {0, 0, 0};
            std::complex<double> D_negativeSign[] = {0, 0, 0};

            DMILikeCouplingParam D_k_values(qVec.x(), qVec.y(), qVec.z(), parameters);
            DMILikeCouplingParam D_k_values_minus(-qVec.x(), -qVec.y(), -qVec.z(), parameters);

            for (Axis ax : std::vector<Axis>{X, Y, Z})
            {
                std::complex<double> D_plus = D_k_values_minus.D[0][ax] + i * D_k_values_minus.D[1][ax];
                std::complex<double> D_minus = D_k_values.D[0][ax] - i * D_k_values.D[1][ax];

                std::complex<double> D_plus_negative_sign = D_k_values.D[0][ax] + i * D_k_values.D[1][ax];
                std::complex<double> D_minus_negative_sign = D_k_values_minus.D[0][ax] - i * D_k_values_minus.D[1][ax];

                for (int branch = 0; branch < 3; branch++)
                {
                    double pol_vec_component = phDisp.at(idx).polVectors[ax][branch];
                    C[branch] += 2.0 * i / sqrt(2 * S) * 3.8636 * D_minus * pol_vec_component * sqrt(1 / (2 * atomicMass * phDisp.at(idx).E[branch]));
                    D[branch] += -2.0 * i / sqrt(2 * S) * 3.8636 * D_plus * pol_vec_component * sqrt(1 / (2 * atomicMass * phDisp.at(idx).E[branch]));

                    C_negativeSign[branch] += 2.0 * i / sqrt(2 * S) * 3.8636 * D_minus_negative_sign * pol_vec_component * sqrt(1 / (2 * atomicMass * phDisp.at(idx).E[branch]));
                    D_negativeSign[branch] += -2.0 * i / sqrt(2 * S) * 3.8636 * D_plus_negative_sign * pol_vec_component * sqrt(1 / (2 * atomicMass * phDisp.at(idx).E[branch]));
                }
            }

            DGrid.at(idx) = {D[0], D[1], D[2]};
            CGrid.at(idx) = {C[0], C[1], C[2]};
            DGrid_negative_sign.at(idx) = {D_negativeSign[0], D_negativeSign[1], D_negativeSign[2]};
            CGrid_negative_sign.at(idx) = {C_negativeSign[0], C_negativeSign[1], C_negativeSign[2]};
        }

        // second order in magnon variables
        for (int j = 0; j < irreducibleBZVectors.size(); j++)
        {
            Eigen::Vector3d kVec = irreducibleBZVectors.at(idx);
            Eigen::Vector3d qVec = irreducibleBZVectors.at(j);

            std::complex<double> gammaPlus[] = {0, 0, 0};
            std::complex<double> gammaMinus[] = {0, 0, 0};
            std::complex<double> gammaZ[] = {0, 0, 0};

            std::complex<double> gammaPlus_negativeSign[] = {0, 0, 0};
            std::complex<double> gammaMinus_negativeSign[] = {0, 0, 0};
            std::complex<double> gammaZ_negativeSign[] = {0, 0, 0};

            int signs_plus_q[3] = {-1, 1, 1};
            int signs_minus_q[3] = {-1, 1, -1};

            Eigen::Vector3d gammaZ_k_prime_plus_q = irreducibleBZVectors.at(get_k_prime(kVec, qVec, signs_plus_q));
            Eigen::Vector3d gammaZ_k_prime_minus_q = irreducibleBZVectors.at(get_k_prime(kVec, -qVec, signs_minus_q));

            for (Axis ax : std::vector<Axis>{X, Y, Z})
            {
                // Coefficients for the gamma+, gamma-, and gammaZ terns with positive q
                std::complex<double> gammaMJxxmu = J_kq(kVec(0), kVec(1), kVec(2), qVec(0), qVec(1), qVec(2), parameters, X, X, ax);
                std::complex<double> gammaMJyymu = J_kq(kVec(0), kVec(1), kVec(2), qVec(0), qVec(1), qVec(2), parameters, Y, Y, ax);
                std::complex<double> gammaMJyxmu = J_kq(kVec(0), kVec(1), kVec(2), qVec(0), qVec(1), qVec(2), parameters, Y, X, ax);
                std::complex<double> gammaMJxymu = J_kq(kVec(0), kVec(1), kVec(2), qVec(0), qVec(1), qVec(2), parameters, X, Y, ax);

                std::complex<double> gammaPJxxmu = J_kq(-kVec(0), -kVec(1), -kVec(2), qVec(0), qVec(1), qVec(2), parameters, X, X, ax);
                std::complex<double> gammaPJyymu = J_kq(-kVec(0), -kVec(1), -kVec(2), qVec(0), qVec(1), qVec(2), parameters, Y, Y, ax);
                std::complex<double> gammaPJyxmu = J_kq(-kVec(0), -kVec(1), -kVec(2), qVec(0), qVec(1), qVec(2), parameters, Y, X, ax);
                std::complex<double> gammaPJxymu = J_kq(-kVec(0), -kVec(1), -kVec(2), qVec(0), qVec(1), qVec(2), parameters, X, Y, ax);

                std::complex<double> gammaZJxxmu = gammaPJxxmu;
                std::complex<double> gammaZJyymu = gammaPJyymu;
                std::complex<double> gammaZJyxmu = gammaPJyxmu;
                std::complex<double> gammaZJxymu = gammaPJxymu;
                std::complex<double> gammaZJzzmu = J_kq(-(kVec(0) - gammaZ_k_prime_plus_q(0)), -(kVec(1) - gammaZ_k_prime_plus_q(1)), -(kVec(2) - gammaZ_k_prime_plus_q(2)), qVec(0), qVec(1), qVec(2), parameters, Z, Z, ax);

                // Coefficients for the gamma+, gamma-, and gammaZ terns with negative q
                std::complex<double> gammaMJxxmu_negativeSign = J_kq(kVec(0), kVec(1), kVec(2), -qVec(0), -qVec(1), -qVec(2), parameters, X, X, ax);
                std::complex<double> gammaMJyymu_negativeSign = J_kq(kVec(0), kVec(1), kVec(2), -qVec(0), -qVec(1), -qVec(2), parameters, Y, Y, ax);
                std::complex<double> gammaMJyxmu_negativeSign = J_kq(kVec(0), kVec(1), kVec(2), -qVec(0), -qVec(1), -qVec(2), parameters, Y, X, ax);
                std::complex<double> gammaMJxymu_negativeSign = J_kq(kVec(0), kVec(1), kVec(2), -qVec(0), -qVec(1), -qVec(2), parameters, X, Y, ax);

                std::complex<double> gammaPJxxmu_negativeSign = J_kq(-kVec(0), -kVec(1), -kVec(2), -qVec(0), -qVec(1), -qVec(2), parameters, X, X, ax);
                std::complex<double> gammaPJyymu_negativeSign = J_kq(-kVec(0), -kVec(1), -kVec(2), -qVec(0), -qVec(1), -qVec(2), parameters, Y, Y, ax);
                std::complex<double> gammaPJyxmu_negativeSign = J_kq(-kVec(0), -kVec(1), -kVec(2), -qVec(0), -qVec(1), -qVec(2), parameters, Y, X, ax);
                std::complex<double> gammaPJxymu_negativeSign = J_kq(-kVec(0), -kVec(1), -kVec(2), -qVec(0), -qVec(1), -qVec(2), parameters, X, Y, ax);

                std::complex<double> gammaZJxxmu_negativeSign = gammaPJxxmu_negativeSign;
                std::complex<double> gammaZJyymu_negativeSign = gammaPJyymu_negativeSign;
                std::complex<double> gammaZJyxmu_negativeSign = gammaPJyxmu_negativeSign;
                std::complex<double> gammaZJxymu_negativeSign = gammaPJxymu_negativeSign;
                std::complex<double> gammaZJzzmu_negativeSign = J_kq(-(kVec(0) - gammaZ_k_prime_minus_q(0)), -(kVec(1) - gammaZ_k_prime_minus_q(1)), -(kVec(2) - gammaZ_k_prime_minus_q(2)), -qVec(0), -qVec(1), -qVec(2), parameters, Z, Z, ax);

                for (int branch = 0; branch < 3; branch++)
                {
                    // positive sign of q
                    gammaMinus[branch] += 3.8636 / (2.0 * S * sqrt(ftN)) * (gammaMJxxmu - gammaMJyymu - i * gammaMJxymu - i * gammaMJyxmu) * sqrt(1 / (2 * atomicMass * phDisp.at(j).E[branch])) * phDisp.at(j).polVectors[ax][branch];
                    gammaPlus[branch] += 3.8636 / (2.0 * S * sqrt(ftN)) * (gammaPJxxmu - gammaPJyymu + i * gammaPJxymu + i * gammaPJyxmu) * sqrt(1 / (2 * atomicMass * phDisp.at(j).E[branch])) * phDisp.at(j).polVectors[ax][branch];
                    gammaZ[branch] += 3.8636 / (2.0 * S * sqrt(ftN)) * (2.0 * gammaZJxxmu + 2.0 * gammaZJyymu - 4.0 * gammaZJzzmu + 2.0 * i * gammaZJyxmu - 2.0 * i * gammaZJxymu) * sqrt(1 / (2 * atomicMass * phDisp.at(j).E[branch])) * phDisp.at(j).polVectors[ax][branch];

                    // debug statement
                    // if (gammaZ[branch].real() > 1E4 || gammaZ[branch].imag() > 1E4)
                    //{
                    //    std::cout << gammaZ[branch] << std::endl;
                    //    std::cout << "Jxxmu: " << Jxxmu << " Jyymu " << Jyymu << "Jyxmu " << Jyxmu << "Jxymu " << Jxymu << "Jzzmu " << Jzzmu << std::endl;
                    //}

                    // neg sign of q
                    gammaMinus_negativeSign[branch] += 3.8636 / (2.0 * S * sqrt(ftN)) * (gammaMJxxmu_negativeSign - gammaMJyymu_negativeSign - i * gammaMJxymu_negativeSign - i * gammaMJyxmu_negativeSign) * sqrt(1 / (2 * atomicMass * phDisp.at(j).E[branch])) * phDisp.at(j).polVectors[ax][branch];
                    gammaPlus_negativeSign[branch] += 3.8636 / (2.0 * S * sqrt(ftN)) * (gammaPJxxmu_negativeSign - gammaPJyymu_negativeSign + i * gammaPJxymu_negativeSign + i * gammaPJyxmu_negativeSign) * sqrt(1 / (2 * atomicMass * phDisp.at(j).E[branch])) * phDisp.at(j).polVectors[ax][branch];
                    gammaZ_negativeSign[branch] += 3.8636 / (2.0 * S * sqrt(ftN)) * (2.0 * gammaZJxxmu_negativeSign + 2.0 * gammaZJyymu_negativeSign - 4.0 * gammaZJzzmu_negativeSign + 2.0 * i * gammaZJyxmu_negativeSign - 2.0 * i * gammaZJxymu_negativeSign) * sqrt(1 / (2 * atomicMass * phDisp.at(j).E[branch])) * phDisp.at(j).polVectors[ax][branch];
                }
            }

            for (int branch = 0; branch < 3; branch++)
            {
                gammaMinusGrid.at(idx).at(j)[branch] = gammaMinus[branch];
                gammaPlusGrid.at(idx).at(j)[branch] = gammaPlus[branch];
                gammaZGrid.at(idx).at(j)[branch] = gammaZ[branch];

                gammaMinusGrid_negativeSign.at(idx).at(j)[branch] = gammaMinus_negativeSign[branch];
                gammaZGrid_negativeSign.at(idx).at(j)[branch] = gammaZ_negativeSign[branch];
                gammaPlusGrid_negativeSign.at(idx).at(j)[branch] = gammaPlus_negativeSign[branch];
            }
        }
    }
}

// Find the reciprocal lattice vector G (in which basis?) that satisfies the condition G = -k + k' + q
// Apply the symmetries to the k' vector and possibly also the k and q vectors???
Eigen::Vector3d IrreducibleBZ::getG_Z(const Eigen::Vector3d &kVec, const Eigen::Vector3d &qVec)
{
    int indices[7] = {0, 1, -1, 2, -2};

    for (auto k_prime : irreducibleBZVectors)
    {
        for (auto k_p_trans : applySymmOp(k_prime))
        {
            for (int i : indices)
            {
                for (int j : indices)
                {
                    for (int k : indices)
                    {
                        Eigen::Vector3d G = i * Eigen::Vector3d(2 * pi, 0, 0) + j * Eigen::Vector3d(0, 2 * pi, 0) + k * Eigen::Vector3d(0, 0, 2 * pi);
                        if (distance(G, k_p_trans - kVec + qVec) < 0.0001)
                        {
                            // std::cout << "G " << G(0) / (2 * pi) << " " << G(1) / (2 * pi) << " " << G(2) / (2 * pi) << " | k_prime: " << k_prime(0) / (2 * pi) << " " << k_prime(1) / (2 * pi) << " " << k_prime(2) / (2 * pi) << std::endl;

                            return G;
                        }
                    }
                }
            }
        }
    }

    std::cout << "Error: Could not find GZ" << std::endl;
    return Eigen::Vector3d(0, 0, 0);
}

int IrreducibleBZ::get_k_prime(const Eigen::Vector3d &kVec, const Eigen::Vector3d &qVec, int signs[3])
{
    int indices[7] = {0, 1, -1, 2, -2};

    for (int n = 0; n < irreducibleBZVectors.size(); n++)
    {
        for (auto k_p_trans : applySymmOp(irreducibleBZVectors.at(n)))
        {
            for (int i : indices)
            {
                for (int j : indices)
                {
                    for (int k : indices)
                    {

                        Eigen::Vector3d G = i * Eigen::Vector3d(2 * pi, 0, 0) + j * Eigen::Vector3d(0, 2 * pi, 0) + k * Eigen::Vector3d(0, 0, 2 * pi);
                        if (distance(G, signs[0] * kVec + signs[1] * k_p_trans - signs[2] * qVec) < 0.0001)
                        {
                            // std::cout << "G " << G(0) / (2 * pi) << " " << G(1) / (2 * pi) << " " << G(2) / (2 * pi) << " | k_prime: " << k_p_trans(0) / (2 * pi) << " " << k_p_trans(1) / (2 * pi) << " " << k_p_trans(2) / (2 * pi) << std::endl;
                            return n;
                        }
                    }
                }
            }
        }
    }
    std::cout << "Error: Could not find k_prime" << std::endl;
    return -1;
}

// Find the reciprocal lattice vector G (in which basis?) that satisfies the condition G = -k - k' + q
// Apply the symmetries to the k' vector and possibly also the k and q vectors???
Eigen::Vector3d IrreducibleBZ::getG_gammaP(const Eigen::Vector3d &kVec, const Eigen::Vector3d &qVec)
{
    int indices[7] = {0, 1, -1, 2, -2};

    for (auto k_prime : irreducibleBZVectors)
    {
        for (auto k_p_trans : applySymmOp(k_prime))
        {
            for (int i : indices)
            {
                for (int j : indices)
                {
                    for (int k : indices)
                    {
                        Eigen::Vector3d G = i * Eigen::Vector3d(2 * pi, 0, 0) + j * Eigen::Vector3d(0, 2 * pi, 0) + k * Eigen::Vector3d(0, 0, 2 * pi);
                        if (distance(G, -k_p_trans - kVec + qVec) < 0.01)
                        {
                            // std::cout << "G " << G(0) / (2 * pi) << " " << G(1) / (2 * pi) << " " << G(2) / (2 * pi) << " | k_prime: " << k_prime(0) / (2 * pi) << " " << k_prime(1) / (2 * pi) << " " << k_prime(2) / (2 * pi) << std::endl;

                            return G;
                        }
                    }
                }
            }
        }
    }

    std::cout << "Error: Could not find GP" << std::endl;
    return Eigen::Vector3d(0, 0, 0);
}

// Find the reciprocal lattice vector G (in which basis?) that satisfies the condition G = k + k' + q
// Apply the symmetries to the k' vector and possibly also the k and q vectors???
Eigen::Vector3d IrreducibleBZ::getG_gammaM(const Eigen::Vector3d &kVec, const Eigen::Vector3d &qVec)
{
    int indices[7] = {0, 1, -1, 2, -2};

    for (auto k_prime : irreducibleBZVectors)
    {
        for (auto k_p_trans : applySymmOp(k_prime))
        {
            for (int i : indices)
            {
                for (int j : indices)
                {
                    for (int k : indices)
                    {
                        Eigen::Vector3d G = i * Eigen::Vector3d(2 * pi, 0, 0) + j * Eigen::Vector3d(0, 2 * pi, 0) + k * Eigen::Vector3d(0, 0, 2 * pi);
                        if (distance(G, k_p_trans + kVec + qVec) < 0.001)
                        {
                            // std::cout << "G " << G(0) / (2 * pi) << " " << G(1) / (2 * pi) << " " << G(2) / (2 * pi) << " | k_prime: " << k_prime(0) / (2 * pi) << " " << k_prime(1) / (2 * pi) << " " << k_prime(2) / (2 * pi) << std::endl;

                            return G;
                        }
                    }
                }
            }
        }
    }

    std::cout << "Error: Could not find GM" << std::endl;
    std::cout << "k: " << kVec(0) / (2 * pi) << " " << kVec(1) / (2 * pi) << " " << kVec(2) / (2 * pi) << " | q: " << qVec(0) / (2 * pi) << " " << qVec(1) / (2 * pi) << " " << qVec(2) / (2 * pi) << std::endl;
    std::cout << "k+q: " << (kVec + qVec)(0) / (2 * pi) << " " << (kVec + qVec)(1) / (2 * pi) << " " << (kVec + qVec)(2) / (2 * pi) << std::endl;
    return Eigen::Vector3d(0, 0, 0);
}

// Saves the coefficients of the magnon-phonon coupling to a file
// Absolute value squared of the coefficients are saved
void IrreducibleBZ::saveCoefficientsAsSqrtAbs(std::string filename)
{

    std::ofstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Error: Could not open file " + filename);
    }

    for (int i = 0; i < irreducibleBZVectors.size(); i++)
    {
        for (int j = 0; j < irreducibleBZVectors.size(); j++)
        {
            file << i << " " << j
                 << " " << std::norm(gammaPlusGrid.at(i).at(j)[0]) << " " << std::norm(gammaPlusGrid.at(i).at(j)[1]) << " " << std::norm(gammaPlusGrid.at(i).at(j)[2])
                 << " " << std::norm(gammaMinusGrid.at(i).at(j)[0]) << " " << std::norm(gammaMinusGrid.at(i).at(j)[1]) << " " << std::norm(gammaMinusGrid.at(i).at(j)[2])
                 << " " << std::norm(gammaZGrid.at(i).at(j)[0]) << " " << std::norm(gammaZGrid.at(i).at(j)[1]) << " " << std::norm(gammaZGrid.at(i).at(j)[2])
                 << " " << std::norm(gammaPlusGrid_negativeSign.at(i).at(j)[0]) << " " << std::norm(gammaPlusGrid_negativeSign.at(i).at(j)[1]) << " " << std::norm(gammaPlusGrid_negativeSign.at(i).at(j)[2])
                 << " " << std::norm(gammaMinusGrid_negativeSign.at(i).at(j)[0]) << " " << std::norm(gammaMinusGrid_negativeSign.at(i).at(j)[1]) << " " << std::norm(gammaMinusGrid_negativeSign.at(i).at(j)[2])
                 << " " << std::norm(gammaZGrid_negativeSign.at(i).at(j)[0]) << " " << std::norm(gammaZGrid_negativeSign.at(i).at(j)[1]) << " " << std::norm(gammaZGrid_negativeSign.at(i).at(j)[2])
                 << " " << std::norm(CGrid.at(i)[0]) << " " << std::norm(CGrid.at(i)[1]) << " " << std::norm(CGrid.at(i)[2])
                 << " " << std::norm(DGrid.at(i)[0]) << " " << std::norm(DGrid.at(i)[1]) << " " << std::norm(DGrid.at(i)[2])
                 << " " << std::norm(CGrid_negative_sign.at(i)[0]) << " " << std::norm(CGrid_negative_sign.at(i)[1]) << " " << std::norm(CGrid_negative_sign.at(i)[2])
                 << " " << std::norm(DGrid_negative_sign.at(i)[0]) << " " << std::norm(DGrid_negative_sign.at(i)[1]) << " " << std::norm(DGrid_negative_sign.at(i)[2])
                 << "\n";
        }
    }

    file.close();
}

void IrreducibleBZ::readCoefficients(std::string filename)
{

    gammaPlusGrid.resize(irreducibleBZVectors.size());
    gammaMinusGrid.resize(irreducibleBZVectors.size());
    gammaZGrid.resize(irreducibleBZVectors.size());

    gammaPlusGrid_negativeSign.resize(irreducibleBZVectors.size());
    gammaMinusGrid_negativeSign.resize(irreducibleBZVectors.size());
    gammaZGrid_negativeSign.resize(irreducibleBZVectors.size());

    CGrid.resize(irreducibleBZVectors.size());
    DGrid.resize(irreducibleBZVectors.size());
    CGrid_negative_sign.resize(irreducibleBZVectors.size());
    DGrid_negative_sign.resize(irreducibleBZVectors.size());

    for (int i = 0; i < irreducibleBZVectors.size(); i++)
    {
        gammaPlusGrid.at(i).resize(irreducibleBZVectors.size());
        gammaMinusGrid.at(i).resize(irreducibleBZVectors.size());
        gammaZGrid.at(i).resize(irreducibleBZVectors.size());

        gammaPlusGrid_negativeSign.at(i).resize(irreducibleBZVectors.size());
        gammaMinusGrid_negativeSign.at(i).resize(irreducibleBZVectors.size());
        gammaZGrid_negativeSign.at(i).resize(irreducibleBZVectors.size());
    }

    std::ifstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Error: Could not open file " + filename);
    }

    // Read the coefficients from the file and initialize the coefficient grid
    // Assuming the file format is: index1 index2 gammaM gammaP gammaZ

    int index1, index2;
    double gammaM[3], gammaP[3], gammaZ[3], gammaM_negativeSign[3], gammaP_negativeSign[3], gammaZ_negativeSign[3], C[3], D[3], C_negativeSign[3], D_negativeSign[3];

    while (file >> index1 >> index2 >> gammaP[0] >> gammaP[1] >> gammaP[2] >> gammaM[0] >> gammaM[1] >> gammaM[2] >> gammaZ[0] >> gammaZ[1] >> gammaZ[2] >> gammaP_negativeSign[0] >> gammaP_negativeSign[1] >> gammaP_negativeSign[2] >> gammaM_negativeSign[0] >> gammaM_negativeSign[1] >> gammaM_negativeSign[2] >> gammaZ_negativeSign[0] >> gammaZ_negativeSign[1] >> gammaZ_negativeSign[2] >> C[0] >> C[1] >> C[2] >> D[0] >> D[1] >> D[2] >> C_negativeSign[0] >> C_negativeSign[1] >> C_negativeSign[2] >> D_negativeSign[0] >> D_negativeSign[1] >> D_negativeSign[2])
    {
        gammaMinusGrid.at(index1).at(index2)[0] = gammaM[0];
        gammaMinusGrid.at(index1).at(index2)[1] = gammaM[1];
        gammaMinusGrid.at(index1).at(index2)[2] = gammaM[2];

        gammaPlusGrid.at(index1).at(index2)[0] = gammaP[0];
        gammaPlusGrid.at(index1).at(index2)[1] = gammaP[1];
        gammaPlusGrid.at(index1).at(index2)[2] = gammaP[2];

        gammaZGrid.at(index1).at(index2)[0] = gammaZ[0];
        gammaZGrid.at(index1).at(index2)[1] = gammaZ[1];
        gammaZGrid.at(index1).at(index2)[2] = gammaZ[2];

        gammaMinusGrid_negativeSign.at(index1).at(index2)[0] = gammaM_negativeSign[0];
        gammaMinusGrid_negativeSign.at(index1).at(index2)[1] = gammaM_negativeSign[1];
        gammaMinusGrid_negativeSign.at(index1).at(index2)[2] = gammaM_negativeSign[2];

        gammaPlusGrid_negativeSign.at(index1).at(index2)[0] = gammaP_negativeSign[0];
        gammaPlusGrid_negativeSign.at(index1).at(index2)[1] = gammaP_negativeSign[1];
        gammaPlusGrid_negativeSign.at(index1).at(index2)[2] = gammaP_negativeSign[2];

        gammaZGrid_negativeSign.at(index1).at(index2)[0] = gammaZ_negativeSign[0];
        gammaZGrid_negativeSign.at(index1).at(index2)[1] = gammaZ_negativeSign[1];
        gammaZGrid_negativeSign.at(index1).at(index2)[2] = gammaZ_negativeSign[2];

        CGrid.at(index1)[0] = C[0];
        CGrid.at(index1)[1] = C[1];
        CGrid.at(index1)[2] = C[2];

        DGrid.at(index1)[0] = D[0];
        DGrid.at(index1)[1] = D[1];
        DGrid.at(index1)[2] = D[2];

        CGrid_negative_sign.at(index1)[0] = C_negativeSign[0];
        CGrid_negative_sign.at(index1)[1] = C_negativeSign[1];
        CGrid_negative_sign.at(index1)[2] = C_negativeSign[2];

        DGrid_negative_sign.at(index1)[0] = D_negativeSign[0];
        DGrid_negative_sign.at(index1)[1] = D_negativeSign[1];
        DGrid_negative_sign.at(index1)[2] = D_negativeSign[2];
    }

    file.close();
    std::cout << "Successfully read and initialized coefficients from file " << filename << std::endl;
}

double deltaDistrApprox(double x)
{
    double a = 1; // 2 < a <  a=5 -> way to much, a=2.2 -> basically always zero, a=3 -> energy escalates slowly

    if (std::abs(x) < a)
    {
        // std::cout << "x: " << x << std::endl;
    }
    return 1.0 / pi * exp(-(x / a) * (x / a)); /// (a * sqrt(pi))
}

void IrreducibleBZ::integrate()
{

    double tMax = 1E1;
    double dt = 1E-5;
    int nMax = int(tMax / dt);

    std::ofstream energy_file("Outputs/time_evolut_energies.txt");
    std::ofstream temp_file("Outputs/time_evolut_temperature.txt");

    // init the occupation for the next current step
    std::vector<std::array<double, 3>> occNumPh_curr = phOccNumbers;
    std::vector<double> occNumMag_curr = magOccNumbers;

    // Main integration loop over time
    for (int n = 0; n < nMax; n++)
    {
        // Loop over all vectors in the irreducible BZ
        // #pragma omp parallel for

        double sumPh[3] = {0, 0, 0};
        double sumMag = 0;

        for (int vec_outer_idx = 0; vec_outer_idx < irreducibleBZVectors.size(); vec_outer_idx++)
        {
            // std::cout << magOccNumbers.at(vec_outer_idx) << std::endl;

            //  Calculate sum
            for (int branch = 0; branch < 3; branch++)
            {

                // FIRST ORDER IN MAGNON VARIABLES

                // M^+(0,-q) = D(q)
                // M^+(0,q) = D(-q)
                // M^-(0,-q) = C(-q)
                // M^-(0,q) = C(q)

                // phonon
                {
                    double N_mk = occNumMag_curr.at(vec_outer_idx);
                    double n_k = occNumPh_curr.at(vec_outer_idx)[branch];
                    double deltaDistr = deltaDistrApprox(phDisp.at(vec_outer_idx).E[branch] - magDisp.at(vec_outer_idx).energy);
                    double M_minus_negSign = CGrid_negative_sign.at(vec_outer_idx)[branch].real();
                    double M_plus_posSign = DGrid_negative_sign.at(vec_outer_idx)[branch].real();
                    sumPh[branch] += 2 * pi * deltaDistr * (M_minus_negSign * N_mk * (n_k + 1) - M_plus_posSign * (N_mk + 1) * n_k);
                }
                // magnon
                {
                    double N_k = occNumMag_curr.at(vec_outer_idx);
                    double n_mk = occNumPh_curr.at(vec_outer_idx)[branch];
                    double deltaDistr = deltaDistrApprox(phDisp.at(vec_outer_idx).E[branch] - magDisp.at(vec_outer_idx).energy);
                    double M_plus_negSign = DGrid.at(vec_outer_idx)[branch].real();
                    double M_minus_posSign = CGrid.at(vec_outer_idx)[branch].real();
                    sumMag += 2 * pi * deltaDistr * (M_plus_negSign * (N_k + 1) * n_mk - M_minus_posSign * N_k * (n_mk + 1));
                }

                // SECOND ORDER IN MAGNON VARIABLES
                // INNER LOOP over all vectors in the irreducible BZ
                for (int vec_inner_idx = 0; vec_inner_idx < irreducibleBZVectors.size(); vec_inner_idx++)
                {
                    // multiplicities for the inner loop over the irr BZ
                    double m = multiplicities.at(vec_inner_idx);

                    // phonon
                    {
                        // handle terms with gamma minus and -q
                        {
                            int rep_index_k_prime = k_prime_representatives_gammaM_plus_q[vec_inner_idx][vec_outer_idx];
                            double N_k_prime = occNumMag_curr.at(rep_index_k_prime);
                            double deltaDistr = deltaDistrApprox(phDisp.at(vec_outer_idx).E[branch] - magDisp.at(vec_inner_idx).energy - magDisp.at(rep_index_k_prime).energy);
                            double n_q = occNumPh_curr.at(vec_outer_idx)[branch];
                            double N_k = occNumMag_curr.at(vec_inner_idx);
                            double m_k_prime = multiplicities.at(rep_index_k_prime);

                            sumPh[branch] += 2 * pi * m * m_k_prime * (n_q + 1) * gammaMinusGrid_negativeSign[vec_inner_idx][vec_outer_idx][branch].real() * N_k * N_k_prime * deltaDistr;
                        }
                        // handle terms with gamma z and -q
                        {
                            int rep_index_k_prime = k_prime_representatives_gammaZ_minus_q[vec_inner_idx][vec_outer_idx];
                            double N_k_prime = occNumMag_curr.at(rep_index_k_prime);
                            double deltaDistr = deltaDistrApprox(phDisp.at(vec_outer_idx).E[branch] + magDisp.at(vec_inner_idx).energy - magDisp.at(rep_index_k_prime).energy);
                            double n_q = occNumPh_curr.at(vec_outer_idx)[branch];
                            double N_k = occNumMag_curr.at(vec_inner_idx);
                            double m_k_prime = multiplicities.at(rep_index_k_prime);

                            sumPh[branch] += 2 * pi * m * m_k_prime * (n_q + 1) * gammaZGrid_negativeSign[vec_inner_idx][vec_outer_idx][branch].real() * (N_k + 1) * N_k_prime * deltaDistr;
                        }
                        // handle terms with gamma P and +q
                        {
                            int rep_index_k_prime = k_prime_representatives_gammaP_plus_q[vec_inner_idx][vec_outer_idx];
                            double N_k_prime = occNumMag_curr.at(rep_index_k_prime);
                            double deltaDistr = deltaDistrApprox(-phDisp.at(vec_outer_idx).E[branch] + magDisp.at(vec_inner_idx).energy + magDisp.at(rep_index_k_prime).energy);
                            double n_q = occNumPh_curr.at(vec_outer_idx)[branch];
                            double N_k = occNumMag_curr.at(vec_inner_idx);
                            double m_k_prime = multiplicities.at(rep_index_k_prime);

                            sumPh[branch] -= 2 * pi * m * m_k_prime * n_q * gammaPlusGrid[vec_inner_idx][vec_outer_idx][branch].real() * (N_k + 1) * (N_k_prime + 1) * deltaDistr;
                        }

                        // handle terms with gamma Z and +q
                        {
                            int rep_index_k_prime = k_prime_representatives_gammaZ_minus_q[vec_inner_idx][vec_outer_idx];
                            double N_k_prime = occNumMag_curr.at(rep_index_k_prime);
                            double n_k = occNumPh_curr.at(vec_outer_idx)[branch];
                            double N_k = occNumMag_curr.at(vec_inner_idx);
                            double m_k_prime = multiplicities.at(rep_index_k_prime);

                            double deltaDistrGammaZ = deltaDistrApprox(-phDisp.at(vec_outer_idx).E[branch] + magDisp.at(vec_inner_idx).energy - magDisp.at(rep_index_k_prime).energy);

                            sumPh[branch] -= 2 * pi * m * m_k_prime * n_k * gammaZGrid[vec_inner_idx][vec_outer_idx][branch].real() * (N_k + 1) * N_k_prime * deltaDistrGammaZ;
                        }
                    }

                    // Magnon
                    {
                        // positive gammaP term
                        {
                            int rep_index_k_prime = k_prime_representatives_gammaP_plus_q[vec_outer_idx][vec_inner_idx];
                            double N_k_prime = occNumMag_curr.at(rep_index_k_prime);
                            double deltaDistr = deltaDistrApprox(magDisp.at(vec_outer_idx).energy + magDisp.at(rep_index_k_prime).energy - phDisp.at(vec_inner_idx).E[branch]);
                            double N_k = occNumMag_curr.at(vec_outer_idx);
                            double n_q = occNumPh_curr.at(vec_inner_idx)[branch];
                            double m_k_prime = multiplicities.at(rep_index_k_prime);

                            sumMag += 2 * pi * m * m_k_prime * (N_k + 1) * 2 * gammaPlusGrid[vec_outer_idx][vec_inner_idx][branch].real() * (N_k_prime + 1) * n_q * deltaDistr;
                        }
                        // positive gammaZ term
                        {
                            int rep_index_k_prime = k_prime_representatives_gammaZ_plus_q[vec_outer_idx][vec_inner_idx];
                            double N_k_prime = occNumMag_curr.at(rep_index_k_prime);
                            double deltaDistr_pmm = deltaDistrApprox(magDisp.at(vec_outer_idx).energy - magDisp.at(rep_index_k_prime).energy - phDisp.at(vec_inner_idx).E[branch]);
                            double deltaDistr_pmp = deltaDistrApprox(magDisp.at(vec_outer_idx).energy - magDisp.at(rep_index_k_prime).energy + phDisp.at(vec_inner_idx).E[branch]);
                            double n_qv = occNumPh_curr.at(vec_inner_idx)[branch];
                            double N_k = occNumMag_curr.at(vec_outer_idx);
                            double n_mqv = n_qv; // phOccNumbers.at(findRepresentative(-irreducibleBZVectors.at(vec_inner_idx)))[branch];
                            double m_k_prime = multiplicities.at(rep_index_k_prime);

                            sumMag += 2 * pi * m * m_k_prime * (N_k + 1) * N_k_prime * (n_qv * deltaDistr_pmm + (n_mqv + 1) * deltaDistr_pmp) * gammaZGrid[vec_outer_idx][vec_inner_idx][branch].real();
                        }
                        // negative gammaM term
                        {
                            int rep_index_k_prime = k_prime_representatives_gammaM_plus_q[vec_outer_idx][vec_inner_idx];
                            double N_k_prime = magOccNumbers.at(rep_index_k_prime);
                            double n_mqv = occNumPh_curr.at(vec_inner_idx)[branch]; // phOccNumbers.at(findRepresentative(-irreducibleBZVectors.at(vec_inner_idx)))[branch];
                            double deltaDistr = deltaDistrApprox(-magDisp.at(vec_outer_idx).energy - magDisp.at(rep_index_k_prime).energy + phDisp.at(vec_inner_idx).E[branch]);
                            double N_k = occNumMag_curr.at(vec_outer_idx);
                            // std::cout << rep_index_k_prime << " " << N_k_prime << " " << n_mqv << " " << deltaDistr << gammaMinusGrid[vec_outer_idx][vec_inner_idx][branch].real() << " " << m << " " << occNumMag_curr << " " << m << "\n " << std::endl;
                            double m_k_prime = multiplicities.at(rep_index_k_prime);

                            sumMag -= 2 * pi * m * m_k_prime * N_k * 2 * gammaMinusGrid[vec_outer_idx][vec_inner_idx][branch].real() * N_k_prime * (n_mqv + 1) * deltaDistr;
                        }
                        // negative gammaZ term
                        {
                            int rep_index_k_prime = k_prime_representatives_gammaZ_plus_q[vec_outer_idx][vec_inner_idx];
                            double N_k_prime = occNumMag_curr.at(rep_index_k_prime);
                            double deltaDistr_mpm = deltaDistrApprox(-magDisp.at(vec_outer_idx).energy + magDisp.at(rep_index_k_prime).energy - phDisp.at(vec_inner_idx).E[branch]);
                            double deltaDistr_mpp = deltaDistrApprox(-magDisp.at(vec_outer_idx).energy + magDisp.at(rep_index_k_prime).energy + phDisp.at(vec_inner_idx).E[branch]);
                            double n_qv = occNumPh_curr.at(vec_inner_idx)[branch];
                            double n_mqv = n_qv; // phOccNumbers.at(findRepresentative(-irreducibleBZVectors.at(vec_inner_idx)))[branch];
                            double N_k = occNumMag_curr.at(vec_outer_idx);
                            double m_k_prime = multiplicities.at(rep_index_k_prime);

                            sumMag -= 2 * pi * m * m_k_prime * N_k * gammaZGrid[vec_outer_idx][vec_inner_idx][branch].real() * (N_k_prime + 1) * (n_qv * deltaDistr_mpm + (n_mqv + 1) * deltaDistr_mpp);
                        }
                    }
                }
            }
        }

        // std::cout << "sumMag " << sumMag << std::endl;
        // std::cout << "sumPh " << sumPh[0] << " " << sumPh[1] << " " << sumPh[2] << std::endl;

        // calculate the total energy
        double magEnergy = 0;
        double phEnergy = 0;
        double totalEnergy = 0;
        for (int i = 0; i < irreducibleBZVectors.size(); i++)
        {
            magEnergy += multiplicities.at(i) * occNumMag_curr.at(i) * magDisp.at(i).energy;
            phEnergy += multiplicities.at(i) * (occNumPh_curr.at(i)[0] * phDisp.at(i).E[0] + occNumPh_curr.at(i)[1] * phDisp.at(i).E[1] + occNumPh_curr.at(i)[2] * phDisp.at(i).E[2]);
        }
        totalEnergy = magEnergy + phEnergy;

        // Print progress
        if (n % 1000 == 0)
        {
            std::cout << "Progress: " << double(n) / nMax << std::endl;
        }

        if (n % 100 == 0)
        {
            // write energy to file
            energy_file << n * dt << " " << magEnergy << " " << phEnergy << " " << totalEnergy << " ";

            for (int i = 0; i < irreducibleBZVectors.size(); i++)
            {
                energy_file << occNumMag_curr.at(i) << " " << occNumPh_curr.at(i)[0] << " " << occNumPh_curr.at(i)[1] << " " << occNumPh_curr.at(i)[2] << " ";
            }

            energy_file << "\n";

            // calculate the temperature and write to file
            temp_file << n * dt;

            // temperature phonon
            for (int i = 0; i < irreducibleBZVectors.size(); i++)
            {
                double sum = 0;
                for (int branch : {0, 1, 2})
                {
                    double thermalEnergyPh = BoseEinsteinInverted(phDisp.at(i).E[branch], occNumPh_curr.at(i)[branch]);
                    temp_file << " " << thermalEnergyPh;
                }
            }

            // temperature magnon
            for (int i = 0; i < irreducibleBZVectors.size(); i++)
            {
                double thermalEnergyMag = BoseEinsteinInverted(magDisp.at(i).energy, occNumMag_curr.at(i));
                temp_file << " " << thermalEnergyMag;
            }

            temp_file << "\n";
        }

        // Finally, update the occupation numbers
        for (int i = 0; i < irreducibleBZVectors.size(); i++)
        {
            occNumMag_curr.at(i) = occNumMag_curr.at(i) + dt * sumMag;
            for (int branch : {0, 1, 2})
            {
                occNumPh_curr.at(i)[branch] = occNumPh_curr.at(i)[branch] + dt * sumPh[branch];
            }
        }
    }
    temp_file.close();
    energy_file.close();
}

void IrreducibleBZ::init_k_prime()
{

    /*
        // init all possible rec lattice vectors
    std::vector<std::vector<Eigen::Vector3d>> G_gammaM_minus_q = std::vector(irreducibleBZVectors.size(), std::vector<Eigen::Vector3d>(irreducibleBZVectors.size()));
    std::vector<std::vector<Eigen::Vector3d>> G_gammaZ_minus_q = std::vector(irreducibleBZVectors.size(), std::vector<Eigen::Vector3d>(irreducibleBZVectors.size()));
    std::vector<std::vector<Eigen::Vector3d>> G_gammaP_minus_q = std::vector(irreducibleBZVectors.size(), std::vector<Eigen::Vector3d>(irreducibleBZVectors.size()));

    std::vector<std::vector<Eigen::Vector3d>> G_gammaM_plus_q = std::vector(irreducibleBZVectors.size(), std::vector<Eigen::Vector3d>(irreducibleBZVectors.size()));
    std::vector<std::vector<Eigen::Vector3d>> G_gammaZ_plus_q = std::vector(irreducibleBZVectors.size(), std::vector<Eigen::Vector3d>(irreducibleBZVectors.size()));
    std::vector<std::vector<Eigen::Vector3d>> G_gammaP_plus_q = std::vector(irreducibleBZVectors.size(), std::vector<Eigen::Vector3d>(irreducibleBZVectors.size()));

#pragma omp parallel for
    for (int i = 0; i < irreducibleBZVectors.size(); i++)
    {
        for (int j = 0; j < irreducibleBZVectors.size(); j++)
        {
            G_gammaM_minus_q[i][j] = getG_gammaM(irreducibleBZVectors.at(i), -irreducibleBZVectors.at(j));
            G_gammaZ_minus_q[i][j] = getG_Z(irreducibleBZVectors.at(i), -irreducibleBZVectors.at(j));
            G_gammaP_minus_q[i][j] = getG_gammaP(irreducibleBZVectors.at(i), -irreducibleBZVectors.at(j));

            G_gammaM_plus_q[i][j] = getG_gammaM(irreducibleBZVectors.at(i), irreducibleBZVectors.at(j));
            G_gammaZ_plus_q[i][j] = getG_Z(irreducibleBZVectors.at(i), irreducibleBZVectors.at(j));
            G_gammaP_plus_q[i][j] = getG_gammaP(irreducibleBZVectors.at(i), irreducibleBZVectors.at(j));
        }
    }

    std::vector<std::vector<Eigen::Vector3d>> k_prime_gammaZ_minus_q, k_prime_gammaM_minus_q, k_prime_gammaP_minus_q, k_prime_gammaZ_plus_q, k_prime_gammaM_plus_q, k_prime_gammaP_plus_q;

    k_prime_gammaM_minus_q = std::vector(irreducibleBZVectors.size(), std::vector<Eigen::Vector3d>(irreducibleBZVectors.size()));
    k_prime_gammaZ_minus_q = std::vector(irreducibleBZVectors.size(), std::vector<Eigen::Vector3d>(irreducibleBZVectors.size()));
    k_prime_gammaP_minus_q = std::vector(irreducibleBZVectors.size(), std::vector<Eigen::Vector3d>(irreducibleBZVectors.size()));

    k_prime_gammaM_plus_q = std::vector(irreducibleBZVectors.size(), std::vector<Eigen::Vector3d>(irreducibleBZVectors.size()));
    k_prime_gammaZ_plus_q = std::vector(irreducibleBZVectors.size(), std::vector<Eigen::Vector3d>(irreducibleBZVectors.size()));
    k_prime_gammaP_plus_q = std::vector(irreducibleBZVectors.size(), std::vector<Eigen::Vector3d>(irreducibleBZVectors.size()));

    k_prime_representatives_gammaM_minus_q = std::vector(irreducibleBZVectors.size(), std::vector<int>(irreducibleBZVectors.size()));
    k_prime_representatives_gammaZ_minus_q = std::vector(irreducibleBZVectors.size(), std::vector<int>(irreducibleBZVectors.size()));
    k_prime_representatives_gammaP_minus_q = std::vector(irreducibleBZVectors.size(), std::vector<int>(irreducibleBZVectors.size()));

    k_prime_representatives_gammaM_plus_q = std::vector(irreducibleBZVectors.size(), std::vector<int>(irreducibleBZVectors.size()));
    k_prime_representatives_gammaZ_plus_q = std::vector(irreducibleBZVectors.size(), std::vector<int>(irreducibleBZVectors.size()));
    k_prime_representatives_gammaP_plus_q = std::vector(irreducibleBZVectors.size(), std::vector<int>(irreducibleBZVectors.size()));

#pragma omp parallel for
    for (int i = 0; i < irreducibleBZVectors.size(); i++)
    {
        for (int j = 0; j < irreducibleBZVectors.size(); j++)
        {
            k_prime_gammaM_minus_q[i][j] = -irreducibleBZVectors.at(j) + irreducibleBZVectors.at(i) + G_gammaM_minus_q[i][j]; // -k +q +  G
            k_prime_gammaZ_minus_q[i][j] = irreducibleBZVectors.at(j) + irreducibleBZVectors.at(i) + G_gammaZ_minus_q[i][j];  // k + q +  G
            k_prime_gammaP_minus_q[i][j] = -irreducibleBZVectors.at(j) - irreducibleBZVectors.at(i) - G_gammaP_minus_q[i][j]; // -k - q - G

            k_prime_gammaM_plus_q[i][j] = -irreducibleBZVectors.at(j) - irreducibleBZVectors.at(i) + G_gammaM_plus_q[i][j]; // -k -q +  G
            k_prime_gammaZ_plus_q[i][j] = irreducibleBZVectors.at(j) - irreducibleBZVectors.at(i) + G_gammaZ_plus_q[i][j];  // k - q +  G
            k_prime_gammaP_plus_q[i][j] = -irreducibleBZVectors.at(j) + irreducibleBZVectors.at(i) - G_gammaP_plus_q[i][j]; // -k + q - G

            k_prime_representatives_gammaM_minus_q[i][j] = findRepresentative(k_prime_gammaM_minus_q[i][j]);
            k_prime_representatives_gammaZ_minus_q[i][j] = findRepresentative(k_prime_gammaZ_minus_q[i][j]);
            k_prime_representatives_gammaP_minus_q[i][j] = findRepresentative(k_prime_gammaP_minus_q[i][j]);

            k_prime_representatives_gammaM_plus_q[i][j] = findRepresentative(k_prime_gammaM_plus_q[i][j]);
            k_prime_representatives_gammaZ_plus_q[i][j] = findRepresentative(k_prime_gammaZ_plus_q[i][j]);
            k_prime_representatives_gammaP_plus_q[i][j] = findRepresentative(k_prime_gammaP_plus_q[i][j]);
        }
    }

    */

    k_prime_representatives_gammaM_minus_q = std::vector(irreducibleBZVectors.size(), std::vector<int>(irreducibleBZVectors.size()));
    k_prime_representatives_gammaZ_minus_q = std::vector(irreducibleBZVectors.size(), std::vector<int>(irreducibleBZVectors.size()));
    k_prime_representatives_gammaP_minus_q = std::vector(irreducibleBZVectors.size(), std::vector<int>(irreducibleBZVectors.size()));

    k_prime_representatives_gammaM_plus_q = std::vector(irreducibleBZVectors.size(), std::vector<int>(irreducibleBZVectors.size()));
    k_prime_representatives_gammaZ_plus_q = std::vector(irreducibleBZVectors.size(), std::vector<int>(irreducibleBZVectors.size()));
    k_prime_representatives_gammaP_plus_q = std::vector(irreducibleBZVectors.size(), std::vector<int>(irreducibleBZVectors.size()));

    int sign_GammaM_plus_q[3] = {1, 1, 1};
    int sign_GammaP_plus_q[3] = {-1, -1, 1};
    int sign_GammaZ_plus_q[3] = {-1, 1, 1};
    int sign_GammaM_minus_q[3] = {1, 1, -1};
    int sign_GammaP_minus_q[3] = {-1, -1, -1};
    int sign_GammaZ_minus_q[3] = {-1, 1, -1};

#pragma omp parallel for
    for (int i = 0; i < irreducibleBZVectors.size(); i++)
    {
        for (int j = 0; j < irreducibleBZVectors.size(); j++)
        {
            auto k = irreducibleBZVectors.at(i);
            auto q = irreducibleBZVectors.at(j);

            k_prime_representatives_gammaM_minus_q[i][j] = get_k_prime(k, -q, sign_GammaM_minus_q);
            k_prime_representatives_gammaZ_minus_q[i][j] = get_k_prime(k, -q, sign_GammaZ_minus_q);
            k_prime_representatives_gammaP_minus_q[i][j] = get_k_prime(k, -q, sign_GammaP_minus_q);
            k_prime_representatives_gammaM_plus_q[i][j] = get_k_prime(k, q, sign_GammaM_plus_q);
            k_prime_representatives_gammaZ_plus_q[i][j] = get_k_prime(k, q, sign_GammaZ_plus_q);
            k_prime_representatives_gammaP_plus_q[i][j] = get_k_prime(k, q, sign_GammaP_plus_q);
        }
    }

    std::cout << "initialized k_prime vectors" << std::endl;
}

void IrreducibleBZ::save_k_primes_to_file(std::string filename)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Error: Could not open file " + filename);
    }
    for (int i = 0; i < irreducibleBZVectors.size(); i++)
    {
        for (int j = 0; j < irreducibleBZVectors.size(); j++)
        {
            file << i << " " << j << " " << k_prime_representatives_gammaM_minus_q[i][j] << " " << k_prime_representatives_gammaZ_minus_q[i][j] << " " << k_prime_representatives_gammaP_minus_q[i][j] << " " << k_prime_representatives_gammaM_plus_q[i][j] << " " << k_prime_representatives_gammaZ_plus_q[i][j] << " " << k_prime_representatives_gammaP_plus_q[i][j] << "\n";
        }
    }
    file.close();
}

void IrreducibleBZ::init_k_primes_from_file(std::string filename)
{
    k_prime_representatives_gammaM_minus_q = std::vector(irreducibleBZVectors.size(), std::vector<int>(irreducibleBZVectors.size()));
    k_prime_representatives_gammaZ_minus_q = std::vector(irreducibleBZVectors.size(), std::vector<int>(irreducibleBZVectors.size()));
    k_prime_representatives_gammaP_minus_q = std::vector(irreducibleBZVectors.size(), std::vector<int>(irreducibleBZVectors.size()));
    k_prime_representatives_gammaP_plus_q = std::vector(irreducibleBZVectors.size(), std::vector<int>(irreducibleBZVectors.size()));
    k_prime_representatives_gammaZ_plus_q = std::vector(irreducibleBZVectors.size(), std::vector<int>(irreducibleBZVectors.size()));
    k_prime_representatives_gammaM_plus_q = std::vector(irreducibleBZVectors.size(), std::vector<int>(irreducibleBZVectors.size()));

    std::ifstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Error: Could not open file " + filename);
    }

    while (!file.eof())
    {
        int i, j, k_prime_gammaM_minus_q, k_prime_gammaZ_minus_q, k_prime_gammaP_minus_q, k_prime_gammaM_plus_q, k_prime_gammaZ_plus_q, k_prime_gammaP_plus_q;
        file >> i >> j >> k_prime_gammaM_minus_q >> k_prime_gammaZ_minus_q >> k_prime_gammaP_minus_q >> k_prime_gammaM_plus_q >> k_prime_gammaZ_plus_q >> k_prime_gammaP_plus_q;
        k_prime_representatives_gammaM_minus_q[i][j] = k_prime_gammaM_minus_q;
        k_prime_representatives_gammaZ_minus_q[i][j] = k_prime_gammaZ_minus_q;
        k_prime_representatives_gammaP_minus_q[i][j] = k_prime_gammaP_minus_q;
        k_prime_representatives_gammaM_plus_q[i][j] = k_prime_gammaM_plus_q;
        k_prime_representatives_gammaZ_plus_q[i][j] = k_prime_gammaZ_plus_q;
        k_prime_representatives_gammaP_plus_q[i][j] = k_prime_gammaP_plus_q;
    }
}

// given a line with a symmetry operatation e.g. "+x,-y,z" return a function that applies this operation to a vector
std::function<Eigen::Vector3d(Eigen::Vector3d)> parseOperation(const std::string &line)
{
    return [line](Eigen::Vector3d vec) -> Eigen::Vector3d
    {
        Eigen::Vector3d transformed = Eigen::Vector3d::Zero();
        std::istringstream iss(line);
        for (int i = 0; i < 3; ++i)
        {
            char signChar = '+', axisChar, commaChar;
            iss >> signChar; // Read the possible sign character
            if (signChar != '-' && signChar != '+')
            {
                // If not a sign, it's an axis character; adjust accordingly
                axisChar = signChar;
                signChar = '+'; // Default to positive if no sign is explicitly given
            }
            else
            {
                // If a sign was read, the next character is the axis character
                iss >> axisChar;
            }

            // After reading sign and axis, attempt to read a comma
            if (i < 2)
            {                     // No comma after the last axis character
                iss >> commaChar; // Read and discard the comma
                if (commaChar != ',')
                {
                    // If the expected comma is not found, something went wrong
                    std::cerr << "Error: Expected comma in operation string." << std::endl;
                    return vec; // Return the original vector unchanged
                }
            }

            double sign = (signChar == '-' ? -1.0 : 1.0);
            switch (axisChar)
            {
            case 'x':
                transformed[i] = sign * vec[0];
                break;
            case 'y':
                transformed[i] = sign * vec[1];
                break;
            case 'z':
                transformed[i] = sign * vec[2];
                break;
            default:
                std::cerr << "Error: Unexpected character in operation string." << std::endl;
                return vec; // Return the original vector unchanged
            }
        }
        return transformed;
    };
}

void IrreducibleBZ::initSymmOp(std::string filepath)
{
    std::vector<std::function<Eigen::Vector3d(Eigen::Vector3d)>> operations;

    std::string line;
    std::ifstream file(filepath);

    // Check if file opened successfully
    if (!file.is_open())
    {
        std::cerr << "Could not open the file." << std::endl;
        return;
    }

    // Read the operations from the file
    while (std::getline(file, line))
    {
        operations.push_back(parseOperation(line));
    }

    this->operations = operations;

    file.close();
}

std::vector<Eigen::Vector3d> IrreducibleBZ::applySymmOp(const Eigen::Vector3d &vec)
{

    std::vector<Eigen::Vector3d> transformed_vectors;
    transformed_vectors.reserve(operations.size());

    for (auto &op : operations)
    {
        Eigen::Vector3d transformedVec = op(vec);
        transformed_vectors.push_back(transformedVec);
    }
    return transformed_vectors;
}
