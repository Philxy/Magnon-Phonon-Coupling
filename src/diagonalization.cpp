#include "../include/diagonalization.h"

Diagonalization::Diagonalization(const std::vector<CouplingParameter> &couplingParameters, PhononDispParam &phDispParam, const MagnonDispParam &magDispParam, const Vector3D &kVec, double atomicMass, double S)
    : couplingParameters(couplingParameters), phDisp(phDispParam), magDisp(magDispParam), k(kVec), atomic_mass(atomicMass), S(S)
{
    this->C = std::vector<std::complex<double>>(3, std::complex<double>(0, 0));
    this->D = std::vector<std::complex<double>>(3, std::complex<double>(0, 0));
    matrixHamiltonian = Eigen::MatrixXcd(8, 8);
}

void Diagonalization::calcCD()
{
    DMILikeCouplingParam D_k_values(k.x, k.y, k.z, couplingParameters);
    DMILikeCouplingParam D_k_values_minus(-k.x, -k.y, -k.z, couplingParameters);
    std::complex<double> i(0, 1);

    for (int branch = 0; branch < 3; branch++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            std::complex<double> D_plus = D_k_values_minus.D[0][axis] + i * D_k_values_minus.D[1][axis];
            std::complex<double> D_minus = D_k_values.D[0][axis] - i * D_k_values.D[1][axis];
            double pol_vec_component = phDisp.polVectors[axis][branch];
            C.at(branch) += 2.0 * i / sqrt(2 * S) * 3.8636 * D_minus * pol_vec_component * sqrt(1 / (2 * atomic_mass * phDisp.E[branch]));
            D.at(branch) += -2.0 * i / sqrt(2 * S) * 3.8636 * D_plus * pol_vec_component * sqrt(1 / (2 * atomic_mass * phDisp.E[branch]));
        }
    }
    return;
}

void Diagonalization::calcMatrixHamiltonian()
{
    Eigen::MatrixXd g = getMatrix_g();

    // set all elements to zero
    for (int row = 0; row < 8; row++)
    {
        for (int col = 0; col < 8; col++)
        {
            matrixHamiltonian(row, col) = std::complex<double>(0, 0);
        }
    }

    for (int branch = 0; branch < 3; branch++)
    {

        /*
        if (branch == 3)
        {
            C.at(branch) = 0;
            D.at(branch) = 0;
        }
        */

        // Diagonal phonon elements
        double E_k_branch = phDisp.E[branch];
        matrixHamiltonian(branch, branch) = E_k_branch;
        matrixHamiltonian(branch + 4, branch + 4) = E_k_branch;

        // off-diagonal components
        // upper right quadrant
        matrixHamiltonian(branch, 7) = D.at(branch);
        matrixHamiltonian(3, 4 + branch) = D.at(branch);
        // lower left quadrant
        matrixHamiltonian(4 + branch, 3) = C.at(branch);
        matrixHamiltonian(7, branch) = C.at(branch);
        // upper left quadrant
        matrixHamiltonian(branch, 3) = C.at(branch);
        matrixHamiltonian(3, branch) = D.at(branch);
        // lower right quadrant
        matrixHamiltonian(4 + branch, 7) = D.at(branch);
        matrixHamiltonian(7, 4 + branch) = C.at(branch);

        /*
        // isolate a mode
        int c = 2;
        if (branch != c)
        {
            matrixHamiltonian(branch, 7) = 0;     // D.at(branch);
            matrixHamiltonian(3, 4 + branch) = 0; // D.at(branch);
            // lower left quadrant
            matrixHamiltonian(4 + branch, 3) = 0; // C.at(branch);
            matrixHamiltonian(7, branch) = 0;     // C.at(branch);
            // upper left quadrant
            matrixHamiltonian(branch, 3) = 0; // C.at(branch);
            matrixHamiltonian(3, branch) = 0; // D.at(branch);
            // lower right quadrant
            matrixHamiltonian(4 + branch, 7) = 0; // D.at(branch);
            matrixHamiltonian(7, 4 + branch) = 0; // C.at(branch);

            matrixHamiltonian(branch, branch) = 0;
            matrixHamiltonian(branch + 4, branch + 4) = 0;
        }
        */
    }

    double E_k_mag = magDisp.energy;
    matrixHamiltonian(3, 3) = E_k_mag;
    matrixHamiltonian(7, 7) = E_k_mag;

    return;
}

Eigen::MatrixXd getMatrix_g()
{
    Eigen::MatrixXd g(8, 8);

    for (int row = 0; row < 8; row++)
    {
        for (int col = 0; col < 8; col++)
        {
            if (row != col)
            {
                g(row, col) = 0;
            }
            if (row == col && row < 4)
            {
                g(row, col) = 1;
            }
            if (row == col && row >= 4)
            {
                g(row, col) = -1;
            }
        }
    }
    return g;
}

void Diagonalization::diagonalize()
{
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver;
    Eigen::MatrixXcd matrix = getMatrix_g() * matrixHamiltonian;
    solver.compute(matrix, true);

    eigenvectors = solver.eigenvectors();
    eigenvectors_inverse = solver.eigenvectors().inverse();

    // Loop over each eigenvalue
    for (int i = 0; i < solver.eigenvalues().size(); ++i)
    {
        double realPart = solver.eigenvalues()[i].real();

        eigenEnergies.push_back(realPart);
    }
}

std::vector<std::vector<double>> diagonalizeHamiltonian(const std::vector<Vector3D> &path, const std::vector<PhononDispParam> &phononDispersion, const std::vector<MagnonDispParam> &magnonDispersion, const std::vector<CouplingParameter> &parameters)
{

    std::complex<double> i(0, 1);
    double S = 1.115;
    double atomic_mass = 55.845;

    std::vector<std::vector<double>> allEigenvectors;

    for (int idx = 0; idx < path.size(); idx++)
    {
        Vector3D k = path.at(idx);
        std::vector<double> eigenenergies; // eigenenergies for current k vector

        std::vector<std::complex<double>> C(3, std::complex<double>(0, 0));
        std::vector<std::complex<double>> D(3, std::complex<double>(0, 0));

        for (int branch = 0; branch < 3; branch++)
        {

            double E_k_branch = phononDispersion.at(idx).E[branch];

            DMILikeCouplingParam D_k_values(k.x, k.y, k.z, parameters);
            DMILikeCouplingParam D_k_values_minus(-k.x, -k.y, -k.z, parameters);

            for (int axis = 0; axis < 3; axis++)
            {
                std::complex<double> D_plus = D_k_values_minus.D[0][axis] + i * D_k_values_minus.D[1][axis];
                std::complex<double> D_minus = D_k_values.D[0][axis] - i * D_k_values.D[1][axis];
                double pol_vec_component = phononDispersion.at(idx).polVectors[axis][branch];
                C.at(branch) += 2.0 * i / sqrt(2 * S) * 3.8636 * D_minus * pol_vec_component * sqrt(1 / (2 * atomic_mass * E_k_branch));
                D.at(branch) += -2.0 * i / sqrt(2 * S) * 3.8636 * D_plus * pol_vec_component * sqrt(1 / (2 * atomic_mass * E_k_branch));
            }
        }

        Eigen::MatrixXcd matrixH(8, 8);
        Eigen::MatrixXd g = getMatrix_g();

        // set all elements to zero
        for (int row = 0; row < 8; row++)
        {
            for (int col = 0; col < 8; col++)
            {
                matrixH(row, col) = std::complex<double>(0, 0);
            }
        }

        for (int branch = 0; branch < 3; branch++)
        {
            // Diagonal phonon elements
            double E_k_branch = phononDispersion.at(idx).E[branch];
            matrixH(branch, branch) = E_k_branch;
            matrixH(branch + 4, branch + 4) = E_k_branch;

            // off-diagonal components

            // upper right quadrant
            matrixH(branch, 7) = D.at(branch);     // D.at(branch);
            matrixH(3, 4 + branch) = D.at(branch); // D.at(branch);
            // lower left quadrant
            matrixH(4 + branch, 3) = C.at(branch); // C.at(branch);
            matrixH(7, branch) = C.at(branch);     // C.at(branch);
            // upper left quadrant
            matrixH(branch, 3) = C.at(branch); // C.at(branch);
            matrixH(3, branch) = D.at(branch); // D.at(branch);
            // lower right quadrant
            matrixH(4 + branch, 7) = D.at(branch); // D.at(branch);
            matrixH(7, 4 + branch) = C.at(branch); // C.at(branch);
        }

        double E_k_mag = magnonDispersion.at(idx).energy;
        matrixH(3, 3) = E_k_mag;
        matrixH(7, 7) = E_k_mag;

        // std::cout << matrixH << std::endl;
        // std::cout << " " << std::endl;

        Eigen::MatrixXcd matrix(8, 8);
        matrix = g * matrixH;

        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver;
        solver.compute(matrix);
        // Loop over each eigenvalue
        for (int i = 0; i < solver.eigenvalues().size(); ++i)
        {
            // Access the real part of the ith eigenvalue
            double realPart = solver.eigenvalues()[i].real();

            // Output the real part of the eigenvalue
            std::cout << realPart << std::endl;

            // Add the real part to your std::vector
            eigenenergies.push_back(realPart);
        }
        allEigenvectors.push_back(eigenenergies);
    }

    return allEigenvectors;
}
