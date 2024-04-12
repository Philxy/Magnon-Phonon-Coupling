#include "../include/diagonalization.h"

Diagonalization::Diagonalization(const std::vector<CouplingParameter> &couplingParameters, const PhononDispParam &phDispParam, const MagnonDispParam &magDispParam, const Vector3D &kVec)
    : couplingParameters(couplingParameters), phDisp(phDispParam), magDisp(magDispParam), k(kVec)
{
    this->C = std::vector<std::complex<double>>(3, std::complex<double>(0, 0));
    this->D = std::vector<std::complex<double>>(3, std::complex<double>(0, 0));
    matrixHamiltonian = Eigen::MatrixXcd(8, 8);
}

void Diagonalization::calcCD()
{
    DMILikeCouplingParam D_k_values(k.x, k.y, k.z, couplingParameters);
    DMILikeCouplingParam D_k_values_minus(-k.x, -k.y, -k.z, couplingParameters);
    const std::complex<double> i(0, 1);

    for (int branch = 0; branch < 3; branch++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            std::complex<double> D_plus = D_k_values_minus.D[0][axis] + i * D_k_values_minus.D[1][axis];
            std::complex<double> D_minus = D_k_values.D[0][axis] - i * D_k_values.D[1][axis];

            double pol_vec_component = phDisp.polVectors[axis][branch];

            C.at(branch) += 2.0 * i / sqrt(2 * S) * 3.8636 * D_minus * pol_vec_component * sqrt(1 / (2 * atomicMass * phDisp.E[branch]));
            D.at(branch) += -2.0 * i / sqrt(2 * S) * 3.8636 * D_plus * pol_vec_component * sqrt(1 / (2 * atomicMass * phDisp.E[branch]));
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
/*
Performs the diagonalization of gH
- If V is a matrix with the eigenvectors as its columns
  and D a diagonal matrix withe the eigenvalues on the diagonal
  then gH V = V D -> V^-1 gH V = D
*/
void Diagonalization::diagonalize()
{
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver;
    Eigen::MatrixXcd matrix = getMatrix_g() * matrixHamiltonian;
    solver.compute(matrix);

    eigenvectors = solver.eigenvectors();
    eigenvectors_inverse = solver.eigenvectors().inverse();

    // Loop over each eigenvalue
    for (int i = 0; i < solver.eigenvalues().size(); ++i)
    {
        double realPart = solver.eigenvalues()[i];

        eigenEnergies.push_back(realPart);
    }
}

int getLeviLevita(int i, int j, int k)
{
    return 0;
}

void Diagonalization::calcAngularMomentum(double time)
{

    Eigen::MatrixXcd matrixL = Eigen::MatrixXcd::Zero(8, 8);
    Eigen::MatrixXcd g = getMatrix_g();
    std::complex<double> i(0, 1);

    double polVec = phDisp.polVectors[0][0] * phDisp.polVectors[2][1] - phDisp.polVectors[2][0] * phDisp.polVectors[0][1]; // phDisp.polVectors[1][2]; //

    assert(matrixL.rows() == 8);
    assert(matrixL.cols() == 8);
    assert(g.rows() == 8);
    assert(g.cols() == 8);
    assert(eigenvectors.rows() == 8);
    assert(eigenvectors.cols() == 8);
    assert(eigenvectors_inverse.rows() == 8);
    assert(eigenvectors_inverse.cols() == 8);

    matrixL(0, 1) = -polVec;
    matrixL(1, 0) = polVec;

    matrixL(4, 1) = polVec;
    matrixL(5, 0) = -polVec;

    matrixL(1, 4) = polVec;
    matrixL(5, 4) = -polVec;

    matrixL(0, 5) = -polVec;
    matrixL(4, 5) = polVec;


    auto D = i / 2.0 * eigenvectors_inverse * matrixL * eigenvectors;

    // std::cout << phDisp.polVectors[0][0] << phDisp.polVectors[1][0] << phDisp.polVectors[2][0] << std::endl;
    // std::cout << phDisp.polVectors[0][1] << phDisp.polVectors[1][1] << phDisp.polVectors[2][1] << std::endl;
    // std::cout << phDisp.polVectors[0][2] << phDisp.polVectors[1][2] << phDisp.polVectors[2][2] << std::endl;
    // std::cout << D << std::endl;
    // std::cout << polVec << std::endl;

    for (int i = 0; i < 8; i++)
    {
        angularMomentum.push_back(D(i,i));
    }

    return;
    /*
    Eigen::Vector3d L;
    std::complex<double> i(0, 1);

    for (int nu = 0; nu < 8; nu++)
    {
        Eigen::Vector3cd sum = Eigen::Vector3cd::Zero();
        Eigen::VectorXcd alpha_nu = eigenvectors_inverse.row(nu);
        Eigen::VectorXcd alpha_dagger_nu = alpha_nu.conjugate();

        for (int lambda = 0; lambda < 3; lambda++)
        {
            for (int lambda_prime = 0; lambda_prime < 3; lambda_prime++)
            {
                double E_lam = phDisp.E[lambda];
                double E_lam_prime = phDisp.E[lambda_prime];
                Eigen::Vector3d polVec_lam(phDisp.polVectors[0][lambda], phDisp.polVectors[1][lambda], phDisp.polVectors[2][lambda]);
                Eigen::Vector3d polVec_lam_prime(phDisp.polVectors[0][lambda_prime], phDisp.polVectors[1][lambda_prime], phDisp.polVectors[2][lambda_prime]);
                Eigen::Vector3d polVecCrossProduct = polVec_lam_prime.cross(polVec_lam);
                std::complex<double> expFunc = std::exp(i * (E_lam - E_lam_prime) * time);
                sum += i * sqrt(E_lam_prime / E_lam) * alpha_nu.dot(alpha_dagger_nu) * polVecCrossProduct * expFunc;
            }
        }
        double Lx = sqrt(sum(0).real() * sum(0).real() + sum(0).imag() * sum(0).imag());
        double Ly = sqrt(sum(1).real() * sum(1).real() + sum(1).imag() * sum(1).imag());
        double Lz = sqrt(sum(2).real() * sum(2).real() + sum(2).imag() * sum(2).imag());

        // L = Eigen::Vector3d(Lx, Ly, Lz);
        L = sum.real();
        angularMomentum.push_back(L);
    }
    */
}

std::vector<std::vector<double>> diagonalizeHamiltonian(const std::vector<Vector3D> &path, const std::vector<PhononDispParam> &phononDispersion, const std::vector<MagnonDispParam> &magnonDispersion, const std::vector<CouplingParameter> &parameters)
{

    std::complex<double> i(0, 1);

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
                C.at(branch) += 2.0 * i / sqrt(2 * S) * 3.8636 * D_minus * pol_vec_component * sqrt(1 / (2 * atomicMass * E_k_branch));
                D.at(branch) += -2.0 * i / sqrt(2 * S) * 3.8636 * D_plus * pol_vec_component * sqrt(1 / (2 * atomicMass * E_k_branch));
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
