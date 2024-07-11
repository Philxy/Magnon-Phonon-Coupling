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
        // std::cout << phDisp.polVectors[0][branch] * phDisp.polVectors[0][branch] + phDisp.polVectors[1][branch] * phDisp.polVectors[1][branch] + phDisp.polVectors[2][branch] * phDisp.polVectors[2][branch] << std::endl;
        // std::cout << phDisp.polVectors[0][branch] << " " << phDisp.polVectors[1][branch] << " " << phDisp.polVectors[2][branch] << std::endl;
        // std::cout << "____\n";

        for (int axis = 0; axis < 3; axis++)
        {
            std::complex<double> D_plus = D_k_values_minus.D[axis][0] + i * D_k_values_minus.D[axis][1];
            std::complex<double> D_minus = D_k_values.D[axis][0] - i * D_k_values.D[axis][1];

            double pol_vec_component = phDisp.polVectors[axis][branch];

            C.at(branch) += 2.0 * i / sqrt(2 * S) * 3.8636 * D_minus * pol_vec_component * sqrt(1 / (2 * atomicMass * phDisp.E[branch]));
            D.at(branch) += -2.0 * i / sqrt(2 * S) * 3.8636 * D_plus * pol_vec_component * sqrt(1 / (2 * atomicMass * phDisp.E[branch]));
        }
    }
    return;
}

void Diagonalization::calcAngularMomentumFromEigenvectors()
{

    std::vector<Eigen::Vector3cd> complex_pol_vectors;

    for (int row = 0; row < 8; row++)
    {
        Eigen::Vector3cd pol_vector;

        for (int col = 4; col < 7; col++)
        {
            pol_vector(col - 4) = eigenvectors_inverse(row, col);
        }
        complex_pol_vectors.push_back(pol_vector);

        for (int col = 0; col < 3; col++)
        {
            pol_vector(col) = eigenvectors_inverse(row, col);
        }
        // complex_pol_vectors.push_back(pol_vector);
    }

    for (Eigen::Vector3cd pol_vec : complex_pol_vectors)
    {
        Eigen::Vector3d real_part = pol_vec.real();
        Eigen::Vector3d imag_part = pol_vec.imag();
        Eigen::Vector3d angul_momentum = 2 * real_part.cross(imag_part);
        angularMomentumFromEigenvectors.push_back(angul_momentum);
    }

    return;
}

void Diagonalization::calcDMILike()
{
    // Dxx, Dxy, Dxz, Dyx, Dyy, Dyz
    DMILike dmiLike_plus_k(k.x, k.y, k.z, couplingParameters);

    for (int mu = 0; mu < 3; mu++)
    {
        dmiLike.push_back(dmiLike_plus_k.D[X][mu]);
    }

    for (int mu = 0; mu < 3; mu++)
    {
        dmiLike.push_back(dmiLike_plus_k.D[Y][mu]);
    }
}

// Calculates the DMI-like coefficients projected onto the polarization vectors, termed C and D or D_minus and D_plus
void Diagonalization::calculateCD(bool dmiOnly)
{

    const std::complex<double> i(0, 1);

    if (dmiOnly)
    {

        DMILike dmiLike_plus_k(k.x, k.y, k.z, couplingParameters);
        DMILike dmiLike_minus_k(-k.x, -k.y, -k.z, couplingParameters);

        std::complex<double> DPlus_minus_k[3];
        std::complex<double> DMinus_plus_k[3];

        for (int mu = 0; mu < 3; mu++)
        {
            DPlus_minus_k[mu] = dmiLike_minus_k.D[X][mu] + i * dmiLike_minus_k.D[Y][mu];
            DMinus_plus_k[mu] = dmiLike_plus_k.D[X][mu] - i * dmiLike_plus_k.D[Y][mu];
        }

        for (int branch = 0; branch < 3; branch++)
        {
            for (int axis = 0; axis < 3; axis++)
            {
                double pol_vec_component = phDisp.polVectors[axis][branch];
                std::complex<double> C_add = 2.0 * i / sqrt(2 * S) * 3.8636 * DMinus_plus_k[axis] * pol_vec_component * sqrt(1 / (2 * atomicMass * phDisp.E[branch]));
                std::complex<double> D_add = -2.0 * i / sqrt(2 * S) * 3.8636 * DPlus_minus_k[axis] * pol_vec_component * sqrt(1 / (2 * atomicMass * phDisp.E[branch]));
                C.at(branch) += C_add;
                D.at(branch) += D_add;
            }
        }
        return;
    }
    else
    {
        // Includes the contribution of the anisotropy and the DM

        std::complex<double> J_xz_plus[3] = {0, 0, 0};
        std::complex<double> J_yz_plus[3] = {0, 0, 0};
        std::complex<double> J_xz_minus[3] = {0, 0, 0};
        std::complex<double> J_yz_minus[3] = {0, 0, 0};

        for (Axis ax : {X, Y, Z})
        {
            J_xz_plus[ax] = J_kq(0, 0, 0, k.x, k.y, k.z, couplingParameters, X, Z, ax);
            J_yz_plus[ax] = J_kq(0, 0, 0, k.x, k.y, k.z, couplingParameters, Y, Z, ax);
            J_xz_minus[ax] = J_kq(0, 0, 0, k.x, k.y, k.z, couplingParameters, X, Z, ax);
            J_yz_minus[ax] = J_kq(0, 0, 0, k.x, k.y, k.z, couplingParameters, Y, Z, ax);
        }

        for (int branch = 0; branch < 3; branch++)
        {
            for (int axis = 0; axis < 3; axis++)
            {
                double pol_vec_component = phDisp.polVectors[axis][branch];
                std::complex<double> C_add = 2.0 / sqrt(2 * S) * 3.8636 * pol_vec_component * sqrt(1 / (2 * atomicMass * phDisp.E[branch])) * (J_xz_plus[axis] - i * J_yz_plus[axis]);
                std::complex<double> D_add = 2.0 / sqrt(2 * S) * 3.8636 * pol_vec_component * sqrt(1 / (2 * atomicMass * phDisp.E[branch])) * (J_xz_minus[axis] + i * J_yz_minus[axis]);
                C.at(branch) += C_add;
                D.at(branch) += D_add;
            }
        }
        return;
    }
}

void Diagonalization::calcAB()
{
    std::complex<double> i(0, 1);
    std::complex<double> A(0, 0);
    std::complex<double> B(0, 0);
    Eigen::Vector3cd chi_k1(phDisp.polVectors[0][0], phDisp.polVectors[1][0], phDisp.polVectors[2][0]);
    Eigen::Vector3cd chi_k2(phDisp.polVectors[0][1], phDisp.polVectors[1][1], phDisp.polVectors[2][1]);
    Eigen::Vector3cd chiP = chi_k1 + i * chi_k2;
    Eigen::Vector3cd chiM = chi_k1 - i * chi_k2;

    DMILike dmiLike(k.x, k.y, k.z, couplingParameters);

    Eigen::Vector3cd DM(dmiLike.D[0][0] - i * dmiLike.D[1][0], dmiLike.D[0][1] - i * dmiLike.D[1][1], dmiLike.D[0][2] - i * dmiLike.D[1][2]);

    A = i * 3.84 * sqrt(1 / (S * atomicMass * phDisp.E[0])) * (DM(0) * chiP(0) + DM(1) * chiP(1) + DM(2) * chiP(2));
    B = i * 3.84 * sqrt(1 / (S * atomicMass * phDisp.E[1])) * (DM(0) * chiM(0) + DM(1) * chiM(1) + DM(2) * chiM(2));
    this->A = A;
    this->B = B;

    // Alternatively (simplification along G-H path):
    // this->A = dmiLike.D[1][0] - dmiLike.D[0][1] + i * (dmiLike.D[0][0] + dmiLike.D[1][1]);
    // this->B = dmiLike.D[1][0] + dmiLike.D[0][1] + i * (dmiLike.D[0][0] - dmiLike.D[1][1]);
}

// Calculates the matrix representation of the Hamiltonian and prepares it for diagonalization using the generalized Bogoliubov transformation as outlined by White et al. (2014)
void Diagonalization::calcMatrixHamiltonian()
{

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

    // Magnon energy
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
  and D a diagonal matrix with the eigenvalues on the diagonal
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

void Diagonalization::calcAngularMomentum()
{

    Eigen::MatrixXcd matrixL = Eigen::MatrixXcd::Zero(8, 8);
    Eigen::MatrixXcd g = getMatrix_g();
    std::complex<double> i(0, 1);

    double polVec = phDisp.polVectors[2][2]; // phDisp.polVectors[1][2]; //

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

    matrixL(0, 5) = -polVec;
    matrixL(1, 4) = polVec;

    matrixL(5, 4) = -polVec;
    matrixL(4, 5) = polVec;

    auto D = i / 2.0 * eigenvectors_inverse * matrixL * eigenvectors;

    // std::cout << phDisp.polVectors[0][0] << phDisp.polVectors[1][0] << phDisp.polVectors[2][0] << std::endl;
    // std::cout << phDisp.polVectors[0][1] << phDisp.polVectors[1][1] << phDisp.polVectors[2][1] << std::endl;
    // std::cout << phDisp.polVectors[0][2] << phDisp.polVectors[1][2] << phDisp.polVectors[2][2] << std::endl;
    // std::cout << D << std::endl;
    // std::cout << polVec << std::endl;

    for (int j = 0; j < 8; j++)
    {
        angularMomentum.push_back(D(j, j));
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
