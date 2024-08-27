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

void Diagonalization::calcPolVectors()
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
        polVectors.push_back(pol_vec);
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
void Diagonalization::calculateCD(Interaction interaction)
{

    const std::complex<double> i(0, 1);

    switch (interaction)
    {
    case DMI:
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

        assert(std::abs(std::conj(C.at(0)) - D.at(0)) < 1e-7 && "coefficients are not conjugate");
        assert(std::abs(std::conj(C.at(1)) - D.at(1)) < 1e-7 && "coefficients are not conjugate");
        assert(std::abs(std::conj(C.at(2)) - D.at(2)) < 1e-7 && "coefficients are not conjugate");

        break;
    }
    case ALL:
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
            J_xz_minus[ax] = J_kq(0, 0, 0, -k.x, -k.y, -k.z, couplingParameters, X, Z, ax);
            J_yz_minus[ax] = J_kq(0, 0, 0, -k.x, -k.y, -k.z, couplingParameters, Y, Z, ax);
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
        break;
    }
    case ANISOTROPY:
    {
        // First, calculate the contribution of all interactions

        std::complex<double> J_xz_plus[3] = {0, 0, 0};
        std::complex<double> J_yz_plus[3] = {0, 0, 0};
        std::complex<double> J_xz_minus[3] = {0, 0, 0};
        std::complex<double> J_yz_minus[3] = {0, 0, 0};

        // Calculate the contribution of the DMI and subtract it
        DMILike dmiLike_plus_k(k.x, k.y, k.z, couplingParameters);
        DMILike dmiLike_minus_k(-k.x, -k.y, -k.z, couplingParameters);

        std::complex<double> DPlus_minus_k[3];
        std::complex<double> DMinus_plus_k[3];

        for (int mu = 0; mu < 3; mu++)
        {
            DPlus_minus_k[mu] = dmiLike_minus_k.D[X][mu] + i * dmiLike_minus_k.D[Y][mu];
            DMinus_plus_k[mu] = dmiLike_plus_k.D[X][mu] - i * dmiLike_plus_k.D[Y][mu];

            J_xz_plus[mu] = J_kq(0, 0, 0, k.x, k.y, k.z, couplingParameters, X, Z, static_cast<Axis>(mu));
            J_yz_plus[mu] = J_kq(0, 0, 0, k.x, k.y, k.z, couplingParameters, Y, Z, static_cast<Axis>(mu));
            J_xz_minus[mu] = J_kq(0, 0, 0, -k.x, -k.y, -k.z, couplingParameters, X, Z, static_cast<Axis>(mu));
            J_yz_minus[mu] = J_kq(0, 0, 0, -k.x, -k.y, -k.z, couplingParameters, Y, Z, static_cast<Axis>(mu));
        }

        for (int branch = 0; branch < 3; branch++)
        {
            for (int axis = 0; axis < 3; axis++)
            {
                double pol_vec_component = phDisp.polVectors[axis][branch];

                std::complex<double> C_add = 2.0 / sqrt(2 * S) * 3.8636 * pol_vec_component * sqrt(1 / (2 * atomicMass * phDisp.E[branch])) * (J_xz_plus[axis] - i * J_yz_plus[axis]);
                std::complex<double> D_add = 2.0 / sqrt(2 * S) * 3.8636 * pol_vec_component * sqrt(1 / (2 * atomicMass * phDisp.E[branch])) * (J_xz_minus[axis] + i * J_yz_minus[axis]);

                std::complex<double> C_subtract = 2.0 * i / sqrt(2 * S) * 3.8636 * DMinus_plus_k[axis] * pol_vec_component * sqrt(1 / (2 * atomicMass * phDisp.E[branch]));
                std::complex<double> D_subtract = -2.0 * i / sqrt(2 * S) * 3.8636 * DPlus_minus_k[axis] * pol_vec_component * sqrt(1 / (2 * atomicMass * phDisp.E[branch]));

                C.at(branch) += C_add;
                D.at(branch) += D_add;

                C.at(branch) += C_subtract;
                D.at(branch) += D_subtract;
            }
        }
        break;
    }
    default:
        break;
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

bool isDiagonal(const Eigen::MatrixXcd &matrix)
{
    double epsilon = 1e-7;
    for (int row = 0; row < matrix.rows(); row++)
    {
        for (int col = 0; col < matrix.cols(); col++)
        {
            if (row != col && std::abs(matrix(row, col)) >= epsilon)
            {
                return false;
            }
        }
    }
    return true;
}

double getLargestNonDiagonalElement(const Eigen::MatrixXcd &matrix)
{
    double largestNonDiagonalElement = 0.0;

    for (int row = 0; row < matrix.rows(); row++)
    {
        for (int col = 0; col < matrix.cols(); col++)
        {
            if (row != col)
            {
                double element = std::abs(matrix(row, col));
                if (element > largestNonDiagonalElement)
                {
                    largestNonDiagonalElement = element;
                }
            }
        }
    }

    return largestNonDiagonalElement;
}

void printMatrix(const Eigen::MatrixXcd &matrix)
{
    int width = 15;    // width for each element
    int precision = 3; // number of decimal places

    for (int row = 0; row < matrix.rows(); row++)
    {
        for (int col = 0; col < matrix.cols(); col++)
        {
            std::cout << std::setw(width)
                      << std::setprecision(precision)
                      << std::fixed
                      << matrix(row, col)
                      << " ";
        }
        std::cout << std::endl;
    }
}

/*
Performs the diagonalization of gH
- If V is a matrix with the eigenvectors as its columns
  and D a diagonal matrix with the eigenvalues on the diagonal
  then gH V = V D -> V^-1 gH V = D
*/
void Diagonalization::diagonalize()
{
    Eigen::MatrixXcd matrix = getMatrix_g() * matrixHamiltonian;
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(matrix);

    eigenvectors = solver.eigenvectors();
    eigenvectors_inverse = solver.eigenvectors().inverse();

    // Loop over each eigenvalue
    for (int i = 0; i < solver.eigenvalues().size(); ++i)
    {
        double realPart = solver.eigenvalues()[i].real();

        eigenEnergies.push_back(realPart);
    }

    // std::cout << "matrix Hamiltonian: H " << std::endl;
    // printMatrix(matrixHamiltonian);

    // std::cout << "supposedly diagonal matrix: Q^dagger H Q " << std::endl;
    // printMatrix((eigenvectors.transpose()).conjugate() * matrixHamiltonian * eigenvectors);
    auto diagMat = (eigenvectors.transpose()).conjugate() * matrixHamiltonian * eigenvectors;
    if (!isDiagonal(diagMat))
    {
        // printMatrix(diagMat);
        std::cout << "k-vector: " << k.x << " " << k.y << " " << k.z << std::endl;
        // set the print precision to 10 decimal places
        std::cout << std::setw(15)
                  << std::setprecision(10)
                  << std::fixed;
        std::cout << "larges non diag element: " << getLargestNonDiagonalElement(diagMat) << std::endl;
        std::cout << "_________________" << std::endl;
    }
    // assert(isDiagonal((eigenvectors.transpose()).conjugate() * matrixHamiltonian * eigenvectors) && "matrix is not diagonal");

    // std::cout << "numerically diagonalized matrix: Q^-1 g H Q" << std::endl;
    // printMatrix(eigenvectors_inverse * matrix * eigenvectors);

    // std::cout << "check if Hamiltonian is self adjoint: H-H^dagger == 0 " << std::endl;
    // printMatrix(matrixHamiltonian.adjoint() - matrixHamiltonian);

    // std::cout << "calculate g^prime:" << std::endl;
    Eigen::MatrixXcd g_prime = eigenvectors_inverse * getMatrix_g() * eigenvectors.adjoint().inverse();
    // printMatrix(g_prime);
    // assert(isDiagonal(g_prime) && "g_prime is not diagonal");

    // std::cout << "__________" << std::endl;

    return;
}

// Returns the square root of a diagonal matrix (element-wise)
Eigen::MatrixXcd squaredMatrix(const Eigen::MatrixXcd &matrix)
{
    Eigen::MatrixXcd matSqrt = Eigen::MatrixXcd::Zero(matrix.rows(), matrix.cols());
    for (int i = 0; i < matrix.rows(); ++i)
    {
        matSqrt(i, i) = std::sqrt(matrix(i, i)); // compute the square root of each diagonal element
    }
    return matSqrt;
}

bool equalityToIdentityMatrix(const Eigen::MatrixXcd &matrix)
{
    double epsilon = 1e-10;
    Eigen::MatrixXcd identity = Eigen::MatrixXcd::Identity(matrix.rows(), matrix.cols());
    for (int i = 0; i < matrix.rows(); ++i)
    {
        for (int j = 0; j < matrix.cols(); ++j)
        {
            if (std::abs(matrix(i, j) - identity(i, j)) > epsilon)
            {
                return false;
            }
        }
    }
    return true;
}

bool isHermitian(const Eigen::MatrixXcd &matrix)
{
    return matrix.isApprox(matrix.adjoint(), 1E-10);
}

void Diagonalization::roundMatrixHamiltonian()
{
    int numDigits = 6;
    double factor = std::pow(10.0, numDigits); // Compute factor for rounding

    for (int row = 0; row < matrixHamiltonian.rows(); row++)
    {
        for (int col = 0; col <= row; col++)
        {
            // Round the real part
            double realPart = std::round(matrixHamiltonian(row, col).real() * factor) / factor;

            // Round the imaginary part
            double imagPart = std::round(matrixHamiltonian(row, col).imag() * factor) / factor;

            // Update the lower triangular part
            matrixHamiltonian(row, col) = std::complex<double>(realPart, imagPart);

            // Update the corresponding upper triangular part to ensure Hermitian property
            if (row != col)
            {
                matrixHamiltonian(col, row) = std::complex<double>(realPart, -imagPart);
            }
        }
    }
}

bool isPositiveDefinite(const Eigen::MatrixXcd &matrix)
{
    if (!isHermitian(matrix))
    {
        std::cout << "positive definite check: matrix is not hermitian" << std::endl;
        return false;
    }

    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(matrix);

    for (auto eigenvalue : solver.eigenvalues())
    {
        if (eigenvalue.real() <= 0)
        {
            std::cout << "positive definite check: eigenvalue smaller than zero" << std::endl;
            // std::cout << "Eigenvalues: " << solver.eigenvalues() << std::endl;
            return false;
        }
    }
    return true;
}

bool isApproxZero(const Eigen::MatrixXcd &matrix, double epsilon = 1e-10)
{
    // Check if all elements of the matrix are approximately zero
    return matrix.cwiseAbs().maxCoeff() <= epsilon;
}

struct EigenPair
{
    double eigenvalue;
    Eigen::MatrixXcd eigenvector;
};

bool compareEigenPairs(const EigenPair &a, const EigenPair &b)
{
    return a.eigenvalue > b.eigenvalue;
}

/// @brief Sorts the eigenvalues and eigenvectors in ascending order
/// @param eigenvectors matrix with eigenvectors as columns
/// @param eigenvalues  diagonal matrix with eigenvalues on the diagonal
void sortEigenvaluesAndVectors(Eigen::MatrixXcd &eigenvectors, Eigen::MatrixXcd &eigenvalues)
{
    // Create a vector of eigenpairs
    std::vector<EigenPair> eigenPairs;
    for (int i = 0; i < eigenvalues.cols(); ++i)
    {
        eigenPairs.push_back({eigenvalues(i, i).real(), eigenvectors.col(i)});
    }

    // Sort eigenpairs by eigenvalue
    std::sort(eigenPairs.begin(), eigenPairs.end(), compareEigenPairs);

    // Create sorted eigenvalue vector and eigenvector matrix
    Eigen::MatrixXcd sortedEigenvalues = Eigen::MatrixXcd::Zero(eigenvalues.rows(), eigenvalues.cols());
    Eigen::MatrixXcd sortedEigenvectors(eigenvectors.rows(), eigenvectors.cols());

    for (int i = 0; i < eigenPairs.size(); ++i)
    {
        sortedEigenvalues(i, i) = eigenPairs[i].eigenvalue;
        sortedEigenvectors.col(i) = eigenPairs[i].eigenvector;
    }

    eigenvalues = sortedEigenvalues;
    eigenvectors = sortedEigenvectors;
}

ComputationInfo Diagonalization::diagonalizeColpa()
{
    // assertions:
    // - matrixHamiltonian is hermitian
    // - matrixHamiltonian is positive definite
    // - U is unitary
    // - L is diagonal and has positive and then negative values on the diagonal
    // - Q^dagger H Q is diagonal

    // roundMatrixHamiltonian();

    if (!isHermitian(matrixHamiltonian))
    {
        std::cout << "Matrix is not hermitian" << std::endl;
        return ComputationInfo::FAILURE;
    }

    if (!isPositiveDefinite(matrixHamiltonian))
    {
        std::cout << "Matrix is not positive definite" << std::endl;
        return ComputationInfo::FAILURE;
    }

    Eigen::LLT<Eigen::MatrixXcd> lltOfMatrixHamiltonian(matrixHamiltonian);

    if (lltOfMatrixHamiltonian.info() != Eigen::Success)
    {
        std::cout << "Cholesky decomposition failed" << std::endl;
        return ComputationInfo::FAILURE;
    }

    Eigen::MatrixXcd K = lltOfMatrixHamiltonian.matrixL().adjoint();
    Eigen::MatrixXcd Kdagger = K.adjoint();

    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(K * getMatrix_g() * Kdagger);

    if (solver.info() != Eigen::Success)
    {
        std::cout << "Eigenvalue decomposition of K g H^dagger failed" << std::endl;
        return ComputationInfo::FAILURE;
    }

    // use QR-decomposition to get unitary matrix
    Eigen::MatrixXcd V = solver.eigenvectors();
    Eigen::HouseholderQR<Eigen::MatrixXcd> qr(V);
    Eigen::MatrixXcd qrR = qr.matrixQR().triangularView<Eigen::Upper>();
    Eigen::MatrixXcd qrQ = qr.householderQ();

    if (!isApproxZero(qrQ.adjoint() * qrQ - Eigen::MatrixXcd::Identity(8, 8), 1E-10))
    {
        std::cout << "Q is not unitary" << std::endl;
        std::cout << "max distance to zero: " << (qrQ * qrQ.adjoint() - Eigen::MatrixXcd::Identity(8, 8)).cwiseAbs().maxCoeff() << std::endl;
        // return ComputationInfo::FAILURE;
    }

    if (!isApproxZero(V - qrQ * qrR, 1E-10))
    {
        std::cout << "QR decomposition not successfull" << std::endl;
        // std::cout << U - Q * R << std::endl;
        std::cout << "max dist to zero " << (V - qrQ * qrR).cwiseAbs().maxCoeff() << std::endl;
        // return ComputationInfo::FAILURE;
    }

    Eigen::MatrixXcd U = qrQ;
    Eigen::MatrixXcd Udagger = U.adjoint();

    Eigen::MatrixXcd L = qrR * solver.eigenvalues().asDiagonal() * qrR.inverse();

    if (!isDiagonal(L))
    {
        std::cout << "L is not diagonal" << std::endl;
        return ComputationInfo::FAILURE;
    }

    if (!isDiagonal(Udagger * K * getMatrix_g() * Kdagger * U))
    {
        std::cout << "Udagger * K * g * Kdagger * U is not diagonal" << std::endl;
        return ComputationInfo::FAILURE;
    }

    sortEigenvaluesAndVectors(U, L);

    Udagger = U.adjoint();

    if (!isApproxZero(Udagger * K * getMatrix_g() * Kdagger * U - L, 1E-10))
    {
        std::cout << "Udagger * K * g * Kdagger * U is not equal to L" << std::endl;
        return ComputationInfo::FAILURE;
    }

    // check if U is unitary
    if (!isApproxZero(Udagger * U - Eigen::MatrixXcd::Identity(8, 8), 1E-10))
    {
        std::cout << "Udagger * U is not equal to the identity matrix" << std::endl;
        return ComputationInfo::FAILURE;
    }

    Eigen::MatrixXcd finalDiagonalMatrix = getMatrix_g() * L;
    Eigen::MatrixXcd T = (K.inverse() * U * squaredMatrix(finalDiagonalMatrix)).inverse();

    // eigenvectors = (K.inverse() * U * squaredMatrix(finalDiagonalMatrix)).inverse();
    // eigenvectors_inverse = eigenvectors.inverse();

    // check if T can diagonalize the Hamiltonian
    if (!isDiagonal(T.adjoint().inverse() * matrixHamiltonian * T.inverse()))
    {
        std::cout << "Qdagger H Q is not diagonal. Largest off-diag elem: " << getLargestNonDiagonalElement(eigenvectors_inverse * matrixHamiltonian * eigenvectors) << std::endl;
        printMatrix(eigenvectors.adjoint() * matrixHamiltonian * eigenvectors);
        return ComputationInfo::FAILURE;
    }

    for (auto eigenvalue : finalDiagonalMatrix.diagonal())
    {
        eigenEnergies.push_back(std::abs(eigenvalue.real()));
    }

    // correct notation the program's rest
    eigenvectors = T.inverse();
    eigenvectors_inverse = T;

    if (!isDiagonal(eigenvectors.adjoint() * matrixHamiltonian * eigenvectors))
    {
        std::cout << "Qdagger H Q is not diagonal" << std::endl;
        return ComputationInfo::FAILURE;
    }

    if (!isApproxZero(eigenvectors.adjoint() * matrixHamiltonian * eigenvectors - finalDiagonalMatrix, 1E-10))
    {
        std::cout << " Qdagger H Q != E " << std::endl;
        return ComputationInfo::FAILURE;
    }

    return SUCCESS;
}


void Diagonalization::calcAngularMomentum()
{

    Eigen::MatrixXcd matrixL = Eigen::MatrixXcd::Zero(8, 8);
    std::complex<double> i(0, 1);

    double polVec = phDisp.polVectors[2][2]; // phDisp.polVectors[1][2]; //

    assert(matrixL.rows() == 8);
    assert(matrixL.cols() == 8);
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

    for (int j = 0; j < 8; j++)
    {
        angularMomentum.push_back(D(j, j));
    }

    return;
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
