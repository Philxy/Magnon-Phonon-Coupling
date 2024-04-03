#include "../include/dispersion.h"

// Calculates the magnon-magnon dispersion relation and writes it to a file based on the isotropic Heisenberg exchange between neighbors given in a file
std::vector<MagnonDispParam> getMagneticDispersion(std::string couplingParameterFile, std::string magnonDispOutputPath, std::vector<Vector3D> path)
{
    // units: The isotropic Heisenberg exchange parameters are assumed to be in mRy and they are converted to meV
    const double conversionFactor = 13.606;
    const double d_aniso = 6.97 * 1E-3;

    std::vector<CouplingParameter> parameters_J_iso = readCouplingParametersIso(couplingParameterFile);

    std::vector<MagnonDispParam> magnonDisp;

    std::ofstream outFileIso(magnonDispOutputPath);
    outFileIso << "kx,ky,kz,J\n";

    for (const Vector3D &k : path)
    {
        std::complex<double> Jk = 1.0 / S * (FTJiso(k.x, k.y, k.z, parameters_J_iso).real() * conversionFactor + d_aniso);
        outFileIso << k.x << "," << k.y << "," << k.z << "," << Jk << "\n";

        MagnonDispParam E;
        E.kx = k.x;
        E.ky = k.y;
        E.kz = k.z;
        E.energy = Jk.real();
        magnonDisp.push_back(E);
    }
    outFileIso.close();

    return magnonDisp;
}

std::vector<MagnonDispParam> getMagneticDispersion(std::string couplingParameterFile, std::vector<Vector3D> path)
{
    // units: The isotropic Heisenberg exchange parameters are assumed to be in mRy and they are converted to meV
    const double conversionFactor = 13.606;
    const double d_aniso = 6.97 * 1E-3;

    std::vector<CouplingParameter> parameters_J_iso = readCouplingParametersIso(couplingParameterFile);

    std::vector<MagnonDispParam> magnonDisp;

    for (const Vector3D &k : path)
    {
        std::complex<double> Jk = 1.0 / S * (FTJiso(k.x, k.y, k.z, parameters_J_iso).real() * conversionFactor + d_aniso);

        MagnonDispParam E;
        E.kx = k.x;
        E.ky = k.y;
        E.kz = k.z;
        E.energy = Jk.real();
        magnonDisp.push_back(E);
    }

    return magnonDisp;
}

// Calculates phonon-phonon dispersion relation
std::vector<PhononDispParam> getPhononDispersion(std::string dynamicMatricesFile, std::string nextNeighbourFile, std::string phononDispOutputPath, std::vector<Vector3D> path)
{
    // Units:
    // Calculate the eigenenergies from the eigenvalues in meV
    // assuming the eigenvalues have units of mRy/(a.u.)^2

    std::vector<CouplingParameter> dyn_matrices = readDynMatrices(dynamicMatricesFile);   // retrieve the dynamical matrices
    std::vector<CouplingParameter> next_neighbors = readNextNeighbors(nextNeighbourFile); // get the nearest neighbors in real space

    std::vector<CouplingParameter> parameters_ph;
    std::vector<PhononDispParam> phononDisp;

    for (const CouplingParameter &nn : next_neighbors)
    {
        Eigen::Matrix3d force_mat = forceMatrix(nn.x, nn.y, nn.z, dyn_matrices);
        CouplingParameter p;
        p.x = nn.x;
        p.y = nn.y;
        p.z = nn.z;
        for (int row = 0; row < 3; row++)
        {
            for (int col = 0; col < 3; col++)
            {
                p.Phi[row][col] = force_mat.row(row).col(col).x();
            }
        }
        parameters_ph.push_back(p);
    }

    std::ofstream outFilePh(phononDispOutputPath);
    outFilePh << "kx,ky,kz,omega1,omega2,omega3,e1x,e1y,e1z,e2x,e2y,e2z,e3x,e3y,e3z\n";

    std::vector<Eigen::Vector3cd> all_eigenvalues;

    for (const Vector3D &k : path)
    {
        // Solve the eigenvalue problem
        Eigen::Matrix3d dynMat_k = dynMat(k.x, k.y, k.z, parameters_ph);
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(dynMat_k);
        Eigen::Vector3cd eigenvalues = solver.eigenvalues();
        Eigen::Matrix3cd eigenvectors = solver.eigenvectors();

        // makeEigenvaluesPositive(eigenvalues, eigenvectors);
        changeEigenVecSign(eigenvectors);

        double E1 = sqrt(abs(eigenvalues.x().real()) / atomicMass) * 14.25133727;
        double E2 = sqrt(abs(eigenvalues.y().real()) / atomicMass) * 14.25133727;
        double E3 = sqrt(abs(eigenvalues.z().real()) / atomicMass) * 14.25133727;

        PhononDispParam phDispParam;
        phDispParam.kx = k.x;
        phDispParam.ky = k.y;
        phDispParam.kz = k.z;

        phDispParam.E[0] = E1;
        phDispParam.E[1] = E2;
        phDispParam.E[2] = E3;

        for (int row = 0; row < 3; row++)
        {
            for (int col = 0; col < 3; col++)
            {
                phDispParam.polVectors[row][col] = eigenvectors.row(row).col(col).x().real();
            }
        }

        phononDisp.push_back(phDispParam);

        // sortEigen(eigenvalues, eigenvectors);

        outFilePh << k.x << "," << k.y << "," << k.z << "," << E1 << "," << E2 << "," << E3 << "," << eigenvectors.col(0).x().real() << "," << eigenvectors.col(0).y().real() << "," << eigenvectors.col(0).z().real() << "," << eigenvectors.col(1).x().real() << "," << eigenvectors.col(1).y().real() << "," << eigenvectors.col(1).z().real() << "," << eigenvectors.col(2).x().real() << "," << eigenvectors.col(2).y().real() << "," << eigenvectors.col(2).z().real() << "\n";
    }

    outFilePh.close();
    return phononDisp;
}

std::vector<PhononDispParam> getPhononDispersion(std::string dynamicMatricesFile, std::string nextNeighbourFile, std::vector<Vector3D> path)
{
    // Units:
    // Calculate the eigenenergies from the eigenvalues in meV
    // assuming the eigenvalues have units of mRy/(a.u.)^2

    std::vector<CouplingParameter> dyn_matrices = readDynMatrices(dynamicMatricesFile);   // retrieve the dynamical matrices
    std::vector<CouplingParameter> next_neighbors = readNextNeighbors(nextNeighbourFile); // get the nearest neighbors in real space

    std::vector<CouplingParameter> parameters_ph;
    std::vector<PhononDispParam> phononDisp;

    for (const CouplingParameter &nn : next_neighbors)
    {
        Eigen::Matrix3d force_mat = forceMatrix(nn.x, nn.y, nn.z, dyn_matrices);
        CouplingParameter p;
        p.x = nn.x;
        p.y = nn.y;
        p.z = nn.z;
        for (int row = 0; row < 3; row++)
        {
            for (int col = 0; col < 3; col++)
            {
                p.Phi[row][col] = force_mat.row(row).col(col).x();
            }
        }
        parameters_ph.push_back(p);
    }

    std::vector<Eigen::Vector3cd> all_eigenvalues;

    for (const Vector3D &k : path)
    {
        // Solve the eigenvalue problem
        Eigen::Matrix3d dynMat_k = dynMat(k.x, k.y, k.z, parameters_ph);
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(dynMat_k);
        Eigen::Vector3cd eigenvalues = solver.eigenvalues();
        Eigen::Matrix3cd eigenvectors = solver.eigenvectors();

        // makeEigenvaluesPositive(eigenvalues, eigenvectors);
        changeEigenVecSign(eigenvectors);

        double E1 = sqrt(abs(eigenvalues.x().real()) / atomicMass) * 14.25133727;
        double E2 = sqrt(abs(eigenvalues.y().real()) / atomicMass) * 14.25133727;
        double E3 = sqrt(abs(eigenvalues.z().real()) / atomicMass) * 14.25133727;

        PhononDispParam phDispParam;
        phDispParam.kx = k.x;
        phDispParam.ky = k.y;
        phDispParam.kz = k.z;

        phDispParam.E[0] = E1;
        phDispParam.E[1] = E2;
        phDispParam.E[2] = E3;

        for (int row = 0; row < 3; row++)
        {

            for (int col = 0; col < 3; col++)
            {
                phDispParam.polVectors[row][col] = eigenvectors.row(row).col(col).x().real();
            }
        }

        phononDisp.push_back(phDispParam);

        // sortEigen(eigenvalues, eigenvectors);
    }

    return phononDisp;
}

// Read phonon dispersion parameters from a file
// Assumptions:
// - k vector given in units of 1/a
// - Energies given in meV
std::vector<PhononDispParam> readPhononDispParams(const std::string &filePath)
{
    std::vector<PhononDispParam> params;
    std::ifstream file(filePath);

    if (!file.is_open())
    {
        std::cerr << "Failed to open file: " << filePath << std::endl;
        return params; // Return an empty vector if the file couldn't be opened
    }

    std::string line;

    std::getline(file, line); // skip the first line

    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        PhononDispParam param;
        // Read kx, ky, kz
        iss >> param.kx >> param.ky >> param.kz;
        // Read energies
        for (int i = 0; i < 3; ++i)
        {
            iss >> param.E[i];
        }
        // Read polarization vectors
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                iss >> param.polVectors[j][i]; // Notice the indexing: [j][i] to correctly fill the columns
            }
        }

        // param.E[0] *= 10;
        // param.E[1] *= 10;
        // param.E[2] *= 10;

        // Shift the energies 
        //for (int i = 0; i < 3; i++)
        //{
        //    param.E[i] += 0.1;
        //}

        params.push_back(param);
    }

    file.close();
    return params;
}

// Given a phonon dispersion, calculate the magnetic dispersion
// using the provided isotropic Heisenberg exchange parameters
// and the k vectors from the phonon dispersion
std::vector<MagnonDispParam> getMagnonDispFromPhononDisp(const std::vector<PhononDispParam> &phDisp, std::string couplingParameterFile, std::string magnonDispOutputPath)
{
    std::vector<MagnonDispParam> magnonDisp;

    std::ofstream outFileIso(magnonDispOutputPath);
    outFileIso << "kx,ky,kz,J\n";

    // units: The isotropic Heisenberg exchange parameters are assumed to be in mRy and they are converted to meV
    const double conversionFactor = 13.606;
    const double d_aniso = 6.97 * 1E-3;
    std::vector<CouplingParameter> parameters_J_iso = readCouplingParametersIso(couplingParameterFile);

    for (const PhononDispParam &ph : phDisp)
    {
        MagnonDispParam mag;
        mag.kx = ph.kx;
        mag.ky = ph.ky;
        mag.kz = ph.kz;

        std::complex<double> Jk = 1.0 / S * (FTJiso(ph.kx, ph.ky, ph.kz, parameters_J_iso).real() * conversionFactor + d_aniso);
        outFileIso << ph.kx << "," << ph.ky << "," << ph.kz << "," << Jk << "\n";

        mag.energy = Jk.real();
        magnonDisp.push_back(mag);
    }

    outFileIso.close();

    return magnonDisp;
}

bool checkZeroEnergy(const PhononDispParam &phononDispersion, const Vector3D &kVec)
{
    for (int i = 0; i < 3; i++)
    {
        if (std::abs(phononDispersion.E[i]) < std::numeric_limits<double>::epsilon())
        {
            std::cout << "Zero energy found at k: " << kVec.x << " " << kVec.y << " " << kVec.z << std::endl;
            return true;
        }
    }
    return false;
}