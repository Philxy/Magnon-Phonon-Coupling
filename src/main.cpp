#include <complex>
// #include <fftw3.h>
#include "../include/couplingParameters.h"
#include "../include/fourierTransform.h"
#include "../include/path.h"
#include "../include/eigen.h"
#include <omp.h> // Include the OpenMP header

int main()
{
    // Calculates phonon-phonon dispersion relation
    std::vector<CouplingParameter> dyn_matrices = readDynMatrices("Parameters/dynMat_16x16x16.txt"); // retrieve the dynamical matrices
    std::vector<CouplingParameter> next_neighbors = readNextNeighbors("Parameters/nn4.txt");         // get the nearest neighbors in real space

    std::ofstream outFileFC("Outputs/FC.txt"); // file to write the force cosntants to
    outFileFC << "x,y,z,Phi_xx,Phi_xy,Phi_xz,Phi_yx,Phi_yy,Phi_yz,Phi_zx,Phi_zy,Phi_zz\n";

    for (const CouplingParameter &nn : next_neighbors)
    {
        Eigen::Matrix3d force_mat = forceMatrix(nn.x, nn.y, nn.z, dyn_matrices);
        outFileFC << nn.x << "," << nn.y << "," << nn.z << "," << force_mat(0, 0) << "," << force_mat(0, 1) << "," << force_mat(0, 2) << "," << force_mat(1, 0) << "," << force_mat(1, 1) << "," << force_mat(1, 2) << "," << force_mat(2, 0) << "," << force_mat(2, 1) << "," << force_mat(2, 2) << "\n";
    }

    outFileFC.close();

    std::vector<CouplingParameter> parameters_ph = readCouplingParametersPh("Outputs/FC.txt");

    std::ofstream outFilePh("Outputs/numbersPh.txt");
    outFilePh << "kx,ky,kz,omega1,omega2,omega3\n";

    for (const Vector3D &k : constructPath(100, 1))
    {
        Eigen::Matrix3d dynMat_k = dynMat(k.x, k.y, k.z, parameters_ph);
        Eigen::EigenSolver<Eigen::Matrix3d> solver(dynMat_k);
        Eigen::Vector3cd eigenvalues = solver.eigenvalues();
        Eigen::Matrix3cd eigenvectors = solver.eigenvectors();

        outFilePh << k.x << "," << k.y << "," << k.z << "," << eigenvalues.x().real() << "," << eigenvalues.y().real() << "," << eigenvalues.z().real() << "\n";
    }
    outFilePh.close();
    return 0;

    // Calculates the magnon-magnon dispersion relation and writes it to a file based on the isotropic Heisenberg exchange between neighbors given in a file
    std::vector<CouplingParameter> parameters_J_iso = readCouplingParametersIso("Parameters/J_bccFe.txt");

    std::ofstream outFileIso("Outputs/numbersIso.txt");
    outFileIso << "kx,ky,kz,J\n";

    for (const Vector3D &k : constructPath(200, 1))
    {
        std::complex<double> Jk = FTJiso(k.x, k.y, k.z, parameters_J_iso);
        outFileIso << k.x << "," << k.y << "," << k.z << "," << Jk << "\n";
    }
    outFileIso.close();

    /// Calculates the magnon-magnon interaction DOS by sampling the BZ
    std::ofstream outFileBZSampling("Outputs/numbersBZSampling.txt");
    outFileBZSampling << "kx,ky,kz,J\n";

    std::vector<Vector3D> grid = sampleBZ(100);
    for (Vector3D k : grid)
    {
        std::complex<double> Jk = FTJiso(k.x, k.y, k.z, parameters_J_iso);
        outFileBZSampling << k.x << "," << k.y << "," << k.z << "," << Jk << "\n";
    }
    outFileBZSampling.close();

    std::vector<CouplingParameter> parameters;  // contains all the parameters
    std::vector<CouplingParameter> parametersX; // contains all the parameters with x displacement
    std::vector<CouplingParameter> parametersY; // contains all the parameters with y displacement
    std::vector<CouplingParameter> parametersZ; // contains all the parameters with z displacement

    std::vector<CouplingParameter> ij_uk_x_parameter = readCouplingParameters("Parameters/SLC_Eisen/Fe_full_tensor_ij-uk_x_displacement.csv", X);
    std::vector<CouplingParameter> ij_uk_y_parameter = readCouplingParameters("Parameters/SLC_Eisen/Fe_full_tensor_ij-uk_y_displacement.csv", Y);
    std::vector<CouplingParameter> ij_uk_z_parameter = readCouplingParameters("Parameters/SLC_Eisen/Fe_full_tensor_ij-uk_z_displacement.csv", Z);

    // Add the parameters for the x,y,z direction to their respective vectors:
    parametersX.insert(parametersX.end(), ij_uk_x_parameter.begin(), ij_uk_x_parameter.end());
    parametersY.insert(parametersY.end(), ij_uk_y_parameter.begin(), ij_uk_y_parameter.end());
    parametersZ.insert(parametersZ.end(), ij_uk_z_parameter.begin(), ij_uk_z_parameter.end());

    std::ofstream outFile("Outputs/numbers.txt");

    if (!outFile.is_open())
    {
        std::cerr << "Failed to open the file for writing." << std::endl;
        return 1; // Exit with an error code
    }

    outFile << "kx,ky,kz,D_k_x_x,D_k_x_y,D_k_x_z,D_k_y_x,D_k_y_y,D_k_y_z\n";

    for (const Vector3D &k : constructPath(100, 1))
    {
        std::complex<double> D_k_x_x = FTD(k.x, k.y, k.z, parametersX, X, X);
        std::complex<double> D_k_x_y = FTD(k.x, k.y, k.z, parametersY, X, Y);
        std::complex<double> D_k_x_z = FTD(k.x, k.y, k.z, parametersZ, X, Z);

        std::complex<double> D_k_y_x = FTD(k.x, k.y, k.z, parametersX, Y, X);
        std::complex<double> D_k_y_y = FTD(k.x, k.y, k.z, parametersY, Y, Y);
        std::complex<double> D_k_y_z = FTD(k.x, k.y, k.z, parametersZ, Y, Z);

        // std::complex<double> D_k_z_x = FT(k.x, k.y, k.z, parameters, Z, X);
        // std::complex<double> D_k_z_y = FT(k.x, k.y, k.z, parameters, Z, Y);
        // std::complex<double> D_k_z_z = FT(k.x, k.y, k.z, parameters, Z, Z);

        outFile << k.x << "," << k.y << "," << k.z << "," << D_k_x_x << "," << D_k_x_y << "," << D_k_x_z << "," << D_k_y_x << "," << D_k_y_y << "," << D_k_y_z << "\n";
    }

    // Close the file
    outFile.close();

    return 0;
}
