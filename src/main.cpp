#include <complex>
// #include <fftw3.h>
#include "../include/couplingParameters.h"
#include "../include/fourierTransform.h"
#include "../include/path.h"
#include "../include/eigen.h"

int main()
{



    ///*
    // Eigen "Hello, World" program to determine eigenvalues and eigenvectors of a given matrix
    using Eigen::Matrix3d;
    using Eigen::Vector3d;

    // Define a 3x3 matrix
    Eigen::Matrix3d matrix;

    matrix(0,0) = 1;
    matrix(0,1) = 2;
    matrix(0,2) = 3;
    matrix(1,0) = 4;
    matrix(1,1) = 5;
    matrix(1,2) = 6;
    matrix(2,0) = 7;
    matrix(2,1) = 8;
    matrix(2,2) = 9;

    // Compute the eigenvalues and eigenvectors
    Eigen::EigenSolver<Eigen::Matrix3d> solver(matrix);

    // Print the eigenvalues
    std::cout << "The eigenvalues of the matrix are:\n"
              << solver.eigenvalues() << std::endl;
    // Print the eigenvectors
    std::cout << "The eigenvectors of the matrix are:\n"
              << solver.eigenvectors() << std::endl;

    return 0;
    //*/


    // Calculates the magnon-magnon dispersion relation and writes it to a file based on the isotropic Heisenberg exchange between neighbors given in a fiel
    std::vector<CouplingParameter> parameters_J_iso = readCouplingParametersIso("Parameters/J_bccFe.txt");

    std::ofstream outFileIso("Outputs/numbersIso.txt");
    outFileIso << "kx,ky,kz,J\n";

    for (Vector3D k : constructPath(200, 1))
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

    for (Vector3D k : constructPath(100, 1))
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
