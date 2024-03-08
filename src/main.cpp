#include <complex>
// #include <fftw3.h>
#include "../include/couplingParameters.h"
#include "../include/fourierTransform.h"
#include "../include/path.h"
#include "../include/eigen.h"
#include "../include/util.h"
#include "../include/dispersion.h"
#include "../include/diagonalization.h"
#include "../include/dynamics.h"

int main()
{


    




    std::vector<Vector3D> path = constructPath(100, 1);

    const double ATOMIC_MASS = 55.845; // atomic mass given in Dalton
    const double S = 1.1;

    // Calculates phonon dispersion relation
    std::vector<PhononDispParam> phononDispersion = getPhononDispersion("Parameters/dynMat_16x16x16.txt", "Parameters/nn4.txt", "Outputs/numbersPhTest.txt", path, ATOMIC_MASS);

    // Calculate magnon dispersion relation
    std::vector<MagnonDispParam> magnonDispersion = getMagneticDispersion("Parameters/J_bccFe.txt", "Outputs/numbersJIso_test.txt", path, S);

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

    // parameters.insert(parameters.end(), ij_uk_x_parameter.begin(), ij_uk_x_parameter.end());
    // parameters.insert(parameters.end(), ij_uk_y_parameter.begin(), ij_uk_y_parameter.end());
    // parameters.insert(parameters.end(), ij_uk_z_parameter.begin(), ij_uk_z_parameter.end());

    parameters.insert(parameters.begin(), ij_uk_x_parameter.begin(), ij_uk_x_parameter.end());
    parameters.insert(parameters.begin(), ij_uk_y_parameter.begin(), ij_uk_y_parameter.end());
    parameters.insert(parameters.begin(), ij_uk_z_parameter.begin(), ij_uk_z_parameter.end());


    Coefficients coeff;

    coeff.initReciprocalLatticeVectors();
    coeff.init(parameters);


    /*
    // Diagonalisation
    std::ofstream outFileEV("Outputs/8x8Eigenenergies.txt");
    std::ofstream outFileCD("Outputs/8x8CD.txt");
    std::ofstream outFileEVectors("Outputs/8x8EVec.txt");

    for (int idx = 0; idx < path.size(); idx++)
    {
        std::cout << "Progress: " << idx/double(path.size()) << "\n"; 

        Diagonalization diag(parameters, phononDispersion.at(idx), magnonDispersion.at(idx), path.at(idx), ATOMIC_MASS, S);
        diag.calcCD();
        diag.calcMatrixHamiltonian();
        diag.diagonalize();

        outFileEV << diag.k.x << "," << diag.k.y << "," << diag.k.z << ",";
        outFileCD << diag.k.x << "," << diag.k.y << "," << diag.k.z << ",";
        outFileEVectors << diag.k.x << "," << diag.k.y << "," << diag.k.z << ",";

        for (int i = 0; i < 7; i++)
        {
            outFileEV << diag.eigenEnergies.at(i) << ",";
        }
        outFileEV << diag.eigenEnergies.at(7) << "\n";

        outFileCD << diag.C.at(0) << "," << diag.C.at(1) << "," << diag.C.at(2) << "," << diag.D.at(0) << "," << diag.D.at(1) << "," << diag.D.at(2) << "\n";

        for (int col = 0; col < 8; col++)
        {
            for (int row = 0; row < 8; row++)
            {
                outFileEVectors << diag.eigenvectors_inverse.col(col).row(row).x();
                if (row == 7 && col == 7)
                {
                    outFileEVectors << "\n";
                }
                else
                {
                    outFileEVectors << ",";
                }
            }
        }
    }
    outFileEVectors.close();
    outFileEV.close();
    outFileCD.close();
    */

    /*
    std::vector<std::vector<double>> allEigenenergies = diagonalizeHamiltonian(path, phononDispersion, magnonDispersion, parameters);

    std::ofstream outFileEV("Outputs/8x8Eigenenergies.txt");

    int counter = 0;
    for (auto eigenenergies : allEigenenergies)
    {

        outFileEV << path.at(counter).x << "," << path.at(counter).y << "," << path.at(counter).z << ",";

        if (eigenenergies.size() == 0)
        {
            continue;
        }
        for (int i = 0; i < 7; i++)
        {
            outFileEV << eigenenergies.at(i) << ",";
        }
        outFileEV << eigenenergies.at(7) << "\n";
        counter++;
    }

    outFileEV.close();

    return 0;
    */

    std::ofstream outFile("Outputs/numbersDTest.txt");

    if (!outFile.is_open())
    {
        std::cerr << "Failed to open the file for writing." << std::endl;
        return 1; // Exit with an error code
    }

    outFile << "kx,ky,kz,D_k_x_x,D_k_x_y,D_k_x_z,D_k_y_x,D_k_y_y,D_k_y_z\n";

    for (Vector3D k : path)
    {
        std::complex<double> D_k_x_x = FTD(k.x, k.y, k.z, parametersX, X, X);
        std::complex<double> D_k_x_y = FTD(k.x, k.y, k.z, parametersY, X, Y);
        std::complex<double> D_k_x_z = FTD(k.x, k.y, k.z, parametersZ, X, Z);

        std::complex<double> D_k_y_x = FTD(k.x, k.y, k.z, parametersX, Y, X);
        std::complex<double> D_k_y_y = FTD(k.x, k.y, k.z, parametersY, Y, Y);
        std::complex<double> D_k_y_z = FTD(k.x, k.y, k.z, parametersZ, Y, Z);

        // calculate the FT of the DMI-like coupling parameters with k vectors having opposite sign
        std::complex<double> D_mk_x_x = FTD(-k.x, -k.y, -k.z, parametersX, X, X);
        std::complex<double> D_mk_x_y = FTD(-k.x, -k.y, -k.z, parametersY, X, Y);
        std::complex<double> D_mk_x_z = FTD(-k.x, -k.y, -k.z, parametersZ, X, Z);

        std::complex<double> D_mk_y_x = FTD(-k.x, -k.y, -k.z, parametersX, Y, X);
        std::complex<double> D_mk_y_y = FTD(-k.x, -k.y, -k.z, parametersY, Y, Y);
        std::complex<double> D_mk_y_z = FTD(-k.x, -k.y, -k.z, parametersZ, Y, Z);

        std::complex<double> i(0, 1);

        // std:: cout << std::conj((D_k_x_x - i * D_k_y_x)) - (D_k_x_x + i * D_k_y_x) << std::endl;

        outFile << k.x << "," << k.y << "," << k.z << "," << D_k_x_x << "," << D_k_x_y << "," << D_k_x_z << "," << D_k_y_x << "," << D_k_y_y << "," << D_k_y_z << "," << D_mk_x_x << "," << D_mk_x_y << "," << D_mk_x_z << "," << D_mk_y_x << "," << D_mk_y_y << "," << D_mk_y_z << "\n";
    }

    outFile.close();

    /*
    // Calculates the magnon-magnon interaction DOS by sampling the BZ
    std::ofstream outFileBZSampling("Outputs/numbersBZSampling.txt");
    outFileBZSampling << "kx,ky,kz,J\n";

    std::vector<Vector3D> grid = sampleBZ(100);
    for (Vector3D k : grid)
    {
        std::complex<double> Jk = FTJiso(k.x, k.y, k.z, parameters_J_iso);
        outFileBZSampling << k.x << "," << k.y << "," << k.z << "," << Jk << "\n";
    }
    outFileBZSampling.close();
    */

    return 0;
}
