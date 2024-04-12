#include <complex>
#include <sstream>
// #include <fftw3.h>
#include "../include/couplingParameters.h"
#include "../include/diagonalization.h"
#include "../include/dynamics.h"

#include <chrono>

int main()
{

    // Read the magnon-phonon coupling parameters from the files

    std::vector<CouplingParameter> parameters, parametersX, parametersY, parametersZ; // vectors containing all the parameters or the parameters with x,y,z displacement

    std::vector<CouplingParameter> ij_uk_x_parameter = readCouplingParameters("Parameters/SLC_Eisen/Fe_full_tensor_ij-uk_x_displacement.csv", X);
    std::vector<CouplingParameter> ij_uk_y_parameter = readCouplingParameters("Parameters/SLC_Eisen/Fe_full_tensor_ij-uk_y_displacement.csv", Y);
    std::vector<CouplingParameter> ij_uk_z_parameter = readCouplingParameters("Parameters/SLC_Eisen/Fe_full_tensor_ij-uk_z_displacement.csv", Z);

    // Add the parameters for the x,y,z direction to their respective vectors:
    parametersX.insert(parametersX.end(), ij_uk_x_parameter.begin(), ij_uk_x_parameter.end());
    parametersY.insert(parametersY.end(), ij_uk_y_parameter.begin(), ij_uk_y_parameter.end());
    parametersZ.insert(parametersZ.end(), ij_uk_z_parameter.begin(), ij_uk_z_parameter.end());
    // Add all parameters to the parameters vector
    parameters.insert(parameters.end(), ij_uk_x_parameter.begin(), ij_uk_x_parameter.end());
    parameters.insert(parameters.end(), ij_uk_y_parameter.begin(), ij_uk_y_parameter.end());
    parameters.insert(parameters.end(), ij_uk_z_parameter.begin(), ij_uk_z_parameter.end());

    ///*

    // _____________ Diagonalize the Hamiltonian _____________

    std::vector<PhononDispParam> phononDispersion = readPhononDispParams("scripts/Data/gH1000/path.txt");
    std::vector<MagnonDispParam> magnonDispersion = getMagnonDispFromPhononDisp(phononDispersion, "Parameters/J_bccFe.txt", "Outputs/numbersJIso.txt");

    // Pre-allocate space for output data, initializing with empty strings or appropriate default values
    std::vector<std::string> outEV(phononDispersion.size());
    std::vector<std::string> outCD(phononDispersion.size());
    std::vector<std::string> outEVectors(phononDispersion.size());
    std::vector<std::string> outAngularMomentum(phononDispersion.size());

#pragma omp parallel for
    for (int idx = 0; idx < phononDispersion.size(); idx++)
    {
        Vector3D kVec(phononDispersion.at(idx).kx, phononDispersion.at(idx).ky, phononDispersion.at(idx).kz);

        // do not diagonalize if the phonon energy is zero
        if (checkZeroEnergy(phononDispersion.at(idx), kVec))
        {
            continue;
        }

        Diagonalization diag(parameters, phononDispersion.at(idx), magnonDispersion.at(idx), kVec);
        diag.calcCD();
        diag.calcMatrixHamiltonian();
        diag.diagonalize();
        diag.calcAngularMomentum(0.0);

        // Some debug output
        // std::cout << "k: \n";
        // std::cout << phononDispersion.at(idx).kx << " " << phononDispersion.at(idx).ky << " " << phononDispersion.at(idx).kz << std::endl;
        // std::cout << "Phonon energies: \n";
        // std::cout << phononDispersion.at(idx).E[0] << " " << phononDispersion.at(idx).E[1] << " " << phononDispersion.at(idx).E[2] << std::endl;
        // std::cout << "Pol vectors: \n";
        // std::cout << phononDispersion.at(idx).polVectors[0][0] << " " << phononDispersion.at(idx).polVectors[0][1] << " " << phononDispersion.at(idx).polVectors[0][2] << std::endl;
        // std::cout << phononDispersion.at(idx).polVectors[1][0] << " " << phononDispersion.at(idx).polVectors[1][1] << " " << phononDispersion.at(idx).polVectors[1][2] << std::endl;
        // std::cout << phononDispersion.at(idx).polVectors[2][0] << " " << phononDispersion.at(idx).polVectors[2][1] << " " << phononDispersion.at(idx).polVectors[2][2] << std::endl;
        // std::cout << "C vector: \n";
        // std::cout << diag.C.at(0) << " " << diag.C.at(1) << " " << diag.C.at(2) << std::endl;
        // std::cout << "D vector: \n";
        // std::cout << diag.D.at(0) << " " << diag.D.at(1) << " " << diag.D.at(2) << std::endl;

        // Construct the strings for each output based on the computation
        std::ostringstream evStream, cdStream, eVecStream, angMomStream;

        for (int i = 0; i < 7; i++)
        {
            evStream << diag.eigenEnergies.at(i) << ",";
        }
        evStream << diag.eigenEnergies.at(7);

        cdStream << diag.C.at(0) << "," << diag.C.at(1) << "," << diag.C.at(2) << ","
                 << diag.D.at(0) << "," << diag.D.at(1) << "," << diag.D.at(2);

        for (int col = 0; col < 8; col++)
        {
            for (int row = 0; row < 8; row++)
            {
                eVecStream << diag.eigenvectors_inverse.col(col).row(row).x();
                if (row < 7 || col < 7)
                {
                    eVecStream << ",";
                }
            }
        }


        for (int nu = 0 ; nu < 8 ; nu++)
        {
            angMomStream << diag.angularMomentum.at(nu) << ",";
        }
        // Assign the constructed strings to the corresponding vectors
        outEV[idx] = evStream.str();
        outCD[idx] = cdStream.str();
        outEVectors[idx] = eVecStream.str();
        outAngularMomentum[idx] = angMomStream.str();
    }

    std::ofstream outFileEV("Outputs/gH_1000/Eigenenergies.txt");
    std::ofstream outFileCD("Outputs/gH_1000/CD.txt");
    std::ofstream outFileEVectors("Outputs/gH_1000/EVec.txt");
    std::ofstream outFileAngularMomentum("Outputs/gH_1000/AngularMomentum.txt");

    for (const auto &line : outEV)
    {
        outFileEV << line << "\n";
    }
    for (const auto &line : outCD)
    {
        outFileCD << line << "\n";
    }
    for (const auto &line : outEVectors)
    {
        outFileEVectors << line << "\n";
    }
    for (const auto &line : outAngularMomentum)
    {
        outFileAngularMomentum << line << "\n";
    }

    outFileAngularMomentum.close();
    outFileEVectors.close();
    outFileEV.close();
    outFileCD.close();

    return 0;
    //*/

    // init irreducible BZ
    IrreducibleBZ irrBZ;
    irrBZ.init("Parameters/8x8x8/irrPoints_edit.txt");
    irrBZ.initMagnonPhononDispFromFile("Parameters/8x8x8/ph_data_edit.txt", "Parameters/J_bccFe.txt", "Parameters/8x8x8/magDisp.txt");
    irrBZ.readMultiplicities("Parameters/8x8x8/multiplicity_edit.txt");
    irrBZ.initSymmOp("Parameters/tu_graz_symm_bcc.txt");

    // irrBZ.init_k_prime();
    // irrBZ.save_k_primes_to_file("Outputs/k_primes.txt");
    irrBZ.init_k_primes_from_file("Outputs/k_primes.txt");

    // irrBZ.initCoefficients(parameters, nFT);
    // irrBZ.saveCoefficientsAsSqrtAbs("Outputs/coefficients.txt");
    irrBZ.readCoefficients("Outputs/coefficients.txt");
    irrBZ.initOccNumbers(2, 5); // 25.85, 30.0 (room temp) 10, 15
    irrBZ.integrate();

    /*
    // Testing calculation of k_primes
    int signs[3] = {1, 1, 1};

    for (auto k : irrBZ.irreducibleBZVectors)
    {
        for (auto q : irrBZ.irreducibleBZVectors)
        {
            int n = irrBZ.get_k_prime(k, q, signs);
            if (n == -1)
            {
                std::cout << "k: " << k(0) / (2 * pi) << " " << k(1) / (2 * pi) << " " << k(2) / (2 * pi) << " | q: " << q(0) / (2 * pi) << " " << q(1) / (2 * pi) << " " << q(2) / (2 * pi) << std::endl;
            }
            else
            {
                std::cout << "Success " << n << std::endl;
            }
        }
    }
    */

    // Testing the symmetry operations
    // Eigen::Vector3d q(1, 2, 3);
    // for (auto op: irrBZ.operations)
    //{
    //    Eigen::Vector3d transf = op(q);
    //    std::cout << "Symmetry operation: " << transf(0) << " " << transf(1) << " " << transf(2) << std::endl;
    //}

    // for (auto k : irrBZ.irreducibleBZVectors)
    //{
    //     for (auto q : irrBZ.irreducibleBZVectors)
    //     {
    //         auto sum = k + q;
    //
    //         int rep = irrBZ.findRepresentative(sum);
    //         if (rep == -1)
    //         {
    //             std::cout << "Error: Could not find representative" << std::endl;
    //             std::cout << "k: " << k(0) << " " << k(1) << " " << k(2) << " | q: " << q(0) << " " << q(1) << " " << q(2) << " | sum: " << sum(0) << " " << sum(1) << " " << sum(2) << std::endl;
    //         }
    //         else
    //         {
    //
    //             std::cout << "Success " << rep << std::endl;
    //         }
    //         std::cout << k(0) / (2 * pi) << " " << k(1) / (2 * pi) << " " << k(2) / (2 * pi) << " | " << q(0) / (2 * pi) << " " << q(1) / (2 * pi) << " " << q(2) / (2 * pi) << std::endl;
    //         std::cout << "Sum: " << sum(0) / (2 * pi) << " " << sum(1) / (2 * pi) << " " << sum(2) / (2 * pi) << " | rep: " << irrBZ.irreducibleBZVectors.at(rep)(0) / (2 * pi) << " " << irrBZ.irreducibleBZVectors.at(rep)(1) / (2 * pi) << " " << irrBZ.irreducibleBZVectors.at(rep)(2) / (2 * pi) << std::endl;
    //     }
    // }
    return 0;

    // Check if the irr points are their own representative
    // for (int i = 0; i < irrBZ.irreducibleBZVectors.size(); i++)
    //{
    //    Eigen::Vector3d vec = irrBZ.irreducibleBZVectors.at(i);
    //    int representative = irrBZ.findRepresentative(vec);
    //    if (representative != i)
    //    {
    //        std::cout << "Test failed" << std::endl;
    //    }
    //    else {
    //        std::cout << "Test passed" << std::endl;
    //    }
    //}

    // print symmetry groups
    // for (auto mat : symGroups.symmetrieOperations)
    //{
    //    std::cout << "Symmetry matrix: " << mat << std::endl;
    //}

    /*
    std::ofstream outputFile("Outputs/symmApplied.txt");

    for (int j = 0; j < irrBZ.irreducibleBZVectors.size(); j++)
    {
        auto k = irrBZ.irreducibleBZVectors.at(j);
        std::vector<Eigen::Vector3d> symApplied = irrBZ.symmetryGroup.applySymmetries(k);
        for (int i = 0; i < symApplied.size(); i++)
        {
            auto k_prime = symApplied.at(i);
            // int representant_idx = irrBZ.findRepresentative(k_prime);
            outputFile << k_prime(0)/(2 * pi) << " " << k_prime(1)/(2 * pi) << " " << k_prime(2)/(2 * pi) << "\n"; // << " " << j << " " << representant_idx << "\n";
        }
    }
    outputFile.close();
    */

    for (auto k : irrBZ.irreducibleBZVectors)
    {
        std::cout << "k: " << k(0) / (2 * pi) << " " << k(1) / (2 * pi) << " " << k(2) / (2 * pi) << std::endl;
    }

    int succ_count = 0;
    int fail_count = 0;

    for (auto k : irrBZ.irreducibleBZVectors)
    {
        for (auto k_prime : irrBZ.irreducibleBZVectors)
        {
            auto q = k - k_prime;
            int representant_idx = irrBZ.findRepresentative(q);
            if (representant_idx == -1)
            {
                fail_count++;
                std::cout << "Error: Could not find representative" << std::endl;
                std::cout << "k: " << k(0) / (2 * pi) << " " << k(1) / (2 * pi) << " " << k(2) / (2 * pi) << " | k_prime: " << k_prime(0) / (2 * pi) << " " << k_prime(1) / (2 * pi) << " " << k_prime(2) / (2 * pi) << " | q: " << q(0) / (2 * pi) << " " << q(1) / (2 * pi) << " " << q(2) / (2 * pi) << std::endl;
            }
            else
            {
                succ_count++;
            }
            // std::cout << "k: " << k(0) << " " << k(1) << " " << k(2) << " | k_prime: " << k_prime(0) << " " << k_prime(1) << " " << k_prime(2) << " | q: " << q(0) << " " << q(1) << " " << q(2) << " | representant: " << representant_idx << std::endl;
        }
    }

    std::cout << "Success: " << succ_count << " | Fail: " << fail_count << std::endl;

    /*
    // Calculate the distance between all irreducible BZ vectors
    std::vector<double> distances;
    for (auto k : irrBZ.irreducibleBZVectors)
    {
        for (auto q : irrBZ.irreducibleBZVectors)
        {
            distances.push_back(distance(k, q));
        }
    }
    std::sort(distances.begin(), distances.end(), std::greater<double>());

    int counter = 0;
    for (auto d : distances)
    {
        counter++;
        std::cout << counter << "  " << d << std::endl;
    }
    return 0;
   */

    /*
    std::vector<double> multiplicities(irrBZ.irreducibleBZVectors.size(), 0.0);
    irrBZ.initMultiplicities();

    std::vector<CouplingParameter> dynmat = readDynMatrices("Parameters/dynMat_8x8x8.txt");

    for (int i = 0; i < dynmat.size(); i++)
    {
        auto D = dynmat.at(i);
        // if (i >= 1)
        //{
        //     continue;
        // }
        const Eigen::Vector3d vec(D.kx, D.ky, D.kz);

        int representative = irrBZ.findRepresentative(vec);
        multiplicities.at(representative) += 1;
        if (representative != -1)
        {
            // std::cout << "success" << std::endl;
            //  std::cout << "vec - repr: " << vec(0) - irrBZ.irreducibleBZVectors.at(representative)(0) << " " << vec(1) - irrBZ.irreducibleBZVectors.at(representative)(1) << " " << vec(2) - irrBZ.irreducibleBZVectors.at(representative)(2) << std::endl;
        }
    }

    for (int idx = 0; idx < multiplicities.size(); idx++)
    {
        std::cout << " | Multiplicity: " << multiplicities.at(idx) << " | vec: " << irrBZ.irreducibleBZVectors.at(idx)(0) << "," << irrBZ.irreducibleBZVectors.at(idx)(1) << "," << irrBZ.irreducibleBZVectors.at(idx)(2) << std::endl;
    }
    */

    /*
    for (auto vec : irrBZ.irreducibleBZVectors)
    {
        Eigen::Vector3d vecTemp = vec;
        int representative = irrBZ.findRepresentative(vecTemp);
        std::cout << "vec - repr: " << vec(0) - irrBZ.irreducibleBZVectors.at(representative)(0) << " " << vec(1) - irrBZ.irreducibleBZVectors.at(representative)(1) << " " << vec(2) - irrBZ.irreducibleBZVectors.at(representative)(2) << std::endl;
    }

    double sum = 0;

    for (int i = 0; i < irrBZ.irreducibleBZVectors.size(); i++)
    {
        std::cout << "vec: " << irrBZ.irreducibleBZVectors.at(i)(0) << "," << irrBZ.irreducibleBZVectors.at(i)(1) << "," << irrBZ.irreducibleBZVectors.at(i)(2) << " | Multiplicity: " << irrBZ.multiplicities.at(i) << std::endl;
        sum += irrBZ.multiplicities.at(i);
    }


    std::cout << sum << std::endl;
    */

    /*
    using namespace std::chrono;

    SamplingGrid samplingGrid;
    samplingGrid.initMonkhorstPack();
    samplingGrid.initMagnonDisp("Parameters/J_bccFe.txt");
    samplingGrid.initPhononDisp("Parameters/dynMat_16x16x16.txt", "Parameters/nn5.txt");
    auto start = high_resolution_clock::now();
    samplingGrid.init(parameters);
    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(stop - start);

    // Output the execution time
    std::cout << "Function execution took " << duration.count() << " milliseconds" << std::endl;

    return 0;
    */

    /*

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

    return 0;

    */

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
