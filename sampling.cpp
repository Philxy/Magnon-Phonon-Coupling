#include "../include/sampling.h"

BZVector::BZVector(int n1, int n2, int n3)
{
    n[0] = n1;
    n[1] = n2;
    n[2] = n3;

    // n1,n2,n3 take values from 0 to N-1

    kx = (2.0 * n1 - N) / (2.0 * N) * b1.x + (2.0 * n2 - N) / (2.0 * N) * b2.x + (2.0 * n3 - N) / (2.0 * N) * b3.x;
    ky = (2.0 * n1 - N) / (2.0 * N) * b1.y + (2.0 * n2 - N) / (2.0 * N) * b2.y + (2.0 * n3 - N) / (2.0 * N) * b3.y;
    kz = (2.0 * n1 - N) / (2.0 * N) * b1.z + (2.0 * n2 - N) / (2.0 * N) * b2.z + (2.0 * n3 - N) / (2.0 * N) * b3.z;
    k = Vector3D(kx, ky, kz);
}

void SamplingGrid::initMonkhorstPack()
{
    for (int n1 = 0; n1 < N; n1++)
    {
        for (int n2 = 0; n2 < N; n2++)
        {
            for (int n3 = 0; n3 < N; n3++)
            {
                BZVector k(n1, n2, n3);
                grid[n1][n2][n3] = k;
            }
        }
    }
}

BZVector SamplingGrid::at(int n1, int n2, int n3)
{
    int corrected_n1 = (n1 >= N) ? n1 - N : n1;
    int corrected_n2 = (n2 >= N) ? n3 - N : n2;
    int corrected_n3 = (n3 >= N) ? n3 - N : n3;

    corrected_n1 = (n1 < 0) ? n1 + N : n1;
    corrected_n2 = (n2 < 0) ? n3 + N : n2;
    corrected_n3 = (n3 < 0) ? n3 + N : n3;

    return grid[corrected_n1][corrected_n2][corrected_n3];
}

void SamplingGrid::initMagnonDisp(std::string couplingParameterFile)
{

    std::vector<Vector3D> path;

    for (int n1 = 0; n1 < N; n1++)
    {
        for (int n2 = 0; n2 < N; n2++)
        {
            for (int n3 = 0; n3 < N; n3++)
            {
                const BZVector &k = grid[n1][n2][n3];
                path.push_back(Vector3D(k.kx, k.ky, k.kz));
            }
        }
    }

    std::vector<MagnonDispParam> magDisp = getMagneticDispersion(couplingParameterFile, path);

    int counter = 0;
    for (int n1 = 0; n1 < N; n1++)
    {
        for (int n2 = 0; n2 < N; n2++)
        {
            for (int n3 = 0; n3 < N; n3++)
            {
                grid[n1][n2][n3].mag = magDisp.at(counter);
                counter++;
            }
        }
    }
}

void SamplingGrid::initPhononDisp(std::string dynamicMatricesFile, std::string nextNeighbourFile)
{
    std::vector<Vector3D> path;

    for (int n1 = 0; n1 < N; n1++)
    {
        for (int n2 = 0; n2 < N; n2++)
        {
            for (int n3 = 0; n3 < N; n3++)
            {
                const BZVector &k = grid[n1][n2][n3];
                path.push_back(Vector3D(k.kx, k.ky, k.kz));
            }
        }
    }

    std::vector<PhononDispParam> phDisp = getPhononDispersion(dynamicMatricesFile, nextNeighbourFile, path);

    int counter = 0;
    for (int n1 = 0; n1 < N; n1++)
    {
        for (int n2 = 0; n2 < N; n2++)
        {
            for (int n3 = 0; n3 < N; n3++)
            {
                grid[n1][n2][n3].ph = phDisp.at(counter);
                counter++;
            }
        }
    }
}

void SamplingGrid::init(const std::vector<CouplingParameter> &parameters)
{

    std::complex<double> i(0, 1);

    int counter = 0;

    for (int n1 = 0; n1 < N; n1++)
    {
        for (int n2 = 0; n2 < N; n2++)
        {
            for (int n3 = 0; n3 < N; n3++)
            {
                for (int m1 = 0; m1 < N; m1++)
                {
                    for (int m2 = 0; m2 < N; m2++)
                    {
                        for (int m3 = 0; m3 < N; m3++)
                        {
                            BZVector k(n1, n2, n3);
                            BZVector q(m1, m2, m3);

                            for (int branch = 0; branch < 3; branch++)
                            {

                                SingleCoefficient gammaM, gammaP, gammaZ;

                                gammaM.G = getG_gammaM(k, q);
                                gammaP.G = getG_gammaP(k, q);
                                gammaZ.G = getG_gammaZ(k, q);

                                for (Axis ax : std::vector<Axis>{X, Y, Z})
                                {
                                    std::complex<double> gammaMJxxmu = J_kq(k.kx, k.ky, k.kz, q.kx, q.ky, q.kz, parameters, X, X, ax);
                                    std::complex<double> gammaMJyymu = J_kq(k.kx, k.ky, k.kz, q.kx, q.ky, q.kz, parameters, Y, Y, ax);
                                    std::complex<double> gammaMJyxmu = J_kq(k.kx, k.ky, k.kz, q.kx, q.ky, q.kz, parameters, Y, X, ax);
                                    std::complex<double> gammaMJxymu = J_kq(k.kx, k.ky, k.kz, q.kx, q.ky, q.kz, parameters, X, Y, ax);

                                    std::complex<double> gammaPJxxmu = J_kq(-k.kx, -k.ky, -k.kz, q.kx, q.ky, q.kz, parameters, X, X, ax);
                                    std::complex<double> gammaPJyymu = J_kq(-k.kx, -k.ky, -k.kz, q.kx, q.ky, q.kz, parameters, Y, Y, ax);
                                    std::complex<double> gammaPJyxmu = J_kq(-k.kx, -k.ky, -k.kz, q.kx, q.ky, q.kz, parameters, Y, X, ax);
                                    std::complex<double> gammaPJxymu = J_kq(-k.kx, -k.ky, -k.kz, q.kx, q.ky, q.kz, parameters, X, Y, ax);

                                    gammaM.coefficient += sqrt(ftN) / (2.0 * S) * (gammaMJxxmu - gammaMJyymu - i * gammaMJxymu - i * gammaMJyxmu) * sqrt(1 / (2 * atomicMass * grid[m1][m2][m3].ph.E[branch])) * grid[m1][m2][m3].ph.polVectors[ax][branch];
                                    gammaP.coefficient += sqrt(ftN) / (2.0 * S) * (gammaPJxxmu - gammaPJyymu + i * gammaPJxymu + i * gammaPJyxmu) * sqrt(1 / (2 * atomicMass * grid[m1][m2][m3].ph.E[branch])) * grid[m1][m2][m3].ph.polVectors[ax][branch];

                                    std::complex<double> Jxxmu = J_kq(-k.kx, -k.ky, -k.kz, q.kx, q.ky, q.kz, parameters, X, X, ax);
                                    std::complex<double> Jyymu = J_kq(-k.kx, -k.ky, -k.kz, q.kx, q.ky, q.kz, parameters, Y, Y, ax);
                                    std::complex<double> Jyxmu = J_kq(-k.kx, -k.ky, -k.kz, q.kx, q.ky, q.kz, parameters, Y, X, ax);
                                    std::complex<double> Jxymu = J_kq(-k.kx, -k.ky, -k.kz, q.kx, q.ky, q.kz, parameters, X, Y, ax);
                                    std::complex<double> Jzzmu = J_kq(gammaZ.G.kx - q.kx, gammaZ.G.ky - q.ky, gammaZ.G.kz - q.kz, q.kx, q.ky, q.kz, parameters, X, Y, ax);
                                    gammaZ.coefficient += sqrt(ftN) / (2.0 * S) * (2.0 * Jxxmu + 2.0 * Jyymu - 4.0 * Jzzmu + 2.0 * i * Jyxmu - 2.0 * i * Jxymu) * sqrt(1 / (2 * atomicMass * grid[m1][m2][m3].ph.E[branch])) * grid[m1][m2][m3].ph.polVectors[ax][branch];
                                }

                                gammaPlusGrid.at(n1, n2, n3, m1, m2, m3, branch) = gammaP;
                                gammaMinusGrid.at(n1, n2, n3, m1, m2, m3, branch) = gammaM;
                                gammaZGrid.at(n1, n2, n3, m1, m2, m3, branch) = gammaZ;

                                counter++;

                                if (counter % 10)
                                {
                                    std::cout << "init progress: " << (double)counter / (N*N*N*N*N*N*3) << std::endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// Given two vectors in the first BZ, this function returns the only lattice vector G satisfying G = k + k^prime + q,
// such that k^\prime is in the first BZ as well.
BZVector getG_gammaM(BZVector k, BZVector q)
{
    int g1, g2, g3;

    for (int m : {-1, 0, 1})
    {
        int nkp1 = (2 * N * m + 3 * N) / 2 - k.n[0] - q.n[0];
        int nkp2 = (2 * N * m + 3 * N) / 2 - k.n[1] - q.n[1];
        int nkp3 = (2 * N * m + 3 * N) / 2 - k.n[2] - q.n[2];

        if (nkp1 >= 0 && nkp1 <= N - 1)
        {
            g1 = nkp1;
        }
        if (nkp2 >= 0 && nkp2 <= N - 1)
        {
            g2 = nkp2;
        }
        if (nkp3 >= 0 && nkp3 <= N - 1)
        {
            g3 = nkp3;
        }
    }

    BZVector G(g1, g2, g3);

    return G;
}

BZVector getG_gammaP(BZVector k, BZVector q)
{
    int g1, g2, g3;

    for (int m : {-1, 0, 1})
    {
        int nkp1 = -((2 * N * m + 3 * N) / 2 + k.n[0] - q.n[0]);
        int nkp2 = -((2 * N * m + 3 * N) / 2 + k.n[1] - q.n[1]);
        int nkp3 = -((2 * N * m + 3 * N) / 2 + k.n[2] - q.n[2]);

        if (nkp1 >= 0 && nkp1 <= N - 1)
        {
            g1 = nkp1;
        }
        if (nkp2 >= 0 && nkp2 <= N - 1)
        {
            g2 = nkp2;
        }
        if (nkp3 >= 0 && nkp3 <= N - 1)
        {
            g3 = nkp3;
        }
    }

    BZVector G(g1, g2, g3);

    return G;
}

BZVector getG_gammaZ(BZVector k, BZVector q)
{
    int g1, g2, g3;

    for (int m : {-1, 0, 1})
    {
        int nkp1 = ((2 * N * m + 3 * N) / 2 + k.n[0] - q.n[0]);
        int nkp2 = ((2 * N * m + 3 * N) / 2 + k.n[1] - q.n[1]);
        int nkp3 = ((2 * N * m + 3 * N) / 2 + k.n[2] - q.n[2]);

        if (nkp1 >= 0 && nkp1 <= N - 1)
        {
            g1 = nkp1;
        }
        if (nkp2 >= 0 && nkp2 <= N - 1)
        {
            g2 = nkp2;
        }
        if (nkp3 >= 0 && nkp3 <= N - 1)
        {
            g3 = nkp3;
        }
    }

    BZVector G(g1, g2, g3);

    return G;
}