#include "../include/dynamics.h"

void Coefficients::initReciprocalLatticeVectors()
{

    for (int n1 = 0; n1 <= 1; n1++)
    {
        for (int n2 = 0; n2 <= 1; n2++)
        {
            for (int n3 = 0; n3 <= 1; n3++)
            {
                reciprocalLatticeVectors.push_back(Vector3D(n1 * b1.x + n2 * b2.x + n3 * b3.x, n1 * b1.y + n2 * b2.y + n3 * b3.y, n1 * b1.z + n2 * b2.z + n3 * b3.z));
            }
        }
    }
}

void Coefficients::init(const std::vector<CouplingParameter> &parameters)
{

    std::complex<double> i(0, 1);

    for (int branch = 0; branch < 3; branch++)
    {
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

                                    gammaMinus[n1][n2][n3][m1][m2][m3][branch] = sqrt(ftN) / (2.0 * S) * (gammaMJxxmu - gammaMJyymu - i * gammaMJxymu - i * gammaMJyxmu) * sqrt(1 / (2 * atomic_mass * samplingGrid.at(m1, m2, m3).ph.E[branch])) * samplingGrid.at(m1, m2, m3).ph.polVectors[ax][branch];
                                    gammaPlus[n1][n2][n3][m1][m2][m3][branch] = sqrt(ftN) / (2.0 * S) * (gammaPJxxmu - gammaPJyymu + i * gammaPJxymu + i * gammaPJyxmu) * sqrt(1 / (2 * atomic_mass * samplingGrid.at(m1, m2, m3).ph.E[branch])) * samplingGrid.at(m1, m2, m3).ph.polVectors[ax][branch];
                                }

                                for (int idxG = 0; idxG < reciprocalLatticeVectors.size(); idxG++)
                                {

                                    Vector3D G = reciprocalLatticeVectors.at(idxG);
                                    for (Axis ax : std::vector<Axis>{X, Y, Z})
                                    {
                                        std::complex<double> Jxxmu = J_kq(-k.kx, -k.ky, -k.kz, q.kx, q.ky, q.kz, parameters, X, X, ax);
                                        std::complex<double> Jyymu = J_kq(-k.kx, -k.ky, -k.kz, q.kx, q.ky, q.kz, parameters, Y, Y, ax);
                                        std::complex<double> Jyxmu = J_kq(-k.kx, -k.ky, -k.kz, q.kx, q.ky, q.kz, parameters, Y, X, ax);
                                        std::complex<double> Jxymu = J_kq(-k.kx, -k.ky, -k.kz, q.kx, q.ky, q.kz, parameters, X, Y, ax);
                                        std::complex<double> Jzzmu = J_kq(G.x - q.kx, G.y - q.ky, G.z - q.kz, q.kx, q.ky, q.kz, parameters, X, Y, ax);
                                        gammaZValues.at(idxG)[n1][n2][n3][m1][m2][m3][branch] = sqrt(ftN) / (2.0 * S) * (2.0 * Jxxmu + 2.0 * Jyymu - 4.0 * Jzzmu + 2.0 * i * Jyxmu - 2.0 * i * Jxymu) * sqrt(1 / (2 * atomic_mass * samplingGrid.at(m1, m2, m3).ph.E[branch])) * samplingGrid.at(m1, m2, m3).ph.polVectors[ax][branch];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}



