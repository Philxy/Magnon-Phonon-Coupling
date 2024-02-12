#include "../include/fourierTransform.h"

std::complex<double> FTD(double kx, double ky, double kz, std::vector<CouplingParameter> params, Axis nu, Axis mu)
{

    std::complex<double> result(0, 0);

    const std::complex<double> i(0.0, 1.0);

    for (CouplingParameter p : params)
    {

        if (p.displacementDirection != mu)
        {
            continue;
        }

        double D = 0;

        switch (nu)
        {
        case X:
            D = -0.5 * (p.J[Z][Y] - p.J[Y][Z]);
            break;
        case Y:
            D = 0.5 * (p.J[X][Z] - p.J[Z][X]);
            break;
        case Z:
            D = 0.5 * (p.J[X][Y] - p.J[Y][X]);
            break;
        default:
            continue;
        }

        result += D * std::exp(+i * (kx * p.x_ij + ky * p.y_ij + kz * p.z_ij) + i * (kx * p.x_ki + ky * p.y_ki + kz * p.z_ki));
    }

    return result;
}

std::complex<double> FTD(double kx1, double ky1, double kz1, double kx2, double ky2, double kz2, std::vector<CouplingParameter> params, Axis nu, Axis mu)
{
    std::complex<double> result(0, 0);

    const std::complex<double> i(0.0, 1.0);

    for (CouplingParameter p : params)
    {

        if (p.displacementDirection != mu)
        {
            continue;
        }

        double D;

        switch (nu)
        {
        case X:
            D = 0.5 * (p.J[Z][Y] - p.J[Y][Z]);
            break;
        case Y:
            D = 0.5 * (p.J[X][Z] - p.J[Z][X]);
        case Z:
            D = 0.5 * (p.J[X][Y] - p.J[Y][X]);
            break;
        default:
            continue;
        }

        // set off-diagonal elements to zero
        if (mu == nu)
        {
            D = 0;
            continue;
        }

        result += D * std::exp(+i * (kx1 * p.x_ki + ky1 * p.y_ki + kz1 * p.z_ki) - i * (kx2 * p.x_jk + ky2 * p.y_jk + kz2 * p.z_jk));
    }

    return result;
}

std::complex<double> FTJiso(double kx, double ky, double kz, std::vector<CouplingParameter> params)
{
    std::complex<double> result;
    const std::complex<double> i(0.0, 1.0);

    for (CouplingParameter p : params)
    {
        result += (1.0 - std::exp(i * (p.x * kx + p.y * ky + p.z * kz))) * p.J_iso;
    }

    return result;
}
