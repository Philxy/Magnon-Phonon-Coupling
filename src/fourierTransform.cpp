#include "../include/fourierTransform.h"

std::complex<double> FTD(double kx, double ky, double kz, const std::vector<CouplingParameter> &params, Axis nu, Axis mu)
{

    std::complex<double> result(0, 0);

    const std::complex<double> i(0.0, 1.0);

    for (int idx = 0; idx < params.size(); idx++)
    {

        const CouplingParameter &p = params.at(idx);

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

// calculates D^{\nu\mu}_{\v{k}}
std::complex<double> DMIFT(double kx, double ky, double kz, const std::vector<CouplingParameter> &params, Axis nu, Axis mu)
{
    std::complex<double> result(0, 0);
    const std::complex<double> i(0.0, 1.0);

    for (int idx = 0; idx < params.size(); idx++)
    {
        const CouplingParameter &p = params.at(idx);

        if (p.displacementDirection != mu)
        {
            continue;
        }

        double D = 0;

        switch (nu)
        {
        case X:
            D = 0.5 * (p.J[Y][Z] - p.J[Z][Y]);
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
        result += D * std::exp(i * (kx * p.x_ij + ky * p.y_ij + kz * p.z_ij) + i * (kx * p.x_ki + ky * p.y_ki + kz * p.z_ki));
    }
    return result;
}

DMILike::DMILike(double kx, double ky, double kz, const std::vector<CouplingParameter> &parameters)
{
    this->kx = kx;
    this->ky = ky;
    this->kz = kz;

    std::vector<Axis> allAxis = {X, Y, Z};

    for (Axis nu : allAxis)
    {
        for (Axis mu : allAxis)
        {
            this->D[nu][mu] = DMIFT(kx, ky, kz, parameters, nu, mu);
        }
    }
}

DMILikeCouplingParam::DMILikeCouplingParam(double kx, double ky, double kz, const std::vector<CouplingParameter> &parameters)
{
    this->kx = kx;
    this->ky = ky;
    this->kz = kz;

    std::vector<Axis> allAxis = {X, Y, Z};

    for (Axis axis : allAxis)
    {
        this->D[axis][X] = FTD(kx, ky, kz, parameters, axis, X);
        this->D[axis][Y] = FTD(kx, ky, kz, parameters, axis, Y);
    }
}

/*
std::complex<double> FTD(double kx1, double ky1, double kz1, double kx2, double ky2, double kz2, const std::vector<CouplingParameter> &params, Axis nu, Axis mu)
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
*/

std::complex<double> FTJ(double kx, double ky, double kz, const std::vector<CouplingParameter> &params, Axis alpha, Axis beta, Axis mu)
{
    std::complex<double> result(0, 0);
    const std::complex<double> i(0.0, 1.0);

    for (CouplingParameter p : params)
    {
        if (mu !=  p.displacementDirection)
        {
            continue;  
        }

        result += p.J[alpha][beta] * std::exp(-i * (kx * p.x_ij + ky * p.y_ij + kz * p.z_ij)) * std::exp(-i * (kx * p.x_ki + ky * p.y_ki + kz * p.z_ki));
    }

    return result;
}

std::complex<double> FTJiso(double kx, double ky, double kz, const std::vector<CouplingParameter> &params)
{
    std::complex<double> result;
    const std::complex<double> i(0.0, 1.0);

    for (CouplingParameter p : params)
    {
        result += (1.0 - std::exp(i * (p.x * kx + p.y * ky + p.z * kz))) * p.J_iso;
    }

    return result;
}

Eigen::Matrix3d dynMat(double kx, double ky, double kz, const std::vector<CouplingParameter> &params)
{
    Eigen::Matrix3cd result;
    std::complex<double> i(0, 1);

    for (CouplingParameter p : params)
    {
        for (int alpha = 0; alpha <= 2; alpha++)
        {
            for (int beta = 0; beta <= 2; beta++)
            {
                result(alpha, beta) += p.Phi[alpha][beta] * std::exp(+i * (kx * p.x + ky * p.y + kz * p.z));
            }
        }
    }

    // Discard the imaginary part of the result and return only the real part
    Eigen::Matrix3d result_real;
    for (int alpha = 0; alpha <= 2; alpha++)
    {
        for (int beta = 0; beta <= 2; beta++)
        {
            result_real(alpha, beta) = result(alpha, beta).real();
        }
    }

    return result_real;
}

Eigen::Matrix3d forceMatrix(double x, double y, double z, const std::vector<CouplingParameter> &params)
{
    Eigen::Matrix3cd result;
    std::complex<double> i(0, 1);

    for (CouplingParameter p : params)
    {
        for (int alpha = 0; alpha <= 2; alpha++)
        {
            for (int beta = 0; beta <= 2; beta++)
            {
                result(alpha, beta) += p.D[alpha][beta] * std::exp(i * (x * p.kx + y * p.ky + z * p.kz));
            }
        }
    }

    // Discard the imaginary part of the result and return only the real part
    Eigen::Matrix3d result_real;
    for (int alpha = 0; alpha <= 2; alpha++)
    {
        for (int beta = 0; beta <= 2; beta++)
        {
            result_real(alpha, beta) = result(alpha, beta).real();
        }
    }

    return result_real;
}

std::complex<double> J_kq(double kx, double ky, double kz, double qx, double qy, double qz, const std::vector<CouplingParameter> &params, Axis alpha, Axis beta, Axis mu)
{
    std::complex<double> result(0, 0);
    const std::complex<double> i(0.0, 1.0);

    for (int idx = 0; idx < params.size(); idx++)
    {

        const CouplingParameter &p = params.at(idx);

        if (p.displacementDirection != mu)
        {
            continue;
        }

        result += p.J[alpha][beta] * std::exp(-i * (kx * p.x_ij + ky * p.y_ij + kz * p.z_ij) - i * (qx * p.x_jk + qy * p.y_jk + qz * p.z_jk));
    }
    return result;
}
