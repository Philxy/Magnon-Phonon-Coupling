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

    // remains to be implemented

    return BZVector(0, 0, 0);
}
