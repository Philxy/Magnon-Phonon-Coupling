#include "path.h"
#include "dispersion.h"
#include "globals.h"

constexpr size_t N = 10; // size of the sampling grid (N,N,N)



struct BZVector
{
    int n[3];              // tree integers uniquely defining the k vector within the first BZ
    double kx, ky, kz;     // k vector components
    Vector3D k;            // k vector
    double n_ph[3], n_mag; // phonon and magnon occupation number

    PhononDispParam ph;  // energy and pol vectors of the phonon
    MagnonDispParam mag; // energy of the magnon mode

    BZVector(int n1, int n2, int n3);
};

struct SamplingGrid
{

    void initMonkhorstPack();

    BZVector at(int n1, int n2, int n3);

private:
    BZVector grid[N][N][N];
};
