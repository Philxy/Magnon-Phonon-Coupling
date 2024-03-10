#include "path.h"
#include "dispersion.h"
#include "globals.h"


struct BZVector
{
    int n[3];              // tree integers uniquely defining the k vector within the first BZ
    double kx, ky, kz;     // k vector components
    Vector3D k;            // k vector
    double n_ph[3], n_mag; // phonon and magnon occupation number

    PhononDispParam ph;  // energy and pol vectors of the phonon
    MagnonDispParam mag; // energy of the magnon mode

    BZVector(int n1, int n2, int n3);
    BZVector(){};
};

BZVector getG_gammaM(BZVector k, BZVector q);
BZVector getG_gammaP(BZVector k, BZVector q);
BZVector getG_gammaZ(BZVector k, BZVector q);

struct SingleCoefficient
{
    std::complex<double> coefficient;
    BZVector G;
};

template <typename T>
struct FixedMultiDimGrid
{
    std::vector<T> grid;

    // Constructor automatically sizes the grid based on the fixed dimensions
    FixedMultiDimGrid()
    {
        grid.resize(N * N * N * N * N * N * 3);
    }

    // Provides access to the grid elements as if it were a 6D array followed by a triplet
    T &at(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t p)
    {
        size_t index = ((((((i * N + j) * N + k) * N + l) * N + m) * N + n) * 3 + p);
        return grid[index];
    }
};

struct SamplingGrid
{

    void initMonkhorstPack();

    void initMagnonDisp(std::string couplingParameterFile, double S);
    void initPhononDisp(std::string dynamicMatricesFile, std::string nextNeighbourFile, double atomic_mass);

    BZVector at(int n1, int n2, int n3);

    double S = 1.1;
    double ftN = 1; // welchen wert hat das???
    double atomic_mass = 55.56;

    void init(const std::vector<CouplingParameter> &parameters); // must be initialized!

private:
    BZVector grid[N][N][N];

    FixedMultiDimGrid<SingleCoefficient> gammaPlusGrid;
    FixedMultiDimGrid<SingleCoefficient> gammaMinusGrid;
    FixedMultiDimGrid<SingleCoefficient> gammaZGrid;

    // SingleCoefficient gammaPlusGrid[N][N][N][N][N][N][3];
    // SingleCoefficient gammaMinusGrid[N][N][N][N][N][N][3];
    // SingleCoefficient gammaZGrid[N][N][N][N][N][N][3];
};
