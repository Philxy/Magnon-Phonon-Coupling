#include "../include/path.h"

// Function to linearly interpolate between two points
Vector3D lerp(const Vector3D &A, const Vector3D &B, double t)
{
    return {
        A.x + (B.x - A.x) * t,
        A.y + (B.y - A.y) * t,
        A.z + (B.z - A.z) * t};
}

// Function to generate N points on the line segment connecting A and B
std::vector<Vector3D> generatePointsOnLine(const Vector3D &A, const Vector3D &B, int N)
{
    std::vector<Vector3D> points;
    if (N <= 1)
    {
        // If N is 1 or less, only add the starting point to the list
        points.push_back(A);
    }
    else
    {
        for (int i = 0; i < N; ++i)
        {
            double t = i / (N - 1.0);
            points.push_back(lerp(A, B, t));
        }
    }
    return points;
}

// Constructs a path in k space along high symmetry lines (for bcc currently)
std::vector<Vector3D> constructPath(int number_of_k_vectors, double lattice_constant)
{
    std::vector<Vector3D> kVectors;

    double a = lattice_constant;
    const double PI = 3.14159265358979323846;

    // Construct the BCC path: Gamma H N Gamma P H | P N

    // Gamma - H --> (0,0,0) - (0,0,2pi/a)
    std::vector<Vector3D> path = generatePointsOnLine(Vector3D(0, 0, 0), Vector3D(0, 0, 2 * PI / a), number_of_k_vectors);
    kVectors.insert(kVectors.end(), path.begin(), path.end());

    return kVectors;

    // H - N --> (0,0,2PI/a) - (0,PI/a,PI/a)
    path = generatePointsOnLine(Vector3D(0, 0, 2 * PI / a), Vector3D(0, PI / a, PI / a), number_of_k_vectors);
    kVectors.insert(kVectors.end(), path.begin(), path.end());

    // N - Gamma --> (0,PI/a,PI/a) - (0,0,0)
    path = generatePointsOnLine(Vector3D(0, PI / a, PI / a), Vector3D(0, 0, 0), number_of_k_vectors);
    kVectors.insert(kVectors.end(), path.begin(), path.end());

    // Gamma - P --> (0,0,0) - (π/a,π/a,π/a)
    path = generatePointsOnLine(Vector3D(0, 0, 0), Vector3D(PI / a, PI / a, PI / a), number_of_k_vectors);
    kVectors.insert(kVectors.end(), path.begin(), path.end());

    // P - H --> (π/a,π/a,π/a) - (0,0,2PI/a)
    path = generatePointsOnLine(Vector3D(PI / a, PI / a, PI / a), Vector3D(0, 0, 2 * PI / a), number_of_k_vectors);
    kVectors.insert(kVectors.end(), path.begin(), path.end());

    // P - N --> (π/a,π/a,π/a) - (0,PI/a,PI/a)
    path = generatePointsOnLine(Vector3D(PI / a, PI / a, PI / a), Vector3D(0, PI / a, PI / a), number_of_k_vectors);
    kVectors.insert(kVectors.end(), path.begin(), path.end());

    return kVectors;
}

/*

// Generates a uniform grid in the Brillouin zone using Monkhorst–Pack sampling
// (for a bcc lattice atm. Change reciprocal lattice vectors to do this for other lattices!?)
std::vector<Vector3D> sampleBZ(int gridSize)
{
    double a = 1;
    const double PI = 3.14159265358979323846;

    std::vector<Vector3D> grid;

    Vector3D b1(2 * PI / a, 2 * PI / a, 0);
    Vector3D b2(0, 2 * PI / a, 2 * PI / a);
    Vector3D b3(2 * PI / a, 0, 2 * PI / a);

    for (int i = 1; i <= gridSize; i++)
    {
        for (int j = 1; j <= gridSize; j++)
        {
            for (int m = 1; m <= gridSize; m++)
            {
                double kx = (2 * i - gridSize - 1) / (2.0 * gridSize) * b1.x + (2 * j - gridSize - 1) / (2.0 * gridSize) * b2.x + (2 * m - gridSize - 1) / (2.0 * gridSize) * b3.x;
                double ky = (2 * i - gridSize - 1) / (2.0 * gridSize) * b1.y + (2 * j - gridSize - 1) / (2.0 * gridSize) * b2.y + (2 * m - gridSize - 1) / (2.0 * gridSize) * b3.y;
                double kz = (2 * i - gridSize - 1) / (2.0 * gridSize) * b1.z + (2 * j - gridSize - 1) / (2.0 * gridSize) * b2.z + (2 * m - gridSize - 1) / (2.0 * gridSize) * b3.z;
                grid.push_back(Vector3D(kx, ky, kz));
            }
        }
    }
    return grid;
}

*/
