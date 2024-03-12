#include "../include/dynamics.h"

// returns true or false depending on whether the k vector is inside the first BZ
// CHECK FOR CORRECTNESS REQUIRED
bool insideFirstBZ(Eigen::Vector3d kVec)
{
    double kx = kVec(0);
    double ky = kVec(1);
    double kz = kVec(2);

    if (kx < b_1(0) / 2 && kx > -b_1(0) / 2 && ky < b_2(1) / 2 && ky > -b_2(1) / 2 && kz < b_3(2) / 2 && kz > -b_3(2) / 2)
    {
        return true;
    }
    else
    {
        return false;
    }
}

Eigen::Vector3d mapToFirstBZ(Eigen::Vector3d kVec)
{

    for (int i = -2; i <= 2; i++)
    {
        for (int j = -2; j <= 2; j++)
        {
            for (int k = -2; k <= 2; k++)
            {
                Eigen::Vector3d kVec_temp = kVec + i * b_1 + j * b_2 + k * b_3;
                if (insideFirstBZ(kVec_temp))
                {
                    return kVec_temp;
                }
            }
        }
    }

    std::cout << "Error: Could not map k vector to first BZ" << std::endl;

    return Eigen::Vector3d(0, 0, 0);
}

double distance(Eigen::Vector3d k1, Eigen::Vector3d k2)
{
    Eigen::Vector3d diff = k1 - k2;
    return diff.norm();
}

std::vector<Eigen::Vector3d> SymmetrieGroup::applySymmetries(Eigen::Vector3d kVec)
{
    std::vector<Eigen::Vector3d> kVecs;
    for (Eigen::Matrix3d symOp : symmetrieOperations)
    {
        kVecs.push_back(symOp * kVec);
    }
    return kVecs;
}

int IrreducibleBZ::findRepresentative(Eigen::Vector3d kVec)
{
    kVec = mapToFirstBZ(kVec); // map k vector to first BZ (here or in the loop?)

    double dist_threshold = 0.001;

    for (Eigen::Vector3d kVec_temp : symmetrieGroup.applySymmetries(kVec))
    {
        for (int i = 0; i < irreducibleBZVectors.size(); i++)
        {
            double dist = distance(kVec, kVec_temp);
            if (dist < dist_threshold)
            {
                return i;
            }
        }
    }

    std::cout << "Error: Could not find representative in the irreducible BZ" << std::endl;
    return -1;
}

void IrreducibleBZ::init(std::string irreducibleBZFile)
{
    const double PI = 3.14159265358979323846;

    std::ifstream file(irreducibleBZFile);
    if (!file.is_open())
    {
        throw std::runtime_error("Error: Could not open file " + irreducibleBZFile);
    }

    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream ss(line);
        double x, y, z;
        char comma;

        if (!(ss >> x >> comma >> y >> comma >> z))
        {
            throw std::runtime_error("Error: Could not parse line " + line);
        }

        Eigen::Vector3d kVec(2 * PI * x, 2 * PI * y, 2 * PI * z);
        irreducibleBZVectors.push_back(kVec);
    }

    file.close();
}

void SymmetrieGroup::init(std::string symmetryFile)
{
    std::ifstream file(symmetryFile);
    if (!file.is_open())
    {
        throw std::runtime_error("Error: Could not open file " + symmetryFile);
    }

    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream ss(line);
        Eigen::Matrix3d matrix;

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                if (!(ss >> matrix(i, j)))
                {
                    throw std::runtime_error("Error: Could not parse line " + line);
                }
            }
        }

        symmetrieOperations.push_back(matrix);
    }

    // add the negative of the symmetrie matrices as well
    std::vector<Eigen::Matrix3d> symmetrieOperationsMinus;
    for (Eigen::Matrix3d symOp : symmetrieOperations)
    {
        symmetrieOperationsMinus.push_back(-symOp);
    }

    symmetrieOperations.insert(symmetrieOperations.end(), symmetrieOperationsMinus.begin(), symmetrieOperationsMinus.end());

    file.close();
}

void IrreducibleBZ::initMagnonDisp(std::string couplingParameterFile)
{

    std::vector<Vector3D> path;

    for (Eigen::Vector3d vec : irreducibleBZVectors)
    {
        path.push_back(Vector3D(vec(0), vec(1), vec(2)));
    }

    std::vector<MagnonDispParam> magDisp = getMagneticDispersion(couplingParameterFile, path);

    this->magDisp = magDisp;
}

void IrreducibleBZ::initPhononDisp(std::string dynamicMatricesFile, std::string nextNeighbourFile)
{

    std::vector<Vector3D> path;

    for (Eigen::Vector3d vec : irreducibleBZVectors)
    {
        path.push_back(Vector3D(vec(0), vec(1), vec(2)));
    }

    std::vector<PhononDispParam> phDisp = getPhononDispersion(dynamicMatricesFile, nextNeighbourFile, path);

    this->phDisp = phDisp;
}
