#include "../include/eigen.h"



Eigen::EigenSolver<Eigen::Matrix3d> getEigenValVec(Eigen::Matrix3d dynMat)
{
    Eigen::EigenSolver<Eigen::Matrix3d> solver(dynMat);

    Eigen::Vector3cd eigenvalues = solver.eigenvalues();
    Eigen::Matrix3cd eigenvectors = solver.eigenvectors();

    return solver;
}
