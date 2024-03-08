#pragma once
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <iostream>

// Kartesian coordinate axes
enum Axis
{
    X,
    Y,
    Z
};

void eigenMatrixToCArray(const Eigen::Matrix3d &eigenMatrix, double array[3][3]);

void sortEigen(Eigen::Vector3cd &eigenvalues, Eigen::Matrix3cd &eigenvectors);

void makeEigenvaluesPositive(Eigen::Vector3cd &eigenvalues, Eigen::Matrix3cd &eigenvectors);

void changeEigenVecSign(Eigen::Matrix3cd &eigenvectors);
