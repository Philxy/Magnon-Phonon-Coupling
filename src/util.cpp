#include "../include/util.h"

void eigenMatrixToCArray(const Eigen::Matrix3d &eigenMatrix, double array[3][3])
{
    for (int i = 0; i < eigenMatrix.rows(); ++i)
    {
        for (int j = 0; j < eigenMatrix.cols(); ++j)
        {
            array[i][j] = eigenMatrix(i, j);
        }
    }
}

void sortEigen(Eigen::Vector3cd &eigenvalues, Eigen::Matrix3cd &eigenvectors)
{

    Eigen::Vector3cd sorted_eigenvalues;
    Eigen::Matrix3cd sorted_eigenvectors;

    for (int i = 0; i < 2; i++)
    {

        int max_index = -1;
        double max_element = -1;
        for (int j = 0; j < 2; j++)
        {
            double curr_element = std::abs(eigenvectors(j, i).real());
            // std::cout << curr_element << std::endl;
            if (curr_element > max_element)
            {
                max_element = curr_element;
                max_index = j;
            }
        }
        sorted_eigenvalues(max_index) = eigenvalues(i);
        sorted_eigenvectors(0, i) = eigenvectors(0, max_index);
        sorted_eigenvectors(1, i) = eigenvectors(1, max_index);
        sorted_eigenvectors(2, i) = eigenvectors(2, max_index);
    }

    eigenvalues = sorted_eigenvalues;
    eigenvectors = sorted_eigenvectors;
}

// Switches sign of eigenvalues if it is negative and adjust the corresponding eigenvectors too.
void makeEigenvaluesPositive(Eigen::Vector3cd &eigenvalues, Eigen::Matrix3cd &eigenvectors)
{
    for (int idx = 0; idx < 3; idx++)
    {
        if (eigenvalues(idx).real() < 0)
        {
            eigenvalues(idx) *= -1;
            eigenvectors.col(idx).x() *= -1;
            eigenvectors.col(idx).y() *= -1;
            eigenvectors.col(idx).z() *= -1;
        }
    }
    return;
}

void changeEigenVecSign(Eigen::Matrix3cd &eigenvectors)
{
    for (int idx = 0; idx < 3; idx++)
    {
        if (eigenvectors.col(idx).x().real() < 0)
        {
            eigenvectors.col(idx).x() *= -1;
            eigenvectors.col(idx).y() *= -1;
            eigenvectors.col(idx).z() *= -1;
        }
    }
    return;
}
