#pragma once
#include <vector>
#include "util.h"


struct Vector3D
{
    double x, y, z;

    // Constructor to initialize the vector
    Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}
};

// Function to linearly interpolate between two points
Vector3D lerp(const Vector3D &A, const Vector3D &B, double t);

// Function to generate N points on the line segment connecting A and B
std::vector<Vector3D> generatePointsOnLine(const Vector3D &A, const Vector3D &B, int N);

std::vector<Vector3D> constructPath(int number_of_k_vectors, double lattice_constant); 


std::vector<Vector3D> sampleBZ(int gridSize);