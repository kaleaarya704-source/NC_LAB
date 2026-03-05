#include "system_linear_eq.hpp"

SystemOfLinearEquation::SystemOfLinearEquation(int n)
    : Matrix(n, n + 1)//Initialize as augmented matrix (n x n+1)
{
}

std::vector<double> SystemOfLinearEquation::getSolution() const//Returns the computed solution vector.
{
    return solution;
}