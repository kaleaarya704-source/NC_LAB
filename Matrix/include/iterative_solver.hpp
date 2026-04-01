#ifndef ITERATIVE_SOLVER_HPP
#define ITERATIVE_SOLVER_HPP

#include "matrix.hpp"
#include <vector>

class IterativeSolver
{
public:
    static bool isDiagonallyDominant(const Matrix &A);

    static std::vector<double> gaussJacobi(const Matrix &A, const std::vector<double> &b, int maxIter, double tol);

    static std::vector<double> gaussSeidel(const Matrix &A, const std::vector<double> &b, int maxIter, double tol);
};

#endif