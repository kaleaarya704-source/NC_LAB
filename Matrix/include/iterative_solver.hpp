#ifndef ITERATIVE_SOLVER_HPP
#define ITERATIVE_SOLVER_HPP

#include "matrix.hpp"
#include <vector>

class IterativeSolver
{
public:
    static bool isDiagonallyDominant(const Matrix &A);

    // Attempts to make A (and corresponding b) diagonally dominant via row swapping.
    // Returns true if successful, false if no valid permutation exists.
    static bool makeDiagonallyDominant(Matrix &A, std::vector<double> &b);

    // omega = 1.0 -> standard Jacobi; omega < 1.0 -> under-relaxation (SOR) for non-DD matrices
    static std::vector<double> gaussJacobi(const Matrix &A, const std::vector<double> &b, int maxIter, double tol, double omega = 1.0);

    static std::vector<double> gaussSeidel(const Matrix &A, const std::vector<double> &b, int maxIter, double tol);
};

#endif