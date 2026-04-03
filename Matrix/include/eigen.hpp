#ifndef EIGEN_HPP
#define EIGEN_HPP

#include "matrix.hpp"
#include <vector>
#include <utility> // std::pair

// EigenAnalysis inherits Matrix.
// Stores the coefficient matrix (n x n) and provides
// a base interface for eigenvalue-related computations.
class EigenAnalysis : public Matrix
{
protected:
    int n; // order of the square matrix

public:
    // Construct from an existing n x n Matrix
    explicit EigenAnalysis(const Matrix &A);

    // Virtual destructor
    virtual ~EigenAnalysis() = default;

    // Returns the order of the matrix
    int getOrder() const;

    // Pure virtual: subclasses must implement their eigen-theorem computation
    virtual void compute() = 0;

    // Pure virtual: print results to stdout
    virtual void printResults() const = 0;
};

#endif
