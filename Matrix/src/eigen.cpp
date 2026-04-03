#include "../include/eigen.hpp"
#include <stdexcept>

using namespace std;

// Construct from an existing n x n Matrix.
// Copies all elements into the inherited Matrix storage.
EigenAnalysis::EigenAnalysis(const Matrix &A)
    : Matrix(A), n(A.getRows())
{
    if (!A.isSquare())
        throw invalid_argument("EigenAnalysis requires a square matrix.");
}

int EigenAnalysis::getOrder() const
{
    return n;
}
