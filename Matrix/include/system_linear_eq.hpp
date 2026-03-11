#ifndef SYSTEM_LINEAR_EQ_HPP
#define SYSTEM_LINEAR_EQ_HPP

#include "matrix.hpp"
#include <vector>

class SystemOfLinearEquation : public Matrix //This is class inherited from Matrix.
{
protected:
    std::vector<double> solution;//This will store the solution vector after solving the system.

public:
    SystemOfLinearEquation(int n);//Constructor initializes the base class as an augmented matrix (n x n+1)
    Matrix &getMatrix();
    virtual void solve() = 0;  // PURE VIRTUAL (Abstraction)
    std::vector<double> getSolution() const;//Returns the computed solution vector.
};

#endif
