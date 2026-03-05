#ifndef SYSTEM_LINEAR_EQ_HPP
#define SYSTEM_LINEAR_EQ_HPP

#include "matrix.hpp"
#include <vector>

class SystemOfLinearEquation : public Matrix //This is class inherited from Matrix.
{
protected:
    std::vector<double> solution;

public:
    SystemOfLinearEquation(int n);

    virtual void solve() = 0;  // PURE VIRTUAL (Abstraction)
    std::vector<double> getSolution() const;
};

#endif