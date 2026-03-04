#ifndef LU_SOLVER_HPP
#define LU_SOLVER_HPP

#include "system_linear_eq.hpp"

class LUSolver : public SystemOfLinearEquation
{
public:
    enum Method
    {
        CROUT,
        DOOLITTLE,
        CHOLESKY
    };

private:
    Method method;

public:
    LUSolver(int n, Method m);

    void solve() override;

private:
    void croutDecomposition();
    void doolittleDecomposition();
    void choleskyDecomposition();

void forwardSubstitution(std::vector<std::vector<double>> &L,
                         std::vector<double> &y);
void backwardSubstitution(std::vector<std::vector<double>> &U,
                          std::vector<double> &y);
};

#endif