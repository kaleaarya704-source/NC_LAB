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
    LUSolver(int n, Method m);//Constructor initializes base class and stores selected method.

    void solve() override; //overrides virtual function solve from base class

private:
    void croutDecomposition();
    void doolittleDecomposition();
    void choleskyDecomposition();

    void forwardSubstitution(std::vector<std::vector<double>> &L,
                             std::vector<double> &y);//solves Ly = b
    void backwardSubstitution(std::vector<std::vector<double>> &U,
                              std::vector<double> &y);//solves Ux = y
    void pivotRows(int k);//performs partial pivoting on column k
};

#endif