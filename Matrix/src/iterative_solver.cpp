#include "../include/iterative_solver.hpp"
#include <cmath>
#include <iostream>

using namespace std;

// ✅ Check Diagonal Dominance
bool IterativeSolver::isDiagonallyDominant(const Matrix &A)
{
    int n = A.getRows();

    for (int i = 0; i < n; i++)
    {
        double diag = abs(A(i, i));
        double sum = 0;

        for (int j = 0; j < n; j++)
        {
            if (i != j)
                sum += abs(A(i, j));
        }

        if (diag < sum)
            return false;
    }
    return true;
}

// ✅ Gauss-Jacobi
vector<double> IterativeSolver::gaussJacobi(const Matrix &A, const vector<double> &b, int maxIter, double tol)
{
    int n = A.getRows();
    vector<double> x(n, 0.0), x_new(n, 0.0);

    for (int k = 0; k < maxIter; k++)
    {
        for (int i = 0; i < n; i++)
        {
            double sum = 0;
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                    sum += A(i, j) * x[j];
            }
            x_new[i] = (b[i] - sum) / A(i, i);
        }

        // check convergence
        double err = 0;
        for (int i = 0; i < n; i++)
            err += abs(x_new[i] - x[i]);

        if (err < tol)
            break;

        x = x_new;
    }

    return x;
}

// ✅ Gauss-Seidel
vector<double> IterativeSolver::gaussSeidel(const Matrix &A, const vector<double> &b, int maxIter, double tol)
{
    int n = A.getRows();
    vector<double> x(n, 0.0);

    for (int k = 0; k < maxIter; k++)
    {
        vector<double> old = x;

        for (int i = 0; i < n; i++)
        {
            double sum = 0;
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                    sum += A(i, j) * x[j];
            }
            x[i] = (b[i] - sum) / A(i, i);
        }

        double err = 0;
        for (int i = 0; i < n; i++)
            err += abs(x[i] - old[i]);

        if (err < tol)
            break;
    }

    return x;
}