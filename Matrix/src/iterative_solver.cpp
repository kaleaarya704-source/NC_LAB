#include "../include/iterative_solver.hpp"
#include <cmath>
#include <iostream>

using namespace std;

//Check Diagonal Dominance
bool IterativeSolver::isDiagonallyDominant(const Matrix &A)
{
    int n = A.getRows();//Get size of matrix

    for (int i = 0; i < n; i++)//loop through each loop
    {
        double diag = abs(A(i, i));//diagonal element
        double sum = 0;//store sum of non-diagonal elements

        for (int j = 0; j < n; j++)//add values of non-diagonal elements
        {
            if (i != j)
                sum += abs(A(i, j));
        }

        if (diag < sum)//if sum of non-diagonal elements is less than diagonal element then matrix is not diagonally dominant
            return false;
    }
    return true;
}

//Gauss-Jacobi
vector<double> IterativeSolver::gaussJacobi(const Matrix &A, const vector<double> &b, int maxIter, double tol)
{
    int n = A.getRows();
    vector<double> x(n, 0.0), x_new(n, 0.0);

    for (int k = 0; k < maxIter; k++)//loop until max iterations
    {
        for (int i = 0; i < n; i++)//for each equation
        {
            double sum = 0;
            for (int j = 0; j < n; j++)//uses old x values
            {
                if (i != j)
                    sum += A(i, j) * x[j];
            }
            x_new[i] = (b[i] - sum) / A(i, i);//jacobi formula
        }

        // check convergence
        double err = 0;
        for (int i = 0; i < n; i++)
            err += abs(x_new[i] - x[i]);

        if (err < tol)//stop early if converged
            break;

        x = x_new;//update for next iteration
    }

    return x;
}

//Gauss-Seidel
vector<double> IterativeSolver::gaussSeidel(const Matrix &A, const vector<double> &b, int maxIter, double tol)
{
    int n = A.getRows();
    vector<double> x(n, 0.0);

    for (int k = 0; k < maxIter; k++)
    {
        vector<double> old = x;//stores previous iteration

        for (int i = 0; i < n; i++)
        {
            double sum = 0;
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                    sum += A(i, j) * x[j];//uses updated value
            }
            x[i] = (b[i] - sum) / A(i, i);//update immediately
        }

        double err = 0;
        for (int i = 0; i < n; i++)
            err += abs(x[i] - old[i]);//compare with previous iteration

        if (err < tol)
            break;
    }

    return x;
}