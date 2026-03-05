#include "lu_solver.hpp"
#include <cmath>
#include <stdexcept>

LUSolver::LUSolver(int n, Method m)
    : SystemOfLinearEquation(n), method(m) //Calls parent constructor,Stores selected method.
{
}

void LUSolver::solve()
{
    if (method == CROUT)
        croutDecomposition();
    else if (method == DOOLITTLE)
        doolittleDecomposition();
    else
        choleskyDecomposition();
}

void LUSolver::forwardSubstitution(std::vector<std::vector<double>> &L,
                                   std::vector<double> &y)//solve this system L y = b
{
    int n = rows;
    std::vector<double> b(n);

    for (int i = 0; i < n; i++)
        b[i] = data[i][cols - 1];

    for (int i = 0; i < n; i++)
    {
        double sum = b[i];

        for (int j = 0; j < i; j++)
            sum -= L[i][j] * y[j];//subtracts already known values.

        y[i] = sum / L[i][i];
    }
}

void LUSolver::backwardSubstitution(std::vector<std::vector<double>> &U,
                                    std::vector<double> &y)//solve this system U x = y
{
    int n = rows;
    solution.assign(n, 0);

    for (int i = n - 1; i >= 0; i--)
    {
        double sum = y[i];

        for (int j = i + 1; j < n; j++)
            sum -= U[i][j] * solution[j];//Subtract known values

        if (fabs(U[i][i]) < 1e-10)
            throw runtime_error("Zero pivot during backward substitution");
        solution[i] = sum / U[i][i];
    }
}

void LUSolver::pivotRows(int k)//this will perform partial pivoting.
{
    int maxRow = k;
    double maxVal = fabs(data[k][k]);

    for (int i = k + 1; i < rows; i++)
    {
        if (fabs(data[i][k]) > maxVal)
        {
            maxVal = fabs(data[i][k]);
            maxRow = i;
        }
    }

    if (maxRow != k)
    {
        for (int j = 0; j < cols; j++)
            std::swap(data[k][j], data[maxRow][j]);
    }
}

void LUSolver::croutDecomposition()
{
    int n = rows;

    //crout will split matrix into A = L * U
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0));
    std::vector<std::vector<double>> U(n, std::vector<double>(n, 0));
    std::vector<double> y(n);

    for (int j = 0; j < n; j++)
    {
        pivotRows(j);//Partial pivoting.
        U[j][j] = 1;//Set diagonal of U to 1 because Crout's method requires it.

        for (int i = j; i < n; i++)//Compute column j of L.
        {
            double sum = 0;
            for (int k = 0; k < j; k++)
                sum += L[i][k] * U[k][j];//compute known values

            L[i][j] = data[i][j] - sum;//compute L[i][j]
        }

        for (int i = j + 1; i < n; i++)//Compute row j of U.
        {
            double sum = 0;
            for (int k = 0; k < j; k++)
                sum += L[j][k] * U[k][i];//compute lower triangular part

            if (fabs(L[j][j]) < 1e-10)
                throw runtime_error("Zero pivot in LU decomposition");

            U[j][i] = (data[j][i] - sum) / L[j][j];//compute upper triangular part
        }
    }

    forwardSubstitution(L, y);//Solves Ly = b
    backwardSubstitution(U, y);//Solves Ux = y
}

void LUSolver::doolittleDecomposition()
{
    int n = rows;

    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0));//Create Lower triangular matrix L.
    std::vector<std::vector<double>> U(n, std::vector<double>(n, 0));//Create Upper triangular matrix U.
    std::vector<double> y(n);

    for (int i = 0; i < n; i++)
    {
        pivotRows(i);
        L[i][i] = 1;//set diagonal of L to 1 because Doolittle's method requires it.

        for (int j = i; j < n; j++)//Iterate upper triangular part
        {
            double sum = 0;
            for (int k = 0; k < i; k++)//Iterate known values
                sum += L[i][k] * U[k][j];//compute upper triangular part

            U[i][j] = data[i][j] - sum;//compute U[i][j]
        }

        for (int j = i + 1; j < n; j++)//Iterate lower triangular part
        {
            double sum = 0;
            for (int k = 0; k < i; k++)//Iterate known values
                sum += L[j][k] * U[k][i];//compute lower triangular part

            if (fabs(U[i][i]) < 1e-10)
                throw runtime_error("Zero pivot in Doolittle decomposition");

            L[j][i] = (data[j][i] - sum) / U[i][i];//compute L[j][i]
        }
    }

    forwardSubstitution(L, y);//Solves Ly = b
    backwardSubstitution(U, y);//Solves Ux = y
}

void LUSolver::choleskyDecomposition()
{
    int n = rows;

    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0));
    std::vector<double> y(n);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= i; j++)//L is lower triangular.
        {
            double sum = 0;

            for (int k = 0; k < j; k++)
                sum += L[i][k] * L[j][k];//formula for Cholesky

            if (i == j)//compute diagonal elements
            {
                double val = data[i][i] - sum;

                if (val <= 0)
                    throw runtime_error("Matrix not positive definite");

                L[i][j] = sqrt(val);
            }
            else
                L[i][j] = (data[i][j] - sum) / L[j][j];
        }
    }

    std::vector<std::vector<double>> U = L;

    for (int i = 0; i < n; i++)//Transpose Upper Part
        for (int j = i + 1; j < n; j++)//Iterate upper triangular
            U[i][j] = L[j][i];//Assign Transpose

    forwardSubstitution(L, y);//Solves Ly = b
    backwardSubstitution(U, y);//Solves Ux = y
}