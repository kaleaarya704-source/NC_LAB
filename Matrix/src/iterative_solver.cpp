#include "../include/iterative_solver.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>

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

// Attempts to make A (and b) diagonally dominant by row swapping.
// Uses a greedy assignment: for each column position i, find an unassigned row
// whose absolute diagonal element (column i) dominates its row sum.
// Returns true if a valid permutation is found and applied, false otherwise.
bool IterativeSolver::makeDiagonallyDominant(Matrix &A, vector<double> &b)
{
    int n = A.getRows();

    // For each row r, compute sum of absolute values of all elements
    vector<double> rowSum(n, 0.0);
    for (int r = 0; r < n; r++)
        for (int c = 0; c < n; c++)
            rowSum[r] += abs(A(r, c));

    // Try to find a permutation: assign row perm[i] to position i
    // such that A(perm[i], i) is dominant in that row
    vector<int> perm(n, -1);   // perm[i] = which original row goes to position i
    vector<bool> used(n, false);

    for (int col = 0; col < n; col++)
    {
        int bestRow = -1;
        double bestDiag = -1.0;

        for (int r = 0; r < n; r++)
        {
            if (used[r])
                continue;

            double diag = abs(A(r, col));
            double offDiagSum = rowSum[r] - diag;

            // Check if placing row r at position col satisfies dominance for this position
            if (diag >= offDiagSum && diag > bestDiag)
            {
                bestDiag = diag;
                bestRow = r;
            }
        }

        if (bestRow == -1)
            return false; // No valid permutation found

        perm[col] = bestRow;
        used[bestRow] = true;
    }

    // Apply permutation: reorder rows of A and elements of b
    int cols = A.getCols();
    // Build reordered copies
    Matrix Anew(n, cols);
    vector<double> bnew(n);

    for (int i = 0; i < n; i++)
    {
        int src = perm[i];
        for (int j = 0; j < cols; j++)
            Anew(i, j) = A(src, j);
        bnew[i] = b[src];
    }

    // Copy back into A and b
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < cols; j++)
            A(i, j) = Anew(i, j);
        b[i] = bnew[i];
    }

    return true;
}

//Gauss-Jacobi with optional SOR (Successive Over-Relaxation)
// omega = 1.0  -> standard Jacobi
// omega < 1.0  -> under-relaxation (stabilises diverging matrices)
vector<double> IterativeSolver::gaussJacobi(const Matrix &A, const vector<double> &b, int maxIter, double tol, double omega)
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
            double x_jacobi = (b[i] - sum) / A(i, i);//standard Jacobi update
            x_new[i] = (1.0 - omega) * x[i] + omega * x_jacobi;//SOR blend
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