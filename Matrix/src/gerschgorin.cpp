#include "../include/gerschgorin.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <stdexcept>

using namespace std;

// ---------------------------------------------------------------
GerschgorinSolver::GerschgorinSolver(const Matrix &A)
    : EigenAnalysis(A),
      unionLo( numeric_limits<double>::max()),
      unionHi(-numeric_limits<double>::max()),
      computed(false)
{
}

// ---------------------------------------------------------------
// Gerschgorin Circle Theorem
//   centre  c_i = A(i,i)
//   radius  r_i = sum_{j != i} |A(i,j)|
//   disc    D_i -> real interval [c_i - r_i , c_i + r_i]
// ---------------------------------------------------------------
void GerschgorinSolver::compute()
{
    discs.clear();
    unionLo =  numeric_limits<double>::max();
    unionHi = -numeric_limits<double>::max();

    for (int i = 0; i < n; i++)
    {
        Disc d;
        d.row    = i;
        d.centre = (*this)(i, i);   // diagonal element

        // radius = sum of absolute off-diagonal values in row i
        d.radius = 0.0;
        for (int j = 0; j < n; j++)
            if (j != i)
                d.radius += fabs((*this)(i, j));

        d.lo = d.centre - d.radius;
        d.hi = d.centre + d.radius;

        discs.push_back(d);

        // Expand union interval
        if (d.lo < unionLo) unionLo = d.lo;
        if (d.hi > unionHi) unionHi = d.hi;
    }

    computed = true;
}

// ---------------------------------------------------------------
void GerschgorinSolver::printResults() const
{
    if (!computed)
        throw runtime_error("GerschgorinSolver: call compute() before printResults().");

    cout << fixed << setprecision(4);
    cout << "\n=== Gerschgorin Circle Theorem Results ===\n";
    cout << "Matrix order: " << n << " x " << n << "\n\n";

    cout << left
         << setw(6)  << "Row"
         << setw(14) << "Centre (c_i)"
         << setw(14) << "Radius (r_i)"
         << setw(14) << "Lo (c-r)"
         << setw(14) << "Hi (c+r)"
         << "\n";
    cout << string(62, '-') << "\n";

    for (const auto &d : discs)
    {
        cout << setw(6)  << (d.row + 1)
             << setw(14) << d.centre
             << setw(14) << d.radius
             << setw(14) << d.lo
             << setw(14) << d.hi
             << "\n";
    }

    cout << string(62, '-') << "\n";
    cout << "\nUnion of all discs (eigenvalue bounding interval):\n";
    cout << "  [ " << unionLo << " , " << unionHi << " ]\n";
    cout << "\nAll eigenvalues of the matrix lie within this interval.\n";
}

// ---------------------------------------------------------------
const vector<GerschgorinSolver::Disc> &GerschgorinSolver::getDiscs() const
{
    return discs;
}

double GerschgorinSolver::getUnionLo() const { return unionLo; }
double GerschgorinSolver::getUnionHi() const { return unionHi; }
