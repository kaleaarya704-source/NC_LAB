#include "../include/gaussian_solver.hpp"
#include <cmath>
#include <stdexcept>

GaussianSolver::GaussianSolver(int n)
    : SystemOfLinearEquation(n)
{
}

void GaussianSolver::solve()
{
    gaussianElimination(true);
    solution = backwardSubstitution();
}




































// #include "../include/matrix.hpp"
// #include <stdexcept>
// #include <cmath>

// void Matrix::gaussianElimination(bool pivoting)
// {
//     const double EPS = 1e-10;

//     for (int k = 0; k < rows - 1; k++)
//     {
//         if (pivoting)
//         {
//             int maxRow = k;
//             double maxVal = fabs(data[k][k]);

//             for (int i = k + 1; i < rows; i++)
//             {
//                 if (fabs(data[i][k]) > maxVal)
//                 {
//                     maxVal = fabs(data[i][k]);
//                     maxRow = i;
//                 }
//             }

//             if (maxVal < EPS)
//                 continue;  //avoid division by zero

//             if (maxRow != k)
//             {
//                 for (int j = 0; j < cols; j++)
//                     swap(data[k][j], data[maxRow][j]);
//             }
//         }
//         else
//         {
//             if (fabs(data[k][k]) < EPS)
//                 continue;
//         }

//         for (int i = k + 1; i < rows; i++)
//         {
//             if (fabs(data[k][k]) < EPS)
//                 continue;

//             double factor = data[i][k] / data[k][k];

//             for (int j = k; j < cols; j++)
//                 data[i][j] -= factor * data[k][j];
//         }
//     }
// }

// vector<double> Matrix::backwardSubstitution() const
// {
//     const double EPS = 1e-10;

//     if (cols != rows + 1)
//         throw runtime_error("Not an augmented matrix.");

//     vector<double> x(rows, 0);

//     for (int i = rows - 1; i >= 0; i--)
//     {
//         double sum = data[i][cols - 1];

//         for (int j = i + 1; j < rows; j++)
//             sum -= data[i][j] * x[j];

//         if (fabs(data[i][i]) < EPS)
//         {
//             x[i] = 0;
//             continue;
//         }

//         x[i] = sum / data[i][i];
//     }

//     return x;
// }