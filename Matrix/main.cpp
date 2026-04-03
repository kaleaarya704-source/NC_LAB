#include<cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "include/iterative_solver.hpp"
#include "include/gaussian_solver.hpp"
#include "include/lu_solver.hpp"
#include "include/matrix.hpp"
#include "include/gerschgorin.hpp"
#include <iostream>
#include <fstream>

using namespace std;

int main()
{
    try
    {
        int n = 225; // size of system
        int choice;

        Matrix augmented(n, n + 1); // create augmented matrix
        augmented.readFromFile("input1.txt");

        Matrix A(n, n); // create coefficient matrix

        for (int i = 0; i < n; i++) // copies coefficient part to A for property checks
            for (int j = 0; j < n; j++)
                A(i, j) = augmented(i, j);

        cout << "\nMatrix Properties:\n";

        cout << "Is Square: " << (A.isSquare() ? "Yes" : "No") << endl;
        cout << "Is Symmetric: " << (A.isSymmetric() ? "Yes" : "No") << endl;
        cout << "Is Identity: " << (A.isIdentity() ? "Yes" : "No") << endl;
        cout << "Is Null: " << (A.isNull() ? "Yes" : "No") << endl;
        cout << "Is Diagonal: " << (A.isDiagonal() ? "Yes" : "No") << endl;
        cout << "Is Diagonally Dominant: " << (A.isDiagonallyDominant() ? "Yes" : "No") << endl;

        cout << "\nChoose Method:\n";
        cout << "1. Gaussian Elimination\n";
        cout << "2. LU - Crout\n";
        cout << "3. LU - Doolittle\n";
        cout << "4. LU - Cholesky\n";
        cout << "5. Gauss-Jacobi\n";
        cout << "6. Gauss-Seidel\n";
        cout << "7. Gerschgorin Circle Theorem (eigenvalue bounds)\n";
        cout << "Enter choice: ";

        cin >> choice;

        vector<double> solution; // for all methods

        if (choice >= 1 && choice <= 4)
        {
            SystemOfLinearEquation *solver = nullptr; // pointer to base class for polymorphism

            if (choice == 1)
                solver = new GaussianSolver(n);
            else if (choice == 2)
                solver = new LUSolver(n, LUSolver::CROUT);
            else if (choice == 3)
                solver = new LUSolver(n, LUSolver::DOOLITTLE);
            else if (choice == 4)
                solver = new LUSolver(n, LUSolver::CHOLESKY);

            solver->readFromFile("input1.txt");
            solver->solve();
            solution = solver->getSolution();

            delete solver; // delete only once, after use
        }
        else if (choice == 5 || choice == 6)
        {
            // Extract b vector
            vector<double> b(n);
            for (int i = 0; i < n; i++)
                b[i] = augmented(i, n); // iterative methods need A and b separately

            // Check diagonal dominance and attempt to fix by row swapping
            if (!IterativeSolver::isDiagonallyDominant(A))
            {
                cout << "Matrix is NOT diagonally dominant. Attempting to reorder rows...\n";

                if (IterativeSolver::makeDiagonallyDominant(A, b))
                {
                    cout << "Row reordering successful. Matrix is now diagonally dominant.\n";
                }
                else
                {
                    cout << "Warning: Could not make matrix diagonally dominant by row swapping.\n";
                    cout << "Applying under-relaxation (SOR) for stability.\n";
                }
            }
            else
            {
                cout << "Matrix is diagonally dominant.\n";
            }

            // Compute max spectral radius estimate: max over rows of (off-diag sum / |diag|)
            // omega = 1/(1 + spectral_radius_estimate) ensures convergence via under-relaxation
            double maxRatio = 0.0;
            for (int i = 0; i < n; i++)
            {
                double diag = fabs(A(i, i));
                double offSum = 0.0;
                for (int j = 0; j < n; j++)
                    if (j != i) offSum += fabs(A(i, j));
                if (diag > 0)
                    maxRatio = max(maxRatio, offSum / diag);
            }
            // omega < 1 (under-relaxation) when not diagonally dominant
            double omega = (maxRatio > 1.0) ? 1.0 / (1.0 + maxRatio) : 1.0;
            int iterCount = (maxRatio > 1.0) ? 50000 : 1000;

            if (choice == 5)
            {
                if (omega < 1.0)
                    cout << "Using SOR omega = " << omega << " for Gauss-Jacobi stability.\n";
                solution = IterativeSolver::gaussJacobi(A, b, iterCount, 1e-6, omega);
            }
            else
            {
                solution = IterativeSolver::gaussSeidel(A, b, iterCount, 1e-6);
            }
        }
        else if (choice == 7)
        {
            // --- Gerschgorin Circle Theorem ---
            // Run for both the 49x49 and 225x225 right-hand-side matrices.
            // The theorem requires only the coefficient matrix A (square),
            // so we read only the LEFT (coefficient) files: 49l.txt and 225left.txt.

            // ---- 49 x 49 ----
            {
                cout << "\n--- Gerschgorin: 49 x 49 matrix (49l.txt) ---\n";
                int sz = 49;
                Matrix A49(sz, sz);
                // 49l.txt stores the full augmented matrix (49 x 50); read only coefficient part
                ifstream f49("49l.txt");
                if (!f49)
                    throw runtime_error("Cannot open 49l.txt");
                for (int i = 0; i < sz; i++)
                {
                    for (int j = 0; j < sz; j++)
                        f49 >> A49(i, j);
                    double dummy;
                    f49 >> dummy; // skip the RHS value in the augmented file
                }
                f49.close();

                GerschgorinSolver gs49(A49);
                gs49.compute();
                gs49.printResults();

                // Write per-disc results to file
                ofstream out49("gerschgorin_49.txt");
                out49 << fixed << setprecision(4);
                out49 << "Gerschgorin Discs for 49x49 matrix\n";
                out49 << left
                      << setw(6)  << "Row"
                      << setw(14) << "Centre"
                      << setw(14) << "Radius"
                      << setw(14) << "Lo"
                      << setw(14) << "Hi" << "\n";
                for (const auto &d : gs49.getDiscs())
                    out49 << setw(6)  << (d.row+1)
                          << setw(14) << d.centre
                          << setw(14) << d.radius
                          << setw(14) << d.lo
                          << setw(14) << d.hi << "\n";
                out49 << "\nUnion interval: [ " << gs49.getUnionLo()
                      << " , " << gs49.getUnionHi() << " ]\n";
                out49.close();
                cout << "Results written to gerschgorin_49.txt\n";
            }

            // ---- 225 x 225 ----
            {
                cout << "\n--- Gerschgorin: 225 x 225 matrix (225left.txt) ---\n";
                int sz = 225;
                Matrix A225(sz, sz);
                // 225left.txt stores only the coefficient part (225 x 225)
                ifstream f225("225left.txt");
                if (!f225)
                    throw runtime_error("Cannot open 225left.txt");
                for (int i = 0; i < sz; i++)
                    for (int j = 0; j < sz; j++)
                        f225 >> A225(i, j);
                f225.close();

                GerschgorinSolver gs225(A225);
                gs225.compute();
                gs225.printResults();

                ofstream out225("gerschgorin_225.txt");
                out225 << fixed << setprecision(4);
                out225 << "Gerschgorin Discs for 225x225 matrix\n";
                out225 << left
                       << setw(6)  << "Row"
                       << setw(14) << "Centre"
                       << setw(14) << "Radius"
                       << setw(14) << "Lo"
                       << setw(14) << "Hi" << "\n";
                for (const auto &d : gs225.getDiscs())
                    out225 << setw(6)  << (d.row+1)
                           << setw(14) << d.centre
                           << setw(14) << d.radius
                           << setw(14) << d.lo
                           << setw(14) << d.hi << "\n";
                out225 << "\nUnion interval: [ " << gs225.getUnionLo()
                       << " , " << gs225.getUnionHi() << " ]\n";
                out225.close();
                cout << "Results written to gerschgorin_225.txt\n";
            }

            return 0; // Gerschgorin does not produce a solution vector
        }
        else
        {
            cout << "Invalid choice\n";
            return 0;
        }

        ofstream out("output1.txt");

        for (int i = 0; i < (int)solution.size(); i++) // write solution to file
        {
            out << "x" << i + 1 << " = " << solution[i] << endl;
        }

        out.close();

        cout << "\nSolution written to output1.txt\n";
    }
    catch (exception &e)
    {
        cout << "Error: " << e.what() << endl;

        ofstream out("output1.txt");
        out << "Computation failed.\n";
        out.close();
    }

    return 0;
}
