#include "gaussian_solver.hpp"
#include "lu_solver.hpp"
#include "matrix.hpp" // MODIFICATION: needed to create temporary matrix
#include <iostream>
#include <fstream>

using namespace std;

int main()
{
    try
    {
        int n = 49; // size of system
        int choice;

        // MODIFICATION: create temporary matrix to read and check properties
        Matrix augmented(n, n + 1);
        augmented.readFromFile("input.txt");

        // extract coefficient matrix
        Matrix A(n, n);

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                A(i, j) = augmented(i, j);

        cout << "\nMatrix Properties:\n";

        cout << "Is Square: " << (A.isSquare() ? "Yes" : "No") << endl;
        cout << "Is Symmetric: " << (A.isSymmetric() ? "Yes" : "No") << endl;
        cout << "Is Identity: " << (A.isIdentity() ? "Yes" : "No") << endl;
        cout << "Is Null: " << (A.isNull() ? "Yes" : "No") << endl;
        cout << "Is Diagonal: " << (A.isDiagonal() ? "Yes" : "No") << endl;
        cout << "Is Diagonally Dominant: "
             << (A.isDiagonallyDominant() ? "Yes" : "No") << endl;

        // MODIFICATION: removed determinant and inverse (not feasible for large matrices)

        cout << "\nChoose Method:\n";
        cout << "1. Gaussian Elimination\n";
        cout << "2. LU - Crout\n";
        cout << "3. LU - Doolittle\n";
        cout << "4. LU - Cholesky\n";
        cout << "Enter choice: ";

        cin >> choice;

        SystemOfLinearEquation *solver = nullptr;

        if (choice == 1)
        {
            solver = new GaussianSolver(n);
        }
        else if (choice == 2)
        {
            solver = new LUSolver(n, LUSolver::CROUT);
        }
        else if (choice == 3)
        {
            solver = new LUSolver(n, LUSolver::DOOLITTLE);
        }
        else if (choice == 4)
        {
            solver = new LUSolver(n, LUSolver::CHOLESKY);
        }
        else
        {
            cout << "Invalid choice\n";
            return 0;
        }

        // MODIFICATION: read matrix into solver AFTER creation
        solver->readFromFile("input.txt");

        solver->solve();

        auto solution = solver->getSolution();

        // MODIFICATION: write solution only to file
        ofstream out("output.txt");

        for (int i = 0; i < solution.size(); i++)
        {
            out << "x" << i + 1 << " = " << solution[i] << endl;
        }

        out.close();

        delete solver;

        cout << "\nSolution written to output.txt\n";
    }
    catch (exception &e)
    {
        cout << "Error: " << e.what() << endl;

        ofstream out("output.txt");
        out << "Computation failed.\n";
        out.close();
    }

    return 0;
}