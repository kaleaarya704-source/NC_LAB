#include<vector>
#include "include/iterative_solver.hpp"
#include "gaussian_solver.hpp"
#include "lu_solver.hpp"
#include "matrix.hpp" 
#include <iostream>
#include <fstream>

using namespace std;

int main()
{
    try
    {
        int n = 225; // size of system
        int choice;

        Matrix augmented(n, n + 1);//create augmented matrix
        augmented.readFromFile("input1.txt");

        Matrix A(n, n);//create coefficient matrix

        for (int i = 0; i < n; i++)//copies coefficient part to A for property checks
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
        cout << "Enter choice: ";

        cin >> choice;

        SystemOfLinearEquation *solver = nullptr;//pointer to base class for polymorphism

        vector<double> solution;  // for iterative methods
        if (choice >= 1 && choice <= 4)
        {
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

            delete solver;
        }
        else if (choice == 5 || choice == 6)
        {
            // Extract b vector
            vector<double> b(n);
            for (int i = 0; i < n; i++)
                b[i] = augmented(i, n);//iterative methods need a nad b separately

            // Check diagonal dominance
            if (!IterativeSolver::isDiagonallyDominant(A))
            {
                cout << "Matrix is NOT diagonally dominant. May not converge.\n";
            }
            else
            {
                cout << "Matrix is diagonally dominant.\n";
            }

            if (choice == 5)
            {
                solution = IterativeSolver::gaussJacobi(A, b, 100, 1e-6);
            }
            else
            {
                solution = IterativeSolver::gaussSeidel(A, b, 100, 1e-6);
            }
        }
        else
        {
            cout << "Invalid choice\n";
            return 0;
        }

        if (solver != nullptr)
        {
            solver->readFromFile("input1.txt");//loads augmented matrix into solver object
            solver->solve();//this will execute selected method

        }

        ofstream out("output1.txt");

        for (int i = 0; i < solution.size(); i++)//write solution to file
        {
            out << "x" << i + 1 << " = " << solution[i] << endl;
        }

        out.close();

        delete solver;

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