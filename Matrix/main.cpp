#include "gaussian_solver.hpp"
#include "lu_solver.hpp"
#include <iostream>
#include <fstream>

using namespace std;

int main()
{
    try
    {
        int n = 49;
        int choice;

        cout << "Choose Method:\n";
        cout << "1. Gaussian Elimination\n";
        cout << "2. LU - Crout\n";
        cout << "3. LU - Doolittle\n";
        cout << "4. LU - Cholesky\n";
        cout << "Enter choice: ";
        cin >> choice;

        SystemOfLinearEquation *solver = nullptr;//solver base class ka pointer but it will store obj of derived class (GaussianSolver or LUSolver) at runtime (runtime polymorphism achieved)

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

        solver->createAugmentedMatrix("225left.txt", "225right.txt", "input1.txt");
        solver->readFromFile("input1.txt");
        solver->solve();

        auto solution = solver->getSolution();

        ofstream out("output1.txt");

        for (int i = 0; i < solution.size(); i++)
        {
            out << "x" << i + 1 << " = " << solution[i] << endl;
        }
        out.close();

        delete solver;

        cout << "Augmented matrix written to input1.txt\n";
        cout << "Solution written to output1.txt\n";
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
