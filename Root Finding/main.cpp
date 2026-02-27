#include <iostream>
#include <cmath>
#include "include/RootFinding.hpp"
#include "include/Matrix.hpp"

using namespace std;

int main()
{
    double tol = 0.0001;

    RootFinding* methods[3];

    methods[0] = new bisectionMethod();
    methods[1] = new newtonRaphson();
    methods[2] = new fixedPointIteration();

    cout << "\nFunction : f(x) = x^3 - x - 2\n\n";

    for (int i = 0; i < 3; i++)
    {
        double root = methods[i]->solve(tol);

        if (isnan(root))
            cout << "Method " << i + 1 << " : Interval not found\n";
        else
            cout << "Method " << i + 1 << " Root = " << root << endl;
    }

    for (int i = 0; i < 3; i++)
        delete methods[i];



    //Matrix
    try
    {
        int r ,c;

        cout << "\nEnter rows and columns for matrices : ";
        cin >> r >> c;

        Matrix A(r,c);
        Matrix B(r,c);

        cout << "\nEnter Matrix A \n";
        A.read();

        cout << " \nEnter Matrix B \n";
        B.read();


        
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
    

    return 0;
}
