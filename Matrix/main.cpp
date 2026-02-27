#include "include/matrix.hpp"
#include <sstream>
#include <vector>

int main()
{
    try
    {
        int n = 225; 
        Matrix m(n, n + 1);

        m.readFromFile("input1.txt");

        m.gaussianElimination(true);

        auto solution = m.backwardSubstitution();

        ofstream out("output1.txt");

        for (double val : solution)
            out << val << endl;

        out.close();

        cout << "Solution written to output1.txt\n";
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
