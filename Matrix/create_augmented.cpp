#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;

int main()
{
    ifstream left("49l.txt");
    ifstream right("49r.txt");
    ofstream out("input.txt");

    if (!left || !right || !out)
    {
        cout << "File error.\n";
        return 1;
    }

    string line;
    vector<double> b;
    double value;

    while (right >> value)
        b.push_back(value);

    int row = 0;

    while (getline(left, line))
    {
        if (line.empty()) continue;

        stringstream ss(line);
        while (ss >> value)
            out << value << " ";

        out << b[row++] << endl;
    }

    cout << "Augmented matrix created.\n";
}
