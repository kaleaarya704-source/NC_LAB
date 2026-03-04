#include <fstream>
#include <iostream>
using namespace std;

int main()
{
    ifstream x("output.txt");
    ifstream b("49r.txt");
    ofstream out("combined.txt");

    double xv, bv;

    while (x >> xv && b >> bv)
        out << xv << " " << bv << endl;

    cout << "combined1.txt created\n";
}