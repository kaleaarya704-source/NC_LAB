#include <fstream>
#include <iostream>
using namespace std;

int main()
{
    ifstream x("output1.txt");
    ifstream b("225right.txt");
    ofstream out("combined1.txt");

    double xv, bv;

    while (x >> xv && b >> bv)
        out << xv << " " << bv << endl;

    cout << "combined1.txt created\n";
}