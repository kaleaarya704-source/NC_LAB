#include <iostream>
#include "myComplex.hpp"
using namespace std;

int main()
{
    myComplex c1(3, 4);
    myComplex c2(1, 2);

    myComplex sum = c1.add(c2);
    cout << "Addition : " << sum.a << " + " << sum.b << "i\n";
    // cout << "Addition : " ; sum.display();

    myComplex mul = c1.multiply(c2);
    cout << "Multiplication : " << mul.a << " + " << mul.b << "i\n";

    myComplex div = c1.divide(c2);
    cout << "Division : " << div.a << " + " << div.b << "i\n";

    myComplex con = c1.conjugate();
    cout << "Conjugate of c1 : " << con.a << " + " << con.b << "i\n";

    cout << "Norm of c1 : " << c1.norm() << endl;

    return 0;
}
