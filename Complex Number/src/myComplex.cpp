#include <iostream>
#include <cmath>
#include "myComplex.hpp"
using namespace std;

// default constructor
myComplex::myComplex()
{
    a = 0;
    b = 0;
}

// parameterised constructor
myComplex::myComplex(float x, float y)
{
    a = x;
    b = y;
}

// addition
myComplex myComplex::add(myComplex c)
{
    myComplex ans;
    ans.a = a + c.a;
    ans.b = b + c.b;
    return ans;
}

// multiplication
myComplex myComplex::multiply(myComplex c)
{
    myComplex ans;
    ans.a = a * c.a - b * c.b;
    ans.b = a * c.b + b * c.a;
    return ans;
}

// division
myComplex myComplex::divide(myComplex c)
{
    myComplex ans;
    float den = c.a * c.a + c.b * c.b;

    ans.a = (a * c.a + b * c.b) / den;
    ans.b = (b * c.a - a * c.b) / den;

    return ans;
}

// conjugate
myComplex myComplex::conjugate()
{
    myComplex ans;
    ans.a = a;
    ans.b = -b;
    return ans;
}

// norm
float myComplex::norm()
{
    return sqrt(a * a + b * b);
}

//display
// void myComplex :: display()
// {
//     cout << a << " + " << b << " i " << endl;
// }
