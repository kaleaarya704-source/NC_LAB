class myComplex
{
public:
    float a, b;   // a = real part, b = imaginary part

    // constructors
    myComplex();
    myComplex(float x, float y);

    // operations
    myComplex add(myComplex c);
    myComplex multiply(myComplex c);
    myComplex divide(myComplex c);
    myComplex conjugate();
    // void display();
    float norm();
};

