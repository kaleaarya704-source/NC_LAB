#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

class Matrix
{
private:
    int rows;
    int cols;
    double **data;

public:
    Matrix(int r, int c);
    Matrix(const Matrix &other);
    ~Matrix();

    void readFromFile(const string &filename);
    void writeToFile(const string &filename) const;

    Matrix add(const Matrix &other) const;
    Matrix subtract(const Matrix &other) const;

    void gaussianElimination(bool pivoting);  
    vector<double> backwardSubstitution() const;

    void display() const;
};

