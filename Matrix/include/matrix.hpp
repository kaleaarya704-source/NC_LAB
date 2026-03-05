#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class Matrix
{
protected:
    int rows;
    int cols;
    double **data;

public:
    Matrix(int r, int c);
    Matrix(const Matrix &other);
    virtual ~Matrix(); // make virtual (important for inheritance)

    // Keep your existing functions
    void readFromFile(const string &filename);
    void writeToFile(const string &filename) const;

    Matrix add(const Matrix &other) const;
    Matrix subtract(const Matrix &other) const;

    // Operator Overloading
    Matrix operator+(const Matrix &other) const;
    Matrix operator-(const Matrix &other) const;

    void display() const;

    int getRows() const;
    int getCols() const;

    void gaussianElimination(bool pivoting); //Transforms matrix into upper triangular form.
    vector<double> backwardSubstitution() const; //Solves triangular matrix.
};

#endif