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
    virtual ~Matrix();

    // File operations
    void readFromFile(const string &filename);
    void writeToFile(const string &filename) const;

    // Basic matrix operations
    Matrix add(const Matrix &other) const;
    Matrix subtract(const Matrix &other) const;
    Matrix multiply(const Matrix &other) const;

    // Operator Overloading
    Matrix operator+(const Matrix &other) const;
    Matrix operator-(const Matrix &other) const;
    Matrix operator*(const Matrix &m) const;
    bool operator==(const Matrix &m) const;

    // Element access
    double &operator()(int i, int j);
    double operator()(int i, int j) const;

    // Input / Output operators
    friend istream &operator>>(istream &in, Matrix &m);
    friend ostream &operator<<(ostream &out, const Matrix &m);

    // Display
    void display() const;

    // Size info
    int getRows() const;
    int getCols() const;

    // Matrix property checks (only for square matrices)
    bool isSquare() const;
    bool isSymmetric() const;
    bool isIdentity() const;
    bool isNull() const;
    bool isDiagonal() const;
    bool isDiagonallyDominant() const;
    bool isTranspose(const Matrix &m) const;

    // Matrix transformations
    Matrix transpose() const;

    // Determinant / Inverse (only for square matrices)
    double determinant() const;
    Matrix inverse() const;

    // Numerical methods operations
    void gaussianElimination(bool pivoting);
    vector<double> backwardSubstitution() const;
};

#endif