#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <fstream>
#include <vector>

class Matrix
{
protected:
    int rows;
    int cols;
    double **data;

public:
    Matrix();
    Matrix(int r, int c);
    Matrix(const Matrix &other);
    virtual ~Matrix();

    // File operations
    void readFromFile(const std::string &filename);
    void writeToFile(const std::string &filename) const;

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
    friend std::istream &operator>>(std::istream &in, Matrix &m);
    friend std::ostream &operator<<(std::ostream &out, const Matrix &m);

    // Display
    void display() const;

    // Size info
    int getRows() const;
    int getCols() const;

    // Matrix property checks
    bool isSquare() const;
    bool isSymmetric() const;
    bool isIdentity() const;
    bool isNull() const;
    bool isDiagonal() const;
    bool isDiagonallyDominant() const;
    bool isTranspose(const Matrix &m) const;

    // Matrix transformations
    Matrix transpose() const;

    // Determinant / Inverse
    double determinant() const;
    Matrix inverse() const;

    // Linear system helpers
    void gaussianElimination(bool pivoting);
    std::vector<double> backwardSubstitution() const;
};

#endif