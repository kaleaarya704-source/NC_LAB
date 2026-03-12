#include "../include/matrix.hpp"
#include <stdexcept>
#include <cmath>
#include <iomanip>

Matrix::Matrix()
{
    rows = 0;
    cols = 0;
    data = nullptr;
}

Matrix::Matrix(int r, int c)
{
    if (r <= 0 || c <= 0)
        throw invalid_argument("Invalid matrix size.");

    rows = r;
    cols = c;

    data = new double *[rows]; // Creates array of row pointers.
    for (int i = 0; i < rows; i++)
        data[i] = new double[cols]; // Creates each row.
}

Matrix::Matrix(const Matrix &other)
{
    rows = other.rows;
    cols = other.cols;

    data = new double *[rows];
    for (int i = 0; i < rows; i++)
    {
        data[i] = new double[cols];
        for (int j = 0; j < cols; j++)
            data[i][j] = other.data[i][j];
    }
}

Matrix::~Matrix()
{
    for (int i = 0; i < rows; i++)
        delete[] data[i];
    delete[] data;
}

void Matrix::readFromFile(const string &filename)
{
    ifstream file(filename);
    if (!file)
        throw runtime_error("Cannot open input file.");

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            file >> data[i][j];

    file.close();
}



void Matrix::writeToFile(const string &filename) const
{
    ofstream file(filename);
    if (!file)
        throw runtime_error("Cannot open output file.");

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
            file << data[i][j] << " ";
        file << endl;
    }
    file.close();
}

Matrix Matrix::add(const Matrix &other) const
{
    if (rows != other.rows || cols != other.cols)
        throw runtime_error("Addition size mismatch.");

    Matrix result(rows, cols);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result.data[i][j] = data[i][j] + other.data[i][j];

    return result;
}

Matrix Matrix::subtract(const Matrix &other) const
{
    if (rows != other.rows || cols != other.cols)
        throw runtime_error("Subtraction size mismatch.");

    Matrix result(rows, cols);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result.data[i][j] = data[i][j] - other.data[i][j];

    return result;
}

Matrix Matrix::multiply(const Matrix &other) const
{
    return (*this) * other;
}

Matrix Matrix::operator+(const Matrix &other) const
{
    return add(other); // reuse existing function
}

Matrix Matrix::operator-(const Matrix &other) const
{
    return subtract(other); // reuse existing function
}

double &Matrix::operator()(int i, int j)//Returns reference to matrix element
{
    return data[i][j];
}

double Matrix::operator()(int i, int j) const//Used when matrix cannot be modified
{
    return data[i][j];
}

Matrix Matrix::operator*(const Matrix &m) const
{
    if (cols != m.rows)
        throw runtime_error("Multiplication size mismatch.");

    Matrix result(rows, m.cols);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < m.cols; j++)
        {
            result.data[i][j] = 0;
            for (int k = 0; k < cols; k++)
                result.data[i][j] += data[i][k] * m.data[k][j];
        }

    return result;
}

bool Matrix::operator==(const Matrix &m) const
{
    if (rows != m.rows || cols != m.cols)
        return false;

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            if (data[i][j] != m.data[i][j])
                return false;

    return true;
}

istream &operator>>(istream &in, Matrix &m)
{
    for (int i = 0; i < m.rows; i++)
        for (int j = 0; j < m.cols; j++)
            in >> m.data[i][j];

    return in;
}

ostream &operator<<(ostream &out, const Matrix &m)
{
    out << fixed << setprecision(4);

    for (int i = 0; i < m.rows; i++)
    {
        for (int j = 0; j < m.cols; j++)
            out << setw(10) << m.data[i][j];

        out << endl;
    }

    return out;
}

bool Matrix::isSquare() const
{
    return rows == cols;
}

bool Matrix::isNull() const
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            if (data[i][j] != 0)
                return false;

    return true;
}

bool Matrix::isIdentity() const
{
    if (!isSquare())
        return false;

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < rows; j++) //only check coefficient part
        {
            if (i == j && data[i][j] != 1)
                return false;
            if (i != j && data[i][j] != 0)
                return false;
        }

    return true;
}

bool Matrix::isSymmetric() const
{
    if (!isSquare())
        return false;

    for (int i = 0; i < rows; i++)
        for (int j = i + 1; j < rows; j++) //ignore augmented column
            if (data[i][j] != data[j][i])
                return false;

    return true;
}

bool Matrix::isDiagonal() const
{
    if (!isSquare())
        return false;

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < rows; j++) 
            if (i != j && data[i][j] != 0)
                return false;

    return true;
}

bool Matrix::isDiagonallyDominant() const
{
    if (!isSquare())
        return false;

    for (int i = 0; i < rows; i++)
    {
        double sum = 0;

        for (int j = 0; j < rows; j++) 
        {
            if (i != j)
                sum += fabs(data[i][j]);
        }

        if (fabs(data[i][i]) < sum)
            return false;
    }

    return true;
}


Matrix Matrix::transpose() const
{
    Matrix result(cols, rows);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result.data[j][i] = data[i][j];

    return result;
}

double Matrix::determinant() const
{
    throw runtime_error("Determinant disabled for large matrices.");
}

Matrix Matrix::inverse() const
{
    throw runtime_error("Inverse computation disabled for large matrices.");
}

void Matrix::gaussianElimination(bool pivoting)
{
    const double EPS = 1e-10;

    for (int k = 0; k < rows - 1; k++)
    {
        if (pivoting)
        {
            int maxRow = k;
            double maxVal = fabs(data[k][k]);

            for (int i = k + 1; i < rows; i++)
            {
                if (fabs(data[i][k]) > maxVal)
                {
                    maxVal = fabs(data[i][k]);
                    maxRow = i;
                }
            }

            if (maxVal < EPS)
                throw runtime_error("Singular matrix");

            if (maxRow != k)
            {
                for (int j = 0; j < cols; j++)
                    swap(data[k][j], data[maxRow][j]);
            }
        }

        if (fabs(data[k][k]) < EPS)
            throw runtime_error("Zero pivot encountered");

        for (int i = k + 1; i < rows; i++)
        {
            double factor = data[i][k] / data[k][k];

            for (int j = k; j < cols; j++)
                data[i][j] -= factor * data[k][j];
        }
    }
}

vector<double> Matrix::backwardSubstitution() const
{
    const double EPS = 1e-10;

    if (cols != rows + 1)
        throw runtime_error("Not an augmented matrix.");

    vector<double> x(rows, 0);

    for (int i = rows - 1; i >= 0; i--)
    {
        double sum = data[i][cols - 1];

        for (int j = i + 1; j < rows; j++)
            sum -= data[i][j] * x[j];

        if (fabs(data[i][i]) < EPS)
        {
            x[i] = 0;
            continue;
        }

        x[i] = sum / data[i][i];
    }

    return x;
}