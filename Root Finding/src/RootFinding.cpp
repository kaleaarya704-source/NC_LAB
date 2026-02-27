#include <cmath>
#include "../include/RootFinding.hpp"

// Function : f(x) = x^3 - x - 2 , f(x) = x^2 - 25 , f(x) = 4x^3 - 3x
double RootFinding::func(double x)
{
    return x * x * x - x - 2;
}

// Derivative
double RootFinding::derivative(double x)
{
    return 3 * x * x - 1;
}

// Fixed Point Function
double RootFinding::g(double x)
{
    return cbrt(x + 2);
}

double bisectionMethod::solve(double tol)
{
    double a = 0, b = 0, m;
    double step = 0.1;
    bool found = false;

    for (double i = 0; i <= 100; i += step)
    {
        if (func(i) * func(i + 1) < 0)
        {
            a = i;
            b = i + 1;
            found = true;
            break;
        }
    }

    if (!found)
        return NAN;

    while ((b - a) >= tol)
    {
        m = (a + b) / 2;

        if (func(a) * func(m) < 0)
            b = m;
        else
            a = m;
    }

    return (a + b) / 2;
}

double newtonRaphson::solve(double tol)
{
    double x0 = 1.0;
    return solve(x0, tol);
}

double newtonRaphson::solve(double x0, double tol)
{
    double x1;

    while (true)
    {
        x1 = x0 - func(x0) / derivative(x0);

        if (fabs(x1 - x0) < tol)
            break;

        x0 = x1;
    }

    return x1;
}

double fixedPointIteration::solve(double tol)
{
    double x0 = 1.0;
    return solve(x0, tol);
}

double fixedPointIteration::solve(double x0, double tol)
{
    double x1;

    while (true)
    {
        x1 = g(x0);

        if (fabs(x1 - x0) < tol)
            break;

        x0 = x1;
    }

    return x1;
}
