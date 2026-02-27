class RootFinding
{
public:
    double func(double x);
    double derivative(double x);
    double g(double x);

    virtual double solve(double tol) = 0;
};

class bisectionMethod : public RootFinding
{
public:
    double solve(double tol) override;
};

class newtonRaphson : public RootFinding
{
public:
    double solve(double tol) override;
    double solve(double x0, double tol);  
};

class fixedPointIteration : public RootFinding
{
public:
    double solve(double tol) override;
    double solve(double x0, double tol);
};

