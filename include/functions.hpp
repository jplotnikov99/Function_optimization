#include <iostream>
#include <cmath>
#include <vector>

typedef std::vector<double> gradient;

enum func_name
{
    test,
    besselK1,
    besselK2
};

class Functions
{
private:
    func_name name;
    double p_value; 
    std::vector<double> c;
    gradient g;
public:
    Functions(const func_name n, const double p);
    // prepares the constants and gradients which we want to optimize for 
    void prepare();
    double res(const double x);
    ~Functions(){};

    // functions that we want to optimize:
    // _exact: exact function
    // _appr: approximate function
    // _grad: gradient with respect to the constants we want to optimize
    double besselK1_exact(double x);
    double besselK1_appr(double x);
    gradient besselK1_grad(double x);
    
};