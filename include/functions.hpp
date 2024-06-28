#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include "utils.hpp"

typedef std::vector<double> gradient;

enum func_name
{
    test,
    besselK1,
    besselK2
};

class Function
{
private:
    func_name name;
    double p_value; 
    std::vector<double> c;
    gradient g;
public:
    Function(const func_name n, const double p);
    // prepares the constants and gradients which we want to optimize for 
    double get_p_value();
    size_t get_N_constants();
    void change_constant(const size_t i, const double new_val);
    void randomize_constants(const double l, const double r);
    void prepare();
    double res(const double x);
    double eps(const double xi, const double xf);
    ~Function(){};

    // functions that we want to optimize:
    // _exact: exact function
    // _appr: approximate function
    // _check: checks conditions the constants have to fulfill
    // _grad: gradient with respect to the constants we want to optimize
    double besselK1_exact(double x);
    double besselK1_appr(double x);
    bool besselK1_check();
    gradient besselK1_grad(double x);
    
};