#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include "utils.hpp"

typedef std::vector<double> vec1d;

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
    vec1d g;

public:
    Function(const func_name n, const double p);
    // prepares the constants and gradients which we want to optimize for
    void prepare();
    double get_p_value();
    std::vector<double> get_coeffs();
    size_t get_N_coeffs();
    void print_coeffs();
    void change_constant(const size_t i, const double new_val);
    void randomize_constants(const double l, const double r);
    double res(const double x);
    double grad(const double x);
    bool is_valid();
    ~Function(){};

    // functions that we want to optimize:
    // _exact: exact function
    // _appr: approximate function
    // _check: checks conditions the constants have to fulfill
    // _grad: gradient with respect to the constants we want to optimize
    double besselK1_exact(const double x);
    double besselK1_appr(const double x);
    bool besselK1_valid();
    vec1d besselK1_grad(const double x);
};