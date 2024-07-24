#pragma once
#include <iostream>

enum int_method
{
    gauss15,
    simps38
};

class Integrator2
{
private:
    const double x_lower, x_upper;
    const int_method method;

public:
    Integrator2(const int_method m, const double xi, const double xu) : method(m), x_lower(xi), x_upper(xu){};
    double kronrod_61(const double l, const double r);
    double adap_gauss_kronrod_15(double l, double r, const double appr, const double err = 1e-10);
    ~Integrator2(){};
};