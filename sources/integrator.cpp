#include "integrator.hpp"

Integrator::Integrator(std::unique_ptr<Function> &function, const int_method m)
{
    F = std::move(function);
    method = m;
}

double Integrator::adap_gauss_kronrod_15(double a, double b)
{
    return F->res(a);
}