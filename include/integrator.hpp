#pragma once

#include <iostream>
#include <memory>
#include "abscissa.hpp"
#include "functions.hpp"

enum int_method
{
   gauss15,
   simps38
};

class Integrator
{
private:
   double xi, xf;
   int_method method;

public:
   std::unique_ptr<Function> F;
   Integrator(std::unique_ptr<Function> &function, const int_method m, const double x_ini, const double x_fin);
   double kronrod_61(const double l, const double r);
   double adap_gauss_kronrod_15(double l, double r, const double appr, const double err = 1e-5);

   double integrate(const double err = 1e-5);
   ~Integrator(){};
};
