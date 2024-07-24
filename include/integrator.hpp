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
   bool is_grad = false;
   int_method method;

public:
   std::unique_ptr<Function> F;
   Integrator(std::unique_ptr<Function> &function, const int_method m, const double x_ini, const double x_fin);
   void switch_to_grad();
   void switch_to_res();
   double kronrod_61(const double l, const double r, const size_t c_i = 0);
   double adap_gauss_kronrod_15(double l, double r, const double appr, const size_t c_i = 0, const double err = 1e-10);

   double integrate(const size_t c_i = 0, const double err = 1e-10);
   ~Integrator(){};
};
