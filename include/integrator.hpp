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
   int_method method;
   std::unique_ptr<Function> F;

public:
   Integrator(std::unique_ptr<Function> &function, const int_method m);
   double adap_gauss_kronrod_15(double a, double b);
   ~Integrator(){};
};
