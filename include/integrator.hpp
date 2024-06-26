#include <iostream>
#include "abscissa.hpp"

template<typename F>
double adap_gauss_kronrod_15(F, double a, double b)
{
   return F(a); 
}