#pragma once

#include "../include/integrator.hpp"
#include "function.hpp"

class dD2dw_exact {
   private:
    double s;
    double x;
    double vw;

   public:
    dD2dw_exact(const Particle_type pt) { s = pt == fermion ? 1 : -1; };
    void set_x(const double new_x);
    void set_vw(const double new_vw);
    double operator()(const double u);
    ~dD2dw_exact() {};
};