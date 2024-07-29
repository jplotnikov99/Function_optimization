#pragma once

#include "function.hpp"
#include "../include/integrator.hpp"

class dD0dw
{
private:
    double s;
    double x;

public:
    dD0dw(const Particle_type pt)
    {
        s = pt == fermion ? 1 : -1;
    };
    void set_x(const double new_x);
    double operator()(const double u);
    ~dD0dw() {};
};

class D0 : public Function
{
private:
    double a1;
    const double a2 = 12 / (M_PI * M_PI);
    dD0dw integrand;
    std::unique_ptr<Integrator> I;
    std::unordered_map<double, double> exact_map;

public:
    D0(const size_t p, const Particle_type pt) : integrand(pt)
    {
        I = std::make_unique<Integrator>(gauss15, 0., 1.);
        a1 = pt == fermion ? 1 : 2; 
        N_coeffs = 11;
        p_value = p;
        for (size_t i = 0; i < N_coeffs; i++)
        {
            c.push_back(0);
        }
    };

    double exact(const double x);
    double approx(const double x);
    double gradient(const double x);
    bool is_valid();

    double operator()(const double x);

    ~D0() {};
};