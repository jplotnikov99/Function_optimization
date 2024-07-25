#pragma once

#include "function.hpp"

class BesselK1 : public Function
{
private:

public:
    BesselK1(const size_t p)
    {
        N_coeffs = 5;
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
    ~BesselK1() {};
};
