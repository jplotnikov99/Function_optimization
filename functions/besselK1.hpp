#pragma once

#include <iostream>
#include "../include/utils.hpp"

class BesselK1
{
private:
    Output_type ot = result;
    const size_t p_value;
    const size_t N_coeffs = 5;
    size_t cur_ci = 0;
    vec1d c;
    vec1d g;

public:
    BesselK1(const size_t p) : p_value(p)
    {
        for (size_t i = 0; i < N_coeffs; i++)
        {
            c.push_back(0);
            g.push_back(0);
        }
    };
    void switch_to_res();
    void switch_to_grad();
    void select_cur_ci(const size_t i);
    double get_p_value();
    vec1d get_coeffs();
    size_t get_N_coeffs();
    void change_coeff(const size_t i, const double new_val);
    void change_all_coeffs(const vec1d &new_vals);
    void print_coeffs();

    double besselK1_exact(const double x);
    double besselK1_appr(const double x);
    double besselK1_grad(const double x);
    bool is_valid();

    double operator()(const double x);
    ~BesselK1() {};
};
