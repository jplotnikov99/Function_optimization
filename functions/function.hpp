#pragma once
#include <iostream>
#include "../include/utils.hpp"

class Function
{
public:
    Output_type ot = result;
    size_t p_value;
    size_t N_coeffs;
    size_t cur_ci = 0;
    vec1d c;

    Function() {};

    void switch_to_res();
    void switch_to_grad();
    void select_cur_ci(const size_t i);
    double get_p_value();
    vec1d get_coeffs();
    size_t get_N_coeffs();
    void change_coeff(const size_t i, const double new_val);
    void change_all_coeffs(const vec1d &new_vals);
    void print_coeffs();

    ~Function() {};
};