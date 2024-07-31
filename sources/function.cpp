#include "function.hpp"

void Function::switch_to_res() { ot = result; }

void Function::switch_to_grad() { ot = gradient; }

void Function::select_cur_ci(const size_t i) { cur_ci = i; }

double Function::get_p_value() { return p_value; }

vec1d Function::get_coeffs() { return c; }

size_t Function::get_N_coeffs() { return N_coeffs; }

void Function::change_coeff(const size_t i, const double new_val) {
    assert(i < N_coeffs);
    c[i] = new_val;
}

void Function::change_all_coeffs(const vec1d &new_vals) {
    assert(new_vals.size() == N_coeffs);
    c = new_vals;
}

void Function::print_coeffs() {
    std::cout << "Current coefficents:\n";
    for (auto &it : c) {
        std::cout << it << "\t";
    }
    std::cout << "\n";
}