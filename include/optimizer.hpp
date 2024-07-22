#pragma once

#include <iostream>
#include <memory>
#include <vector>
#include <cassert>
#include "integrator.hpp"

class Optimizer
{
private:
    std::unique_ptr<Integrator> I;
    std::string save_file;
    vstring header;
    double min_epsilon = 1e100;
    size_t N_coeffs;
    size_t N_bins = 50, N_spaces = 1;
    vec2d coefficent_spaces;
    vec2d weights;
    std::vector<double> opt_coeffs;

public:
    Optimizer(std::unique_ptr<Integrator> &Inte, const vec1d &lower, const vec1d &upper,
              const std::string file_name);
    
    vec1d get_cur_coeffs();
    void change_coeff(const size_t c_i, const double new_val);
    void set_opt_coeffs();
    vec1d get_opt_coeffs();

    void add_space(const vec1d &lower, const vec1d &upper);
    void reset_space();
    void print_space();

    double epsilon();
    vec1d grad_epsilon(const int coeff = -1);
    double get_min_epsilon();

    // clear weights and update the vector size of the weights for new coefficent space
    void update_grid();

    void print_grid_row(const size_t i);
    void print_grid();

    // returns the space in which the random coefficients were chosen
    size_t randomize_coeffs();

    // return the bin for each coefficient, the space they are in and their epsilon value
    vec1d set_weight(const size_t space, const vec1d &constants, const double eps);

    // makes new search spaces with grids that have the smallest weight
    void make_new_spaces(const vec2d &grids);

    // returns the best N_new_spaces bins for the coefficient and their epsilon value
    void monte_carlo(const size_t N, const size_t N_new_spaces = 1);

    // does a montecarlo search with N_points and N_loops keeping the best N_new_spaces each time
    void repeated_monte_carlo(const size_t N_points, const size_t N_loops, const size_t N_new_spaces);

    vec1d gradient_descent(const int coeff = -1);

    ~Optimizer(){};
};