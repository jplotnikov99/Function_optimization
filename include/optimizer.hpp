#pragma once

#include <iostream>
#include <memory>
#include <vector>
#include <cassert>
#include "integrator.hpp"

typedef std::vector<double> vec1d;
typedef std::vector<vec1d> vec2d;

class Optimizer
{
private:
    std::unique_ptr<Integrator> I;
    size_t N_coeffs;
    double min_epsilon = 1e100;
    size_t N_bins = 10;
    vec2d disjoint_spaces;
    vec2d weights;
    std::vector<double> opt_coeffs;

public:
    Optimizer(std::unique_ptr<Integrator> &Inte, const vec1d &lower, const vec1d &upper);
    double get_min_epsilon();
    std::vector<double> get_opt_coeffs();
    double epsilon();
    void update_grid();
    void print_grid_row(const size_t i);
    void print_grid();
    size_t randomize_coeffs();
    void set_weight(const size_t space, const vec1d &constants, const double eps);
    void monte_carlo(const size_t N);
    void eliminate_weak_grids(const size_t keepers);
    ~Optimizer(){};
};