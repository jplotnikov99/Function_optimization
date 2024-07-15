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
    double xi, xf;
    double min_epsilon = 1e100;
    size_t N_grids = 10;
    vec2d weights;
    vec2d N_weights;
    std::vector<double> opt_c;

public:
    Optimizer(std::unique_ptr<Integrator> &Inte, const double start, const double finish);
    double get_min_epsilon();
    std::vector<double> get_opt_c();
    double epsilon();
    void init_grid();
    void print_grid_row(const size_t i);
    void print_grid();
    void update_grid(const vec1d &lower, const vec1d &upper, const vec1d &constants, const double eps);
    void monte_carlo(const vec1d lower, const vec1d upper, const size_t N);
    ~Optimizer(){};
};