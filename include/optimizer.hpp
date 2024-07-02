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
    double xi, xf;
    double min_epsilon = 1e100;
public:
    Optimizer(std::unique_ptr<Integrator> &Inte, const double start, const double finish);
    double get_min_epsilon();
    double epsilon();
    void monte_carlo(const std::vector<double> lower, const std::vector<double> upper, const size_t N);
    ~Optimizer(){};
};