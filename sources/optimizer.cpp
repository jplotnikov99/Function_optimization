#include "optimizer.hpp"

Optimizer::Optimizer(std::unique_ptr<Integrator> &inte, const double start, const double finish)
{
    I = std::move(inte);
    xi = start;
    xf = finish;
}

double Optimizer::get_min_epsilon()
{
    return min_epsilon;
}

double Optimizer::epsilon()
{
    double p = I->F->get_p_value();
    return pow(I->adap_gauss_kronrod_15(xi, xf), 1 / p);
}

void Optimizer::monte_carlo(const double lower, const double upper, const size_t N)
{
    double cur_epsilon;
    for (size_t i = 0; i < N; i++)
    {
        I->F->randomize_constants(lower, upper);
        cur_epsilon = epsilon();
        if (cur_epsilon < min_epsilon)
            min_epsilon = cur_epsilon;
    }
}