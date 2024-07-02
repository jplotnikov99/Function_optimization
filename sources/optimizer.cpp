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
    return pow(I->integrate(xi, xf), 1 / p);
}

void Optimizer::monte_carlo(const std::vector<double> lower, const std::vector<double> upper, const size_t N)
{
    assert(lower.size() == upper.size());
    bool passed;
    double cur_epsilon;

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < lower.size(); j++)
            I->F->change_constant(j, generate_random(lower.at(j), upper.at(j)));

        passed = I->F->is_valid();
        if (passed)
        {
            cur_epsilon = epsilon();
            if (cur_epsilon < min_epsilon)
                min_epsilon = cur_epsilon;
        }
        else
        {
            i--;
        }
    }
}