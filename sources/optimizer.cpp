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

vec1d Optimizer::get_opt_c()
{
    return opt_c;
}

double Optimizer::epsilon()
{
    double p = I->F->get_p_value();
    return pow(I->integrate(xi, xf), 1 / p);
}

void Optimizer::init_grid()
{
    size_t NC = I->F->get_N_c();
    for (size_t i = 0; i < NC; i++)
    {
        weights.push_back({});
        N_weights.push_back({});
        for (size_t j = 0; j < N_grids; j++)
        {
            weights.at(i).push_back(0.);
            N_weights.at(i).push_back(0.);
        }
    }
}

void Optimizer::print_grid_row(const size_t i)
{
    assert(i < I->F->get_N_c());
    for(auto it : weights.at(i))
    {
        std::cout << it << "\t";
    }
    std::cout << "\n";
}

void Optimizer::print_grid()
{
    for(size_t i = 0; i < I->F->get_N_c(); i++)
        print_grid_row(i);
}

void Optimizer::update_grid(const vec1d &lower, const vec1d &upper, const vec1d &constants, const double eps)
{
    size_t NC = I->F->get_N_c();
    assert(constants.size() == NC);
    size_t index{0};

    for (size_t i = 0; i < NC; i++)
    {
        index = std::ceil((constants.at(i) - lower.at(i)) * (double)N_grids / (upper.at(i) - lower.at(i))) - 1;

        weights[i][index] *= (double)N_weights[i][index];
        N_weights[i][index]++;
        weights[i][index] = (weights[i][index] + eps) / (double)N_weights[i][index];
    }
}

void Optimizer::monte_carlo(const vec1d lower, const vec1d upper, const size_t N)
{
    assert(lower.size() == upper.size());
    bool passed;
    double cur_epsilon;
    vec1d cur_c;

    init_grid();

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < lower.size(); j++)
            I->F->change_constant(j, generate_random(lower.at(j), upper.at(j)));

        passed = I->F->is_valid();
        if (passed)
        {
            cur_c = I->F->get_constants();
            cur_epsilon = epsilon();
            update_grid(lower, upper, cur_c, cur_epsilon);

            if (cur_epsilon < min_epsilon)
            {
                opt_c = cur_c;
                min_epsilon = cur_epsilon;
            }
        }
        else
        {
            i--;
        }
    }
}