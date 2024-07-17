#include "optimizer.hpp"

Optimizer::Optimizer(std::unique_ptr<Integrator> &inte, const vec1d &lower, const vec1d &upper)
{
    I = std::move(inte);
    N_coeffs = I->F->get_N_c();
    assert(N_coeffs == lower.size());
    assert(N_coeffs == upper.size());

    disjoint_spaces.push_back(lower);
    disjoint_spaces.push_back(upper);
}

double Optimizer::get_min_epsilon()
{
    return min_epsilon;
}

vec1d Optimizer::get_opt_coeffs()
{
    return opt_coeffs;
}

double Optimizer::epsilon()
{
    double p = I->F->get_p_value();
    return pow(I->integrate(), 1 / p);
}

void Optimizer::update_grid()
{
    /* weights data structure example for 2 bins, 3 coeffs, 2 disjoint_spaces
        XX XX
        XX XX
        XX XX

        XX XX
        XX XX
        XX XX
    */
    const size_t N_disjoint_spaces = disjoint_spaces.size() / 2;
    weights.clear();
    for (size_t i = 0; i < N_coeffs * N_disjoint_spaces; i++)
    {
        weights.push_back({});
        for (size_t j = 0; j < N_bins; j++)
        {
            weights.at(i).push_back(0.);
        }
    }
}

void Optimizer::print_grid_row(const size_t i)
{
    assert(i < I->F->get_N_c());
    for (auto it : weights.at(i))
    {
        std::cout << it << "\t";
    }
    std::cout << "\n";
}

void Optimizer::print_grid()
{
    for (size_t i = 0; i < I->F->get_N_c(); i++)
        print_grid_row(i);
}

size_t Optimizer::randomize_coeffs()
{
    const size_t N_disjoint_spaces = disjoint_spaces.size() / 2;
    const size_t space = std::floor(generate_random(0., (double)N_disjoint_spaces));

    for (size_t i = 0; i < N_coeffs; i++)
        I->F->change_constant(i, generate_random(disjoint_spaces[2 * space][i], disjoint_spaces[2 * space + 1][i]));

    return space;
}

void Optimizer::set_weight(const size_t space, const vec1d &constants, const double eps)
{
    size_t index{};

    for (size_t i = 0; i < N_coeffs; i++)
    {
        // index = (xi-x0)*Nb/(xf-xi) - 1 + space*Nc
        index = std::ceil((constants.at(i) - disjoint_spaces[2 * space][i]) * (double)N_bins /
                          (disjoint_spaces[2 * space + 1][i] - disjoint_spaces[2 * space][i])) -
                1 + space * N_coeffs;
        weights[i][index] = weights[i][index] == 0 ? eps : std::min(weights[i][index], eps);
    }
}

void Optimizer::monte_carlo(const size_t N)
{
    bool passed;
    double cur_epsilon;
    vec1d cur_c;

    update_grid();

    for (size_t i = 0; i < N; i++)
    {
        size_t space = randomize_coeffs();
        passed = I->F->is_valid();
        if (passed)
        {
            cur_c = I->F->get_constants();
            cur_epsilon = epsilon();
            set_weight(space, cur_c, cur_epsilon);

            if (cur_epsilon < min_epsilon)
            {
                opt_coeffs = cur_c;
                min_epsilon = cur_epsilon;
            }
        }
        else
        {
            i--;
        }
    }
}

void Optimizer::eliminate_weak_grids(const size_t keepers)
{
    assert(keepers < N_bins);
}