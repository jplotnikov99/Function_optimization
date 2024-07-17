#include "optimizer.hpp"

Optimizer::Optimizer(std::unique_ptr<Integrator> &inte, const vec1d &lower, const vec1d &upper)
{
    I = std::move(inte);
    N_coeffs = I->F->get_N_coeffs();
    assert(N_coeffs == lower.size());
    assert(N_coeffs == upper.size());

    disjoint_spaces.push_back(lower);
    disjoint_spaces.push_back(upper);
}

void Optimizer::add_space(const vec1d &lower, const vec1d &upper)
{
    disjoint_spaces.push_back(lower);
    disjoint_spaces.push_back(upper);
    N_spaces++;
}

double Optimizer::epsilon()
{
    double p = I->F->get_p_value();
    return pow(I->integrate(), 1 / p);
}

double Optimizer::get_min_epsilon()
{
    return min_epsilon;
}

vec1d Optimizer::get_opt_coeffs()
{
    return opt_coeffs;
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
    weights.clear();
    for (size_t i = 0; i < N_coeffs * N_spaces; i++)
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
    for (auto it : weights.at(i))
    {
        std::cout << it << "\t";
    }
    std::cout << "\n";
}

void Optimizer::print_grid()
{
    for (size_t i = 0; i < N_spaces; i++)
    {
        std::cout << "---------------------------------------" << std::endl;
        for (size_t j = i * N_coeffs; j < (i + 1) * N_coeffs; j++)
            print_grid_row(j);
    }
}

size_t Optimizer::randomize_coeffs()
{
    const size_t space = std::floor(generate_random(0., (double)N_spaces));

    for (size_t i = 0; i < N_coeffs; i++)
        I->F->change_constant(i, generate_random(disjoint_spaces[2 * space][i], disjoint_spaces[2 * space + 1][i]));

    return space;
}

vec1d Optimizer::set_weight(const size_t space, const vec1d &constants, const double eps)
{
    vec1d res(N_coeffs + 2);
    size_t i1 = N_coeffs * space;
    size_t i2;

    for (size_t i = 0; i < N_coeffs; i++)
    {
        // index = (xi-x0)*Nb/(xf-xi) - 1
        i2 = std::ceil((constants.at(i) - disjoint_spaces[2 * space][i]) * (double)N_bins /
                       (disjoint_spaces[2 * space + 1][i] - disjoint_spaces[2 * space][i])) -
             1;
        weights[i1 + i][i2] = weights[i1 + i][i2] == 0 ? eps : std::min(weights[i1 + i][i2], eps);
        res[i] = i2;
    }
    res[N_coeffs] = space;
    res[N_coeffs + 1] = eps;

    return res;
}

vec2d Optimizer::monte_carlo(const size_t N, const size_t N_new_spaces)
{
    vec2d res(N_new_spaces, vec1d(N_coeffs + 2, 1e100));
    double cur_epsilon;
    vec1d cur_coeffs;
    vec1d cur_weight;
    for (size_t i = 0; i < N; i++)
    {
        size_t space = randomize_coeffs();
        if (I->F->is_valid())
        {
            cur_coeffs = I->F->get_coeffs();
            cur_epsilon = epsilon();
            cur_weight = set_weight(space, cur_coeffs, cur_epsilon);
            if (cur_epsilon < min_epsilon)
            {
                opt_coeffs = cur_coeffs;
                min_epsilon = cur_epsilon;
            }
            size_t l = N_new_spaces;
            for (int j = N_new_spaces - 1; j >= 0; j--)
            {
                if (cur_weight[N_coeffs + 1] < res[j][N_coeffs + 1])
                {
                    l = j;
                }
                else
                {
                    break;
                }
            }
            if (l != N_new_spaces)
                res[l] = cur_weight;
        }
        else
        {
            i--;
        }
    }
    return res;
}

void Optimizer::make_new_spaces(const vec2d &grids)
{
    for(auto it : grids)
    {
        
    }
}
