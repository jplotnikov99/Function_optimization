#include "optimizer.hpp"

Optimizer::Optimizer(std::unique_ptr<Integrator> &inte, const vec1d &lower, const vec1d &upper,
                     const std::string file_name)
{
    save_file = file_name;
    for (size_t i = 0; i < N_coeffs; i++)
    {
        header.push_back("c" + std::to_string(i));
    }
    header.push_back("eps");

    I = std::move(inte);
    N_coeffs = I->F->get_N_coeffs();
    assert(N_coeffs == lower.size());
    assert(N_coeffs == upper.size());

    coefficent_spaces.push_back(lower);
    coefficent_spaces.push_back(upper);
    update_grid();
}

vec1d Optimizer::get_cur_coeffs()
{
    return I->F->get_coeffs();
}

void Optimizer::change_coeff(const size_t c_i, const double new_val)
{
    I->F->change_coeff(c_i, new_val);
}

void Optimizer::set_opt_coeffs()
{
    I->F->change_all_coeffs(opt_coeffs);
}

vec1d Optimizer::get_opt_coeffs()
{
    return opt_coeffs;
}

void Optimizer::add_space(const vec1d &lower, const vec1d &upper)
{
    coefficent_spaces.push_back(lower);
    coefficent_spaces.push_back(upper);
    N_spaces++;
}

void Optimizer::reset_space()
{
    min_epsilon = 1e100;
    opt_coeffs.clear();
    coefficent_spaces.clear();
    N_spaces = 0;
}

void Optimizer::print_space()
{
    for (auto &it : coefficent_spaces)
    {
        for (auto &jt : it)
        {
            std::cout << jt << "\t";
        }
        std::cout << "\n";
    }
}

double Optimizer::epsilon()
{
    I->switch_to_res();
    double p = I->F->get_p_value();
    return pow(I->integrate(), 1 / p);
}

vec1d Optimizer::grad_epsilon(const int coeff)
{
    vec1d res;
    double p = I->F->get_p_value();
    I->switch_to_res();
    double outer = pow(I->integrate(), 1 / p - 1);
    I->switch_to_grad();

    if (coeff != -1)
    {
        return {outer * I->integrate((size_t)coeff)};
    }

    for (size_t i = 0; i < N_coeffs; i++)
    {
        res.push_back(I->integrate(i));
    }
    return outer * res;
}

double Optimizer::get_min_epsilon()
{
    return min_epsilon;
}

void Optimizer::update_grid()
{
    /* weights data structure example for 2 bins, 3 coeffs, 2 coefficent_spaces
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
        I->F->change_coeff(i, generate_random(coefficent_spaces[2 * space][i], coefficent_spaces[2 * space + 1][i]));

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
        i2 = std::ceil((constants.at(i) - coefficent_spaces[2 * space][i]) * (double)N_bins /
                       (coefficent_spaces[2 * space + 1][i] - coefficent_spaces[2 * space][i])) -
             1;
        weights[i1 + i][i2] = weights[i1 + i][i2] == 0 ? eps : std::min(weights[i1 + i][i2], eps);
        res[i] = i2;
    }
    res[N_coeffs] = space;
    res[N_coeffs + 1] = eps;

    return res;
}

void Optimizer::make_new_spaces(const vec2d &grids)
{
    vec1d lo, up;
    vec2d new_spaces;
    double xN, x0, del;
    for (auto it : grids)
    {
        for (size_t i = 0; i < N_coeffs; i++)
        {
            xN = coefficent_spaces[2 * it[N_coeffs + 1] + 1][i];
            x0 = coefficent_spaces[2 * it[N_coeffs + 1]][i];
            del = (xN - x0) / N_bins;
            lo.push_back(x0 + it[i] * del);
            up.push_back(x0 + (it[i] + 1) * del);
        }
        new_spaces.push_back(lo);
        new_spaces.push_back(up);
        lo.clear();
        up.clear();
    }
    coefficent_spaces = new_spaces;
    N_spaces = grids.size();
    update_grid();
}

void Optimizer::monte_carlo(const size_t N, const size_t N_new_spaces)
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
    make_new_spaces(res);
}

void Optimizer::repeated_monte_carlo(const size_t N_points, const size_t N_loops, const size_t N_new_spaces)
{
    double best_eps = min_epsilon;
    vstring header;

    for (size_t i = 0; i < N_loops; i++)
    {
        monte_carlo(N_points, N_new_spaces);
        std::cout << "Iteration: " << i + 1 << ". Min epsilon value: " << get_min_epsilon() << "\n";
    }
    if (best_eps > min_epsilon)
    {
        vec1d data = opt_coeffs;
        data.push_back(min_epsilon);
        save_data(save_file, header, data);
    }
}

vec1d Optimizer::gradient_descent(const int coeff)
{
    const double RATE = 0.01, ACCURACY = 1e-4;
    double cur_eps = epsilon(), old_eps;
    const size_t MAX_IT = 10000;
    size_t CUR_IT = 0;
    vec1d cur_coeff, grad;
    do
    {
        old_eps = cur_eps;
        cur_coeff = I->F->get_coeffs();
        grad = grad_epsilon(coeff);

        if (coeff != -1)
        {
            I->F->change_coeff(coeff, cur_coeff[coeff] - RATE * grad[0]);
        }
        else
        {
            for (size_t i = 0; i < N_coeffs; i++)
            {
                I->F->change_all_coeffs(cur_coeff - RATE * grad);
            }
        }
        cur_eps = epsilon();
        CUR_IT++;
    } while ((fabs((cur_eps - old_eps) / old_eps) > ACCURACY) && (CUR_IT < MAX_IT) && (cur_eps < old_eps));
    cur_coeff.push_back(cur_eps);
    save_data(save_file, header, cur_coeff);
    return cur_coeff;
}