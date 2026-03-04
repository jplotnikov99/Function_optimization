#pragma once

#include <cassert>
#include <iostream>
#include <memory>
#include <vector>

#include "integrator.hpp"
#include "utils.hpp"

class Optimizer {
   private:
    std::unique_ptr<Integrator> I;
    std::string save_file;
    vstring header;
    double min_epsilon = 1e100;
    size_t N_coeffs;
    size_t N_bins = 10, N_spaces = 1;
    vec2d coefficent_spaces;
    vec2d weights;
    std::vector<double> opt_coeffs;

   public:
    Optimizer(std::unique_ptr<Integrator> &Inte, const size_t N,
              const vec1d &lower, const vec1d &upper,
              const std::string file_name);

    vec1d get_opt_coeffs();

    void add_space(const vec1d &lower, const vec1d &upper);
    void reset_space();
    void print_space();

    // clear weights and update the vector size of the weights for new
    // coefficent space
    void update_grid();

    void print_grid_row(const size_t i);
    void print_grid();

    // returns the space in which the random coefficients were chosen
    template <class FUNC>
    size_t randomize_coeffs(FUNC &f);

    // return the bin for each coefficient, the space they are in and their
    // epsilon value
    vec1d set_weight(const size_t space, const vec1d &constants,
                     const double eps);

    // makes new search spaces with grids that have the smallest weight
    void make_new_spaces(const vec2d &grids);

    template <class FUNC>
    double epsilon(FUNC &f);

    template <class FUNC>
    vec1d grad_epsilon(FUNC &f, const int coeff = -1);
    double get_min_epsilon();

    // returns the best N_new_spaces bins for the coefficient and their epsilon
    // value
    template <class FUNC>
    void monte_carlo(FUNC &f, const size_t N, const size_t N_new_spaces = 1);

    // does a montecarlo search with N_points and N_loops keeping the best
    // N_new_spaces each time
    template <class FUNC>
    vec1d repeated_monte_carlo(FUNC &f, const size_t N_points,
                               const size_t N_loops, const size_t N_new_spaces);

    template <class FUNC>
    vec1d gradient_descent(FUNC &f, const int coeff = -1);

    template <class FUNC>
    vec1d descent_best_direction(FUNC &f);

    ~Optimizer() {};
};

template <class FUNC>
size_t Optimizer::randomize_coeffs(FUNC &f) {
    const size_t space = std::floor(generate_random(0., (double)N_spaces));

    for (size_t i = 0; i < N_coeffs; i++)
        f.change_coeff(i, generate_random(coefficent_spaces[2 * space][i],
                                          coefficent_spaces[2 * space + 1][i]));

    return space;
}

template <class FUNC>
double Optimizer::epsilon(FUNC &f) {
    f.switch_to_res();
    double p = f.get_p_value();
    return pow(I->integrate(f), 1 / p);
}

template <class FUNC>
vec1d Optimizer::grad_epsilon(FUNC &f, const int coeff) {
    vec1d res;
    double p = f.get_p_value();
    f.switch_to_res();
    double outer = pow(I->integrate(f), 1 / p - 1);
    f.switch_to_grad();

    if (coeff != -1) {
        f.select_cur_ci((size_t)coeff);
        return {outer * I->integrate(f)};
    }
    for (size_t i = 0; i < N_coeffs; i++) {
        f.select_cur_ci((size_t)i);
        res.push_back(I->integrate(f));
    }
    return outer * res;
}

template <class FUNC>
void Optimizer::monte_carlo(FUNC &f, const size_t N,
                            const size_t N_new_spaces) {
    vec2d res(N_new_spaces, vec1d(N_coeffs + 2, 1e100));
    double cur_epsilon;
    vec1d cur_coeffs;
    vec1d cur_weight;
    for (size_t i = 0; i < N; i++) {
        size_t space = randomize_coeffs(f);
        if (f.is_valid()) {
            cur_coeffs = f.get_coeffs();
            cur_epsilon = epsilon(f);
            cur_weight = set_weight(space, cur_coeffs, cur_epsilon);
            if (cur_epsilon < min_epsilon) {
                opt_coeffs = cur_coeffs;
                min_epsilon = cur_epsilon;
            }
            size_t l = N_new_spaces;
            for (int j = N_new_spaces - 1; j >= 0; j--) {
                if (cur_weight[N_coeffs + 1] < res[j][N_coeffs + 1]) {
                    l = j;
                } else {
                    break;
                }
            }
            if (l != N_new_spaces) res[l] = cur_weight;
        } else {
            i--;
        }
    }
    make_new_spaces(res);
}

template <class FUNC>
vec1d Optimizer::repeated_monte_carlo(FUNC &f, const size_t N_points,
                                      const size_t N_loops,
                                      const size_t N_new_spaces) {
    double best_eps = min_epsilon;

    for (size_t i = 0; i < N_loops; i++) {
        monte_carlo(f, N_points, N_new_spaces);
        std::cout << "Iteration: " << i + 1
                  << ". Min epsilon value: " << get_min_epsilon() << "\n";
    }
    if (best_eps > min_epsilon) {
        vec1d data = opt_coeffs;
        data.push_back(min_epsilon);
        save_data(save_file, header, data);
    }
    return opt_coeffs;
}

template <class FUNC>
vec1d Optimizer::gradient_descent(FUNC &f, const int coeff) {
    const double RATE = 0.0001, ACCURACY = 1e-4;
    double cur_eps = epsilon(f), old_eps;
    const size_t MAX_IT = 1000;
    size_t CUR_IT = 0;
    vec1d cur_coeff, grad;
    do {
        old_eps = cur_eps;
        cur_coeff = f.get_coeffs();
        grad = grad_epsilon(f, coeff);

        if (coeff != -1) {
            f.change_coeff(coeff, cur_coeff[coeff] - RATE * grad[0]);
        } else {
            f.change_all_coeffs(cur_coeff - RATE * grad);
        }
        cur_eps = epsilon(f);
        CUR_IT++;
    } while ((fabs((cur_eps - old_eps) / old_eps) > ACCURACY) &&
             (CUR_IT < MAX_IT) && (cur_eps < old_eps));
    std::cout << "After gradient descent: " << cur_eps << "\n";
    cur_coeff.push_back(cur_eps);
    save_data(save_file, header, cur_coeff);
    return cur_coeff;
}

template <class FUNC>
vec1d Optimizer::descent_best_direction(FUNC &f) {
    vec1d cur_coeff, best_coeff;
    double cur_eps, best_eps = epsilon(f);
    for (size_t i = 0; i < N_coeffs; i++) {
        cur_coeff = gradient_descent(f, i);
        cur_eps = cur_coeff.back();
        cur_coeff.pop_back();

        if (cur_eps < best_eps) {
            best_coeff = cur_coeff;
            best_eps = cur_eps;
        }
        f.change_all_coeffs(cur_coeff);
    }
    std::cout << best_eps << "\n";
    f.change_all_coeffs(best_coeff);
    return best_coeff;
}