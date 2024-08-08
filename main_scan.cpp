#include <iomanip>
#include <iostream>
#include <memory>

#include "functions/D0.hpp"
#include "include/epsilon.hpp"
#include "include/integrator.hpp"
#include "include/optimizer.hpp"

int main() {
    clock_t begin_time = clock();

    // srand((unsigned)time(NULL));
    std::cout << std::setprecision(16);
    D0 d0(2, boson);
    std::unique_ptr<Integrator> I2 =
        std::make_unique<Integrator>(gauss15, 0., 200.);
    I2->switch_debug();
    std::unique_ptr<Eps> E = std::make_unique<Eps>(I2);

    const double a = 3;
    std::vector<double> lo1 = {-a, -a, 0, -a, 0, -a, 0, 0, -a, 0, 0};
    std::vector<double> up1 = {a, a, a, a, a, a, a, a, a, a, a};
    std::unique_ptr<Optimizer> O = std::make_unique<Optimizer>(
        E, d0.get_N_coeffs(), lo1, up1, "D0_test.dat");
    vec1d best_generation(d0.get_N_coeffs()), cur_generation;
    best_generation.push_back(1e100);

    for (size_t i = 0; i < 1000; i++) {
        std::cout << "Generation: " << i << "\n";
        std::cout << "-------------------------------------\n";
        cur_generation = O->repeated_monte_carlo(d0, 1000, 3, 2);
        d0.change_all_coeffs(cur_generation);
        cur_generation = O->gradient_descent(d0);
        if (cur_generation.back() < best_generation.back()) {
            best_generation = cur_generation;
        }
        O->reset_space();
        O->add_space(lo1, up1);
    }
    best_generation.pop_back();
    d0.change_all_coeffs(best_generation);
    O->gradient_descent(d0);

    std::cout << "Computation time:\n"
              << float(clock() - begin_time) / CLOCKS_PER_SEC << "\n";

    return 0;
}