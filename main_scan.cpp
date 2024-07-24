#include <iostream>
#include <memory>
#include <iomanip>
#include "include/functions.hpp"
#include "include/integrator.hpp"
#include "include/optimizer.hpp"

int main()
{
    clock_t begin_time = clock();

    srand((unsigned)time(NULL));
    std::cout << std::setprecision(4);
    std::unique_ptr<Function> F = std::make_unique<Function>(besselK2, 4);
    size_t NC = F->get_N_coeffs();
    std::unique_ptr<Integrator> integrator = std::make_unique<Integrator>(F, gauss15, 0, 200);

    std::vector<double> lo1 = {1.8, 0, -2, -2, 0, 2};
    std::vector<double> up1 = {2.2, 2, 2, 2, 2, 3};
    std::unique_ptr<Optimizer> O = std::make_unique<Optimizer>(integrator, lo1, up1, "besselK2_run1.dat");

    vec1d best_generation(NC), cur_generation;
    best_generation.push_back(1e100);

    for (size_t i = 0; i < 1000; i++)
    {
        clock_t begin_time = clock();
        std::cout << "Generation: " << i << "\n";
        std::cout << "-------------------------------------\n";
        O->repeated_monte_carlo(1000, 4, 1);
        O->set_opt_coeffs();
        // cur_generation = O->gradient_descent();
        // if (cur_generation.back() < best_generation.back())
        //{
        //     best_generation = cur_generation;
        // }
        O->reset_space();
        O->add_space(lo1, up1);
        std::cout << "Computation time:\n"
                  << float(clock() - begin_time) / CLOCKS_PER_SEC << "\n";
    }

    std::cout << "Computation time:\n"
              << float(clock() - begin_time) / CLOCKS_PER_SEC << "\n";

    return 0;
}