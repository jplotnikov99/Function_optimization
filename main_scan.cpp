#include <iostream>
#include <memory>
#include <iomanip>
#include "functions/besselK1.hpp"
#include "include/integrator.hpp"
#include "include/optimizer.hpp"

int main()
{
    clock_t begin_time = clock();

    // srand((unsigned)time(NULL));
    std::cout << std::setprecision(4);
    BesselK1 B1(4);
    std::unique_ptr<Integrator> I2 = std::make_unique<Integrator>(gauss15, 0., 150.);
    std::vector<double> lo1 = {1.8, 0, -2, -2, 0};
    std::vector<double> up1 = {2.2, 2, 2, 2, 2};
    std::unique_ptr<Optimizer> O2 = std::make_unique<Optimizer>(I2, B1.get_N_coeffs(), lo1, up1, "test.dat");
    vec1d best_generation(B1.get_N_coeffs()), cur_generation;
    best_generation.push_back(1e100);

    for (size_t i = 0; i < 10; i++)
    {
        std::cout << "Generation: " << i << "\n";
        std::cout << "-------------------------------------\n";
        cur_generation = O2->repeated_monte_carlo(B1, 10000, 2, 1);
        B1.change_all_coeffs(cur_generation);
        cur_generation = O2->gradient_descent(B1);
        if (cur_generation.back() < best_generation.back())
        {
            best_generation = cur_generation;
        }
        O2->reset_space();
        O2->add_space(lo1, up1);
    }

    std::cout << "Computation time:\n"
              << float(clock() - begin_time) / CLOCKS_PER_SEC << "\n";

    return 0;
}