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
    std::unique_ptr<Function> F = std::make_unique<Function>(besselK1, 3);

    std::unique_ptr<Integrator> I = std::make_unique<Integrator>(F, gauss15, 0, 100);

    std::vector<double> lo1 = {0, -2, -2, -2};
    std::vector<double> up1 = {2, 2, 2, 2};
    std::vector<double> lo2 = {1, -1, -1, -1};
    std::vector<double> up2 = {3, 3, 3, 3};

    std::unique_ptr<Optimizer> O = std::make_unique<Optimizer>(I, lo1, up1);
    O->add_space(lo2, up2);
    O->update_grid();

    vec2d res = O->monte_carlo(1000, 4);

    std::vector<double> opt_c = O->get_opt_coeffs();

    std::cout << O->get_min_epsilon() << std::endl;

    std::cout << "Computation time:\n"
              << float(clock() - begin_time) / CLOCKS_PER_SEC << "\n";

    return 0;
}