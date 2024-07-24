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
    std::vector<double> lo1 = {0., 0., 0., 0., 0., 0.};
    std::vector<double> up1 = {0., 0., 0., 0., 0., 0.};
    std::unique_ptr<Optimizer> O = std::make_unique<Optimizer>(integrator, lo1, up1, "besselK2_run1.dat");

    vec1d cur_generation;

    O->change_coeff(0, 1.97556);
    O->change_coeff(1, 0.273943);
    O->change_coeff(2, -0.765224);
    O->change_coeff(3, 0.727862);
    O->change_coeff(4, 1.40001);
    O->change_coeff(5, 2.85774);

    std::cout << O->epsilon() << "\n";

    O->gradient_descent();

    std::cout << O->epsilon() << "\n";

    std::cout << "Computation time:\n"
              << float(clock() - begin_time) / CLOCKS_PER_SEC << "\n";

    return 0;
}