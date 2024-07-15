#include <iostream>
#include <memory>
#include <iomanip>
#include "include/functions.hpp"
#include "include/integrator.hpp"
#include "include/optimizer.hpp"

int main()
{

    srand((unsigned)time(NULL));
    std::cout << std::setprecision(3);
    std::unique_ptr<Function> F = std::make_unique<Function>(besselK1, 3);

    std::unique_ptr<Integrator> I = std::make_unique<Integrator>(F, gauss15);

    std::unique_ptr<Optimizer> O = std::make_unique<Optimizer>(I, 0, 100);

    std::vector<double> lo = {0, -2, -2, -2};
    std::vector<double> up = {2, 2, 2, 2};

    O->monte_carlo(lo, up, 10000);

    std::cout << O->get_min_epsilon() << "\n";

    std::vector<double> opt_c = O->get_opt_c();
    
    O->print_grid();

    return 0;
}