#include <iostream>
#include <memory>
#include "include/functions.hpp"
#include "include/integrator.hpp"
#include "include/optimizer.hpp"

int main()
{

    srand((unsigned)time(NULL));
    std::unique_ptr<Function> F = std::make_unique<Function>(besselK1, 2);

    std::unique_ptr<Integrator> I = std::make_unique<Integrator>(F, gauss15);

    std::unique_ptr<Optimizer> O = std::make_unique<Optimizer>(I, 0, 100);

    std::vector<double> lo = {0, -2, -2, -2};
    std::vector<double> up = {2, 2, 2, 2};

    O->monte_carlo(lo, up, 100000);

    std::cout << O->get_min_epsilon() << "\n";

    return 0;
}