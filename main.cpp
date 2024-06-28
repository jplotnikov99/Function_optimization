#include <iostream>
#include <memory>
#include "include/functions.hpp"
#include "include/integrator.hpp"
#include "include/optimizer.hpp"

int main(){

    srand((unsigned) time(NULL));
    std::unique_ptr<Function> F = std::make_unique<Function>(besselK1, 2);

    F->randomize_constants(-2,2);
    std::cout << F->res(1) << std::endl;
    exit(1);
    
    std::unique_ptr<Integrator> I = std::make_unique<Integrator>(F, gauss15);

    std::unique_ptr<Optimizer> O = std::make_unique<Optimizer>(I, 0, 500);

    O->monte_carlo(-2,2,1);

    std::cout << O->get_min_epsilon() << "\n";

    return 0;
}