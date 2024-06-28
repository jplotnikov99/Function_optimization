#include <iostream>
#include <memory>
#include "include/functions.hpp"
#include "include/integrator.hpp"

int main(){
    std::unique_ptr<Function> B1 = std::make_unique<Function>(besselK1, 2);
    
    Integrator I(B1, gauss15);

    std::cout << I.adap_gauss_kronrod_15(700,2) << std::endl;
    return 0;
}