#include <iostream>
#include "include/functions.hpp"

int main(){
    Functions B1(besselK1, 2);

    std::cout << B1.res(0.00001) << std::endl;
    return 0;
}