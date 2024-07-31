#include <iomanip>
#include <iostream>
#include <memory>

#include "functions/besselK2.hpp"
#include "include/integrator.hpp"
#include "include/optimizer.hpp"

int main() {
    clock_t begin_time = clock();
    srand((unsigned)time(NULL));

    BesselK2 B2(4);
    std::cout << std::setprecision(4);
    std::unique_ptr<Integrator> integrator =
        std::make_unique<Integrator>(gauss15, 0, 200);
    std::vector<double> lo1 = {0., 0., 0., 0., 0., 0.};
    std::vector<double> up1 = {0., 0., 0., 0., 0., 0.};
    std::unique_ptr<Optimizer> O = std::make_unique<Optimizer>(
        integrator, B2.get_N_coeffs(), lo1, up1, "test.dat");

    vec1d cur_generation;

    B2.change_coeff(0, 1.97556);
    B2.change_coeff(1, 0.273943);
    B2.change_coeff(2, -0.765224);
    B2.change_coeff(3, 0.727862);
    B2.change_coeff(4, 1.40001);
    B2.change_coeff(5, 2.85774);

    std::cout << O->epsilon(B2) << "\n";

    O->gradient_descent(B2);

    std::cout << O->epsilon(B2) << "\n";

    std::cout << "Computation time:\n"
              << float(clock() - begin_time) / CLOCKS_PER_SEC << "\n";

    return 0;
}