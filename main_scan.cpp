#include <iomanip>
#include <iostream>
#include <memory>

#include "functions/D0.hpp"
#include "functions/Q8o1.hpp"
#include "include/integrator.hpp"
#include "include/optimizer.hpp"

int main() {
    clock_t begin_time = clock();

    srand((unsigned)time(NULL));
    std::cout << std::setprecision(16);
    Q8o1 q8o1(2);

    std::ofstream yfile("../Q8o1_test.dat");
    double N = 200;
    double lo = log10(1e-3), hi = log10(20);
    double step = (hi - lo) / N;
    q8o1.set_vw(0.98);
    for (double i = 0; i <= N; i++) {
        double x = pow(10, lo + i * step);
        yfile << x << "\t" << q8o1.approx(x) << "\t" << q8o1.exact(x) << "\n";
    }
    yfile.close();

    exit(1);

    std::unique_ptr<Integrator> I2 =
        std::make_unique<Integrator>(gauss15, 0., 150.);
    const double a = 3;
    std::vector<double> lo1 = {0, 0, -a, -a, -a, -a};
    std::vector<double> up1 = {a, a, a, a, a, a};
    std::unique_ptr<Optimizer> O = std::make_unique<Optimizer>(
        I2, q8o1.get_N_coeffs(), lo1, up1, "Q8o1_test.dat");
    vec1d best_generation(q8o1.get_N_coeffs()), cur_generation;

    best_generation.push_back(1e100);

    for (size_t i = 0; i < 1000; i++) {
        std::cout << "Generation: " << i << "\n";
        std::cout << "-------------------------------------\n";
        cur_generation = O->repeated_monte_carlo(q8o1, 1000, 3, 8);
        q8o1.change_all_coeffs(cur_generation);
        // cur_generation = O->gradient_descent(d0);
        if (cur_generation.back() < best_generation.back()) {
            best_generation = cur_generation;
        }
        O->reset_space();
        O->add_space(lo1, up1);
    }
    exit(1);

    /* D0 d0(2, boson);
    std::unique_ptr<Integrator> I2 =
        std::make_unique<Integrator>(gauss15, 0, 150.);
    const double a = 3;
    std::vector<double> lo1 = {-a, -a, 0, -a, 0, -a, 0, 0, -a, 0, 0};
    std::vector<double> up1 = {a, a, a, a, a, a, a, a, a, a, a};
    std::unique_ptr<Optimizer> O = std::make_unique<Optimizer>(
        I2, d0.get_N_coeffs(), lo1, up1, "D0_test.dat");
    vec1d best_generation(d0.get_N_coeffs()), cur_generation;
    best_generation.push_back(1e100);

    for (size_t i = 0; i < 1000; i++) {
        std::cout << "Generation: " << i << "\n";
        std::cout << "-------------------------------------\n";
        cur_generation = O->repeated_monte_carlo(d0, 1000, 3, 8);
        d0.change_all_coeffs(cur_generation);
        cur_generation = O->gradient_descent(d0);
        if (cur_generation.back() < best_generation.back()) {
            best_generation = cur_generation;
        }
        O->reset_space();
        O->add_space(lo1, up1);
    }
    best_generation.pop_back();
    d0.change_all_coeffs(best_generation);
    O->gradient_descent(d0); */

    std::cout << "Computation time:\n"
              << float(clock() - begin_time) / CLOCKS_PER_SEC << "\n";

    return 0;
}