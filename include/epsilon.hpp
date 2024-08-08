#pragma once
#include <memory>

#include "integrator.hpp"
#include "utils.hpp"

class Eps {
   private:
    std::unique_ptr<Integrator> I;

   public:
    Eps(std::unique_ptr<Integrator> &Inte) {
        I = std::move(Inte);
    };
    template <class FUNC>
    double epsilon(FUNC &f);

    template <class FUNC>
    vec1d grad_epsilon(FUNC &f, const int coeff = -1);

    ~Eps() {};
};

template <class FUNC>
double Eps::epsilon(FUNC &f) {
    f.switch_to_res();
    double p = f.get_p_value();
    return pow(I->integrate(f), 1 / p);
}

template <class FUNC>
vec1d Eps::grad_epsilon(FUNC &f, const int coeff) {
    vec1d res;
    double p = f.get_p_value();
    f.switch_to_res();
    double outer = pow(I->integrate(f), 1 / p - 1);
    f.switch_to_grad();

    if (coeff != -1) {
        f.select_cur_ci((size_t)coeff);
        return {outer * I->integrate(f)};
    }
    for (size_t i = 0; i < f.get_N_coeffs(); i++) {
        f.select_cur_ci((size_t)i);
        res.push_back(I->integrate(f));
    }
    return outer * res;
}