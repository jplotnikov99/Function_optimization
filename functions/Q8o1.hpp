#pragma once

#include <cmath>

#include "../include/integrator.hpp"
#include "function.hpp"

class Q8o1int1 {
   private:
    double x, vw, gamw, w, pwt;

   public:
    double pre;
    Q8o1int1() {};

    void set_x(const double x_in);
    void set_vw(const double vw_in);
    void set_w(const double w_in);
    double operator()(const double y);

    ~Q8o1int1() {};
};

class Q8o1int2 {
   private:
    Q8o1int1 integrand;

   public:
    Q8o1int2() {};

    void set_x(const double x_in);
    void set_vw(const double vw_in);
    double operator()(const double u);

    ~Q8o1int2() {};
};

class Q8o1 : public Function {
   private:
    Q8o1int2 integrand;
    const double tc = 1.;
    double gamw, vw;
    double cvw1, cvw2, cvw3, cvw4, cvw5;
    double regularizer;
    std::unordered_map<double, double> exact_map;

   public:
    Q8o1(const size_t p) {
        N_coeffs = 6;
        p_value = p;
        for (size_t i = 0; i < N_coeffs; i++) {
            c.push_back(0);
        }
    };

    void set_vw(const double vw_in);

    double exact(const double x);
    double approx(const double x);
    bool is_valid();
    double operator()(const double x);

    ~Q8o1() {};
};