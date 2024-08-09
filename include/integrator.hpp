#pragma once
#include <iostream>

#include "abscissa.hpp"

enum int_method { gauss15, simps38 };

class Integrator {
   private:
    bool debug = false;
    const double x_lower, x_upper;
    const int_method method;

    template <class FUNC>
    double kronrod_61(FUNC &f, const double l, const double r);

    template <class FUNC>
    double adap_gauss_kronrod_15(FUNC &f, const double l, const double r,
                                 const double appr, const double err = 1e-10,
                                 const size_t depth = 0);

    double simpson_est(const double l, const double r, const double *f) {
        return (r - l) / 24 *
               (f[0] + f[9] + 3 * (f[1] + f[2] + f[4] + f[5] + f[7] + f[8]) +
                2 * (f[3] * f[6]));
    }

    template <class FUNC>
    double adap_simpson38(FUNC &f, const double l, const double r, double *f0,
                          const double appr, const double err = 1e-10,
                          const size_t depth = 0);

   public:
    Integrator(const int_method m, const double xi, const double xu)
        : method(m), x_lower(xi), x_upper(xu) {};

    void switch_debug() { debug = !debug; };

    template <class FUNC>
    double integrate(FUNC &f, const double err = 1e-10);

    ~Integrator() {};
};

template <class FUNC>
double Integrator::kronrod_61(FUNC &f, const double l, const double r) {
    double m = 0.5 * (r + l);
    double h = 0.5 * (r - l);
    double res{0};

    for (size_t i = 0; i < 30; i++) {
        double dx = h * kronx_61[i];
        res += wkron_61[i] * (f(m + dx) + f(m - dx));
    }
    res += wkron_61[30] * f(m);
    return h * res;
}

template <class FUNC>
double Integrator::adap_gauss_kronrod_15(FUNC &f, const double l,
                                         const double r, const double appr,
                                         const double err, const size_t depth) {
    double I1, I2, y[15];
    double h = (r - l) / 2;
    for (int i = 0; i < 15; i++) {
        y[i] = f((kronx_15[i] + 1) * h + l);
    }
    I1 = h * (0.022935322010529224963732008058970 * (y[0] + y[14]) +
              0.063092092629978553290700663189204 * (y[1] + y[13]) +
              0.104790010322250183839876322541518 * (y[2] + y[12]) +
              0.140653259715525918745189590510238 * (y[3] + y[11]) +
              0.169004726639267902826583426598550 * (y[4] + y[10]) +
              0.190350578064785409913256402421014 * (y[5] + y[9]) +
              0.204432940075298892414161999234649 * (y[6] + y[8]) +
              0.209482141084727828012999174891714 * y[7]);
    if ((I1 == 0) || (depth > 16)) {
        return I1;
    }
    I2 = h * (0.129484966168869693270611432679082 * (y[1] + y[13]) +
              0.279705391489276667901467771423780 * (y[3] + y[11]) +
              0.381830050505118944950369775488975 * (y[5] + y[9]) +
              0.417959183673469387755102040816327 * y[7]);

    if (std::abs((I1 - I2)) < err * std::abs(appr)) {
        return I1;
    }
    double m = (l + r) / 2;
    return adap_gauss_kronrod_15(f, l, m, appr, err, depth + 1) +
           adap_gauss_kronrod_15(f, m, r, appr, err, depth + 1);
}

template <class FUNC>
double Integrator::adap_simpson38(FUNC &f, const double l, const double r,
                                  double *f0, const double appr,
                                  const double err, const size_t depth) {
    double I1, I2, f1[4];
    double m = (r + l) / 2.;
    double h = (r - l) / 8.;
    double Ia = h * (f0[0] + 3 * f0[1] + 3 * f0[2] + f0[3]);
    f1[0] = f(m);
    f1[1] = f0[2];
    f1[2] = f((l + 5 * r) / 6);
    f1[3] = f0[3];
    f0[3] = f1[0];
    f0[2] = f0[1];
    f0[1] = f((5 * l + r) / 6);
    I1 = h / 2 * (f0[0] + 3 * f0[1] + 3 * f0[2] + f0[3]);
    I2 = h / 2 * (f1[0] + 3 * f1[1] + 3 * f1[2] + f1[3]);
    double Ib = I1 + I2;

    if ((std::abs(Ia - Ib) < err * std::abs(appr)) || (depth > 16)) {
        return Ib;
    }
    return adap_simpson38(f, l, m, f0, appr, err, depth + 1) +
           adap_simpson38(f, m, r, f1, appr, err, depth + 1);
}

template <class FUNC>
double Integrator::integrate(FUNC &f, const double err) {
    switch (method) {
        case gauss15: {
            double h = x_upper - x_lower;
            double approx = kronrod_61(f, x_lower, x_lower + h * 1e-3);
            approx += kronrod_61(f, x_lower + h * 1e-3, x_upper);
            return adap_gauss_kronrod_15(f, x_lower, x_upper, approx, err);
        } break;
        case simps38: {
            double f_approx[10];
            double dx = (x_upper - x_lower) / 9.;
            for (size_t i = 0; i < 10; i++) {
                f_approx[i] = f(x_lower + dx * (double)i);
            }
            double approx = simpson_est(x_lower, x_upper, f_approx);
            double f0[4] = {f_approx[0], f_approx[3], f_approx[6], f_approx[9]};
            return adap_simpson38(f, x_lower, x_upper, f0, approx, err);
        }

        default:
            exit(1);
            break;
    }
}