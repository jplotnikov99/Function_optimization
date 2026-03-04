#include "Q8o1.hpp"

void Q8o1int1::set_x(const double x_in) { x = x_in; }
void Q8o1int1::set_vw(const double vw_in) {
    vw = vw_in;
    gamw = 1 / std::sqrt(1 - vw_in * vw_in);
}
void Q8o1int1::set_w(const double u_in) {
    w = x + (1 - u_in) / u_in;
    double expw = std::exp(-w / 2.);
    pwt = std::sqrt(w * w - x * x);
    pre = -pwt / std::pow(1 / expw + expw, 2);
}

void Q8o1int2::set_vw(const double vw_in) { integrand.set_vw(vw_in); }

void Q8o1int2::set_x(const double x_in) { integrand.set_x(x_in); }

double Q8o1int1::operator()(const double y) {
    const double pzt = gamw * (y * pwt - w * vw);
    const double Et = gamw * (w - vw * y * pwt);
    const double Vx =
        1. / (1. + x * x / (pzt * pzt)) / std::sqrt(1. - pzt * pzt / (Et * Et));
    return Vx / pzt;
}

double Q8o1int2::operator()(const double u) {
    integrand.set_w(u);
    if (std::abs(integrand.pre) == 0.) return 0.;
    return integrand.pre * adap_gauss_kronrod_15(integrand, -1, 1, 1e-6) /
           (u * u);
}

void Q8o1::set_vw(const double vw_in) {
    exact_map.clear();
    integrand.set_vw(vw_in);
    vw = vw_in;
    gamw = 1 / std::sqrt(1 - vw_in * vw_in);
    double vw2 = vw_in * vw_in;
    double vw4 = vw2 * vw2;
    double vw6 = vw4 * vw2;
    double vw8 = vw4 * vw4;
    double vw10 = vw8 * vw2;
    cvw1 = -6996285. - 38505000. * vw2 - 62801664. * vw4 + 209737728. * vw6 -
           118063104. * vw8 + 27525120. * vw10;
    cvw2 = 7416744. - 5174016. * vw2 + 4884480. * vw4 - 13467648. * vw6 +
           3932160. * vw8;
    cvw3 = -2877696. + 4245504. * vw2 - 1474560. * vw4 + 786432. * vw6;
    cvw4 = 817152. - 1343488. * vw2 + 262144. * vw4;
    cvw5 = -98304. + 262144. * vw2;
    regularizer = -6.087376104758189e-6 * cvw1;
};

double Q8o1::exact(const double x) {
    double res = 0.;
    if (exact_map.count(x) == 0) {
        integrand.set_x(x);
        double res = -3 / (2. * M_PI * M_PI * gamw) * pow(tc, -2) *
                     adap_gauss_kronrod_15(integrand, 0., 1., 1e-6);
        exact_map[x] = res;
        return res;
    } else {
        return exact_map[x];
    }
}

double Q8o1::approx(const double x) {
    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;
    double x5 = x4 * x;
    return 3. * vw * sqrt(1. - vw * vw) * exp(-x) /
           (262144. * sqrt(2.) * pow(M_PI, 1.5) * tc * tc *
            sqrt(regularizer * regularizer + pow(x, 11))) *
           (cvw1 + cvw2 * x + cvw3 * x2 + cvw4 * x3 + cvw5 * x4 - 262144. * x5);
}

bool Q8o1::is_valid() {
    for (double x = 0.01; x < 10; x += 0.01) {
        double x2 = x * x;
        double x3 = x2 * x;
        double x4 = x3 * x;
        double x5 = x4 * x;
        double test = c[0] + c[1] * x + c[2] * x2 + c[3] * x3 + c[4] * x4 +
                      c[5] * x5 + pow(x, 11);
        if (test < 0) return false;
    }

    return true;
}

double Q8o1::operator()(const double x) {
    double ex = exact(x);
    double ap = approx(x);

    switch (ot) {
        case Output_type::result:
            return pow(fabs(ap / ex - 1), p_value);
            break;

        /* case Output_type::gradient: {
            double deriv = gradient(x);
            return pow(ap / ex - 1., p_value - 1) * deriv / ex;
        } */
        default:
            exit(1);
    }
}