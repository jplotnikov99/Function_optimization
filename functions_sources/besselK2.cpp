#include "besselK2.hpp"

double BesselK2::exact(const double x)
{
    return std::cyl_bessel_k(2, x);
}

double BesselK2::approx(const double x)
{
    return exp(-x) * (1 + c[5] / x + 2 / (x * x)) /
           pow((c[1] * pow(x, c[0] / 5) + c[2] * pow(x, c[0] / 4) + c[3] * pow(x, c[0] / 3) + c[4] * pow(x, c[0] / 2) + pow(2 * x / M_PI, c[0]) + 1), 1 / (2 * c[0]));
}

double BesselK2::gradient(const double x)
{
    vec1d res;
    vec1d temp;

    double den = c[1] * pow(x, c[0] / 5) + c[2] * pow(x, c[0] / 4) + c[3] * pow(x, c[0] / 3) +
                 c[4] * pow(x, c[0] / 2) + pow(2 * x / M_PI, c[0]) + 1;
    double e = exp(-x);
    double pre = (1 + c[5] / x + 2 / (x * x));

    double dFdc0 = 1 / (2 * c[0] * c[0]) * e * pre / pow(den, 1 / (2 * c[0])) *
                   (-c[0] * (log(x) * (12 * c[1] * pow(x, c[0] / 5) + 15 * c[2] * pow(x, c[0] / 4) + 20 * c[3] * pow(x, c[0] / 3) + 30 * c[4] * pow(x, c[0] / 2)) + 60 * pow(2 * x / M_PI, c[0]) * log(2 * x / M_PI)) / (60 * den) + log(den));
    double dFdc5 = e / x / (pow(den, 1 / (2 * c[0])));
    den = pow(den, 1 / (2 * c[0]) + 1);
    pre = -e * pre / (den * 2 * c[0]);
    double dFdc1 = pre * pow(x, c[0] / 5);
    double dFdc2 = pre * pow(x, c[0] / 4);
    double dFdc3 = pre * pow(x, c[0] / 3);
    double dFdc4 = pre * pow(x, c[0] / 2);

    res = {dFdc0, dFdc1, dFdc2, dFdc3, dFdc4, dFdc5};

    return res[cur_ci];
}

bool BesselK2::is_valid()
{
    bool passed = true;
    double den;
    if (c[0] < 0)
        return false;
    for (double x = 1e-2; x < 10; x += 1e-2)
    {
        den = c[1] * pow(x, c[0] / 5) + c[2] * pow(x, c[0] / 4) + c[3] * pow(x, c[0] / 3) +
              c[4] * pow(x, c[0] / 2) + pow(2 * x / M_PI, c[0]) + 1;
        if (den < 0)
            passed = false;
    }
    return passed;
}

double BesselK2::operator()(const double x)
{
    double ex = exact(x);
    double ap = approx(x);

    switch (ot)
    {
    case Output_type::result:
        return pow(fabs(ap / ex - 1), p_value);
        break;

    case Output_type::gradient:
    {
        double deriv = gradient(x);
        return pow(ap / ex - 1., p_value - 1) * deriv / ex;
    }
    default:
        exit(1);
    }
}