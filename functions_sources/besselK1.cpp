#include "besselK1.hpp"

void BesselK1::switch_to_res()
{
    ot = result;
}

void BesselK1::switch_to_grad()
{
    ot = gradient;
}

void BesselK1::select_cur_ci(const size_t i)
{
    cur_ci = i;
}

double BesselK1::get_p_value()
{
    return p_value;
}

vec1d BesselK1::get_coeffs()
{
    return c;
}

size_t BesselK1::get_N_coeffs()
{
    return N_coeffs;
}

void BesselK1::change_coeff(const size_t i, const double new_val)
{
    assert(i < N_coeffs);
    c[i] = new_val;
}

void BesselK1::change_all_coeffs(const vec1d &new_vals)
{
    assert(new_vals.size() == N_coeffs);
    c = new_vals;
}

void BesselK1::print_coeffs()
{
    std::cout << "Current coefficents:\n";
    for (auto &it : c)
    {
        std::cout << it << "\t";
    }
    std::cout << "\n";
}

double BesselK1::besselK1_exact(const double x)
{
    return std::cyl_bessel_k(1, x);
}

double BesselK1::besselK1_appr(const double x)
{
    return exp(-x) * (1 + 1 / x) /
           pow((c[1] * pow(x, c[0] / 5) + c[2] * pow(x, c[0] / 4) + c[3] * pow(x, c[0] / 3) + c[4] * pow(x, c[0] / 2) + pow(2 * x / M_PI, c[0]) + 1), 1 / (2 * c[0]));
}

double BesselK1::besselK1_grad(const double x)
{
    vec1d res;
    vec1d temp;

    double den = c[1] * pow(x, c[0] / 5) + c[2] * pow(x, c[0] / 4) + c[3] * pow(x, c[0] / 3) +
                 c[4] * pow(x, c[0] / 2) + pow(2 * x / M_PI, c[0]) + 1;
    double e = exp(-x);
    double pre = (1 + 1 / x);

    double dFdc0 = 1 / (2 * c[0] * c[0]) * e * pre / pow(den, 1 / (2 * c[0])) *
                   (-c[0] * (log(x) * (12 * c[1] * pow(x, c[0] / 5) + 15 * c[2] * pow(x, c[0] / 4) + 20 * c[3] * pow(x, c[0] / 3) + 30 * c[4] * pow(x, c[0] / 2)) + 60 * pow(2 * x / M_PI, c[0]) * log(2 * x / M_PI)) / (60 * den) + log(den));
    den = pow(den, 1 / (2 * c[0]) + 1);
    pre = -e * pre / (den * 2 * c[0]);
    double dFdc1 = pre * pow(x, c[0] / 5);
    double dFdc2 = pre * pow(x, c[0] / 4);
    double dFdc3 = pre * pow(x, c[0] / 3);
    double dFdc4 = pre * pow(x, c[0] / 2);

    res = {dFdc0, dFdc1, dFdc2, dFdc3, dFdc4};
    return res[cur_ci];
}

bool BesselK1::is_valid()
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

double BesselK1::operator()(const double x)
{
    double exact = besselK1_exact(x);
    double approx = besselK1_appr(x);

    switch (ot)
    {
    case result:
        return pow(fabs(approx / exact - 1), p_value);
        break;

    case gradient:
    {
        double deriv = besselK1_grad(x);
        return pow(approx / exact - 1., p_value - 1) * deriv / exact;
    }
    default:
        exit(1);
    }
}