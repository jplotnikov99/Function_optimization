#include "D0.hpp"

void dD0dw::set_x(const double new_x)
{
    x = new_x;
}

double dD0dw::operator()(const double u)
{
    const double pre = 6 / (M_PI * M_PI);
    double w = x + (1 - u) / u;
    double expw = exp(w / 2);
    return 1 / (u * u) * pre * w * sqrt(w * w - x * x) / ((expw + s / expw) * (expw + s / expw));
}

double D0::exact(const double x)
{
    double res = 0.;
    if (exact_map.count(x) == 0)
    {
        if (x > 100)
        {
            res = 6 * x * x / (M_PI * M_PI) * std::cyl_bessel_k(2, x);
        }
        else
        {

            integrand.set_x(x);
            res += I->integrate(integrand, 1e-6);
        }
        exact_map[x] = res;
        return res;
    }
    else
    {
        return exact_map[x];
    }
}

double D0::approx(const double x)
{
    return exp(-x) * (a1 - a2 + (c[0] * pow(x, c[9]) + c[8] * pow(x, c[10])) / (c[1] * pow(x, c[2]) + c[3] * pow(x, c[4]) + c[5] * pow(x, c[6]) + c[7])) +
           6 * x * x / (M_PI * M_PI) * std::cyl_bessel_k(2, x);
}

double D0::gradient(const double x)
{
    vec1d grad(11);
    double den = c[1] * pow(x, c[2]) + c[3] * pow(x, c[4]) + c[5] * pow(x, c[6]) + c[7];
    double num = (c[8] * pow(x, c[10]) + c[0] * pow(x, c[9]));
    double expx = exp(-x);
    double logx = log(x);
    grad[0] = expx * pow(x, c[9]) / den;
    grad[1] = -expx * pow(x, c[2]) * num / (den * den);
    grad[2] = grad[1] * c[1] * logx;
    grad[3] = -expx * pow(x, c[4]) * num / (den * den);
    grad[4] = grad[3] * c[3] * logx;
    grad[5] = -expx * pow(x, c[6]) * num / (den * den);
    grad[6] = grad[5] * c[5] * logx;
    grad[7] = -expx * num / (den * den);
    grad[8] = expx * pow(x, c[10]) / den;
    grad[9] = grad[0] * c[0] * logx;
    grad[10] = grad[8] * c[8] * logx;
    return grad[cur_ci];
}

bool D0::is_valid()
{
    return true;
    double a = std::max(c[2], c[4]);
    double b = std::max(c[4], c[6]);
    double highest_power = std::max(a, b);
    if ((highest_power == c[2]) && (c[1] > 0))
    {
        return true;
    }
    if ((highest_power == c[4]) && (c[3] > 0))
    {
        return true;
    }
    if ((highest_power == c[6]) && (c[7] > 0))
    {
        return true;
    }
    return false;
}

double D0::operator()(const double x)
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