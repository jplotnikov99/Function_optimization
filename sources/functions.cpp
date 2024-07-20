#include "functions.hpp"

Function::Function(const func_name n, const double p)
{
    name = n;
    p_value = p;
    prepare();
}

void Function::prepare()
{
    size_t NC;
    switch (name)
    {
    case besselK1:
        NC = 5;
        break;

    default:
        break;
    };

    c.clear();
    for (size_t i = 0; i < NC; i++)
    {
        c.push_back(0);
    }
}

double Function::get_p_value()
{
    return p_value;
}

vec1d Function::get_coeffs()
{
    return c;
}

size_t Function::get_N_coeffs()
{
    return c.size();
}

void Function::print_coeffs()
{
    std::cout << "Current coefficents:\n";
    for (auto &it : c)
    {
        std::cout << it << "\t";
    }
    std::cout << "\n";
}

void Function::change_constant(const size_t i, const double new_val)
{
    assert(i < c.size());
    c[i] = new_val;
}

void Function::randomize_constants(const double l, const double r)
{
    for (auto &it : c)
    {
        it = generate_random(l, r);
    }
}

bool Function::is_valid()
{
    switch (name)
    {
    case besselK1:
        return besselK1_valid();
        break;

    default:
        exit(1);
        break;
    }
}

double Function::grad(const double x, const size_t c_i)
{
    assert(c_i < c.size());
    assert((p_value % 2) == 0);

    double exact{0};
    double approx{0};
    double deriv;

    switch (name)
    {
    case besselK1:
        exact = besselK1_exact(x);
        approx = besselK1_appr(x);
        deriv = besselK1_grad(x, c_i);
        break;
    case test:
        break;

    default:
        break;
    }
    double outer = pow(approx / exact - 1., p_value - 1);
    return outer * deriv / exact;
}

double Function::res(const double x, const bool is_grad, const size_t c_i)
{
    if (is_grad)
        return grad(x, c_i);

    double exact{0};
    double appr{0};

    switch (name)
    {
    case besselK1:
        exact = besselK1_exact(x);
        appr = besselK1_appr(x);
        break;
    case test:
        break;

    default:
        break;
    }
    return pow(fabs(appr / exact - 1), p_value);
}

double Function::besselK1_exact(const double x)
{
    return std::cyl_bessel_k(1, x);
}

double Function::besselK1_appr(const double x)
{
    return exp(-x) * (1 + 1 / x) /
           pow((c[1] * pow(x, c[0] / 5) + c[2] * pow(x, c[0] / 4) + c[3] * pow(x, c[0] / 3) + c[4] * pow(x, c[0] / 2) + pow(2 * x / M_PI, c[0]) + 1), 1 / (2 * c[0]));
}

bool Function::besselK1_valid()
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

double Function::besselK1_grad(const double x, const size_t c_i)
{
    vec1d temp, res;

    double den = c[1] * pow(x, c[0] / 5) + c[2] * pow(x, c[0] / 4) + c[3] * pow(x, c[0] / 3) +
                 c[4] * pow(x, c[0] / 2) + pow(2 * x / M_PI, c[0]) + 1;
    double e = exp(-x);
    double pre = (1 + 1 / x);

    double dFdc0 = 1 / (2 * c[0] * c[0]) * e * pre / pow(den, 1 / (2 * c[0])) *
                   (-c[0] * (log(x) * (12 * c[1] * pow(x, c[0] / 5) + 15 * c[2] * pow(x, c[0] / 4) + 20 * c[3] * pow(x, c[0] / 3) + 30 * c[4] * pow(x, c[0] / 2)) + 60 * pow(2 * x / M_PI, c[0]) * log(2 * x / M_PI)) / (60 * den) + log(den));
    den = pow(den, 1 / (2 * c[0]) + 1);
    pre = e * pre / (den * 2 * c[0]);
    double dFdc1 = pre * pow(x, c[0] / 5);
    double dFdc2 = pre * pow(x, c[0] / 4);
    double dFdc3 = pre * pow(x, c[0] / 3);
    double dFdc4 = pre * pow(x, c[0] / 2);

    res = {dFdc0, dFdc1, dFdc2, dFdc3, dFdc4};

    return res[c_i];
}
