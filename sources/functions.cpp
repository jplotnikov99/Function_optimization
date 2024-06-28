#include "functions.hpp"

Function::Function(const func_name n, const double p)
{
    name = n;
    p_value = p;
    prepare();
}

double Function::get_p_value()
{
    return p_value;
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

void Function::prepare()
{
    size_t NC;
    switch (name)
    {
    case besselK1:
        NC = 4;
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

double Function::besselK1_exact(const double x)
{
    return std::cyl_bessel_k(1, x);
}

double Function::besselK1_appr(const double x)
{
    const double c0 = 1.984;
    return exp(-x) * (1 + 1 / x) /
           pow((c[0] * pow(x, c0 / 5) + c[1] * pow(x, c0 / 4) + c[2] * pow(x, c0 / 3) + c[3] * pow(x, c0 / 2) + pow(2 * x / M_PI, c0) + 1), 1 / (2 * c0));
}

bool Function::besselK1_check()
{
    
}

double Function::res(const double x)
{
    double exact{0};
    double appr{0};

    switch (name)
    {
    case besselK1:
        exact = besselK1_exact(x);
        appr = besselK1_appr(x);
        std::cout << appr << std::endl;
        break;
    case test:
        break;

    default:
        break;
    }
    return pow(fabs(appr / exact - 1), p_value);
}