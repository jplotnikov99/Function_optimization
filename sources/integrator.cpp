#include "integrator.hpp"

Integrator::Integrator(std::unique_ptr<Function> &function, const int_method m, const double x_ini, const double x_fin)
{
    F = std::move(function);
    method = m;
    xi = x_ini;
    xf = x_fin;
}

void Integrator::switch_to_grad()
{
    is_grad = true;
}

void Integrator::switch_to_res()
{
    is_grad = false;
}

double Integrator::kronrod_61(const double l, const double r, const size_t c_i)
{
    double m = 0.5 * (r + l);
    double h = 0.5 * (r - l);
    double res{0};

    for (size_t i = 0; i < 30; i++)
    {
        double dx = h * kronx_61[i];
        res += wkron_61[i] * (F->res(m + dx, is_grad, c_i) + F->res(m - dx, is_grad, c_i));
    }
    res += wkron_61[30] * F->res(m, is_grad, c_i);
    return h * res;
}

double Integrator::adap_gauss_kronrod_15(double l, double r, const double appr, const size_t c_i, const double err)
{
    double I1, I2, y[15];
    double h = (r - l) / 2;
    for (int i = 0; i < 15; i++)
    {
        y[i] = F->res((kronx_15[i] + 1) * h + l, is_grad, c_i);
    }
    I1 = h * (0.022935322010529224963732008058970 * (y[0] + y[14]) + 0.063092092629978553290700663189204 * (y[1] + y[13]) + 0.104790010322250183839876322541518 * (y[2] + y[12]) + 0.140653259715525918745189590510238 * (y[3] + y[11]) + 0.169004726639267902826583426598550 * (y[4] + y[10]) + 0.190350578064785409913256402421014 * (y[5] + y[9]) + 0.204432940075298892414161999234649 * (y[6] + y[8]) + 0.209482141084727828012999174891714 * y[7]);
    if ((I1 == 0))
    {
        return I1;
    }
    I2 = h * (0.129484966168869693270611432679082 * (y[1] + y[13]) + 0.279705391489276667901467771423780 * (y[3] + y[11]) + 0.381830050505118944950369775488975 * (y[5] + y[9]) + 0.417959183673469387755102040816327 * y[7]);

    if (fabs((I1 - I2)) < err * fabs(appr))
    {
        return I1;
    }
    double m = (l + r) / 2;
    return adap_gauss_kronrod_15(l, m, appr, c_i, err) + adap_gauss_kronrod_15(m, r, appr, c_i, err);
}

double Integrator::integrate(const size_t c_i, const double err)
{
    switch (method)
    {
    case gauss15:
    {
        double approx = kronrod_61(xi, xi + 1e-3, c_i);
        approx += kronrod_61(xi + 1e-3, xf, c_i);
        return adap_gauss_kronrod_15(xi, xf, approx, c_i, err);
    }
    break;

    default:
        exit(1);
        break;
    }
}