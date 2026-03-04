#include "D2.hpp"

void dD2dw_exact::set_x(const double new_x) { x = new_x; }

void dD2dw_exact::set_vw(const double new_vw) { vw = new_vw; }

double dD2dw_exact::operator()(const double u) {
    double pre = 6 / (M_PI * M_PI * vw * vw * vw);
    double w = x + (1 - u) / u;
    double expw = exp(w / 2);
    return pre * w *
           (vw * (2 * vw * vw - 1) * sqrt(w * w - x * x) +
            (vw * vw - 1) * (vw * vw - 1) * w *
                std::atanh(vw * sqrt(1 - x * x / (w * w)))) /
           ((expw + s / expw) * (expw + s / expw));
}