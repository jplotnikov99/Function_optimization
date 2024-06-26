#include <iostream>

class BesselK1
{
private:
    double c[10];
    const size_t N = 10;
public:
    BesselK1(){};
    double operator()(double x);
    ~BesselK1(){};
};