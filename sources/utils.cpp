#include "utils.hpp"

double generate_random(const double a, const double b)
{
    double random = rand();
    random /= (double)RAND_MAX;
    return a + (b - a) * random;
}