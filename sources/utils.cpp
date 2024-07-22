#include "utils.hpp"

double generate_random(const double a, const double b)
{
    double random = rand();
    random /= (double)RAND_MAX;
    return a + (b - a) * random;
}

vec1d operator+(const vec1d a, const vec1d b)
{
    assert(a.size() == b.size());
    vec1d res;
    for (size_t i = 0; i < a.size(); i++)
        res.push_back(a[i] + b[i]);
    return res;
}

vec1d operator-(const vec1d a, const vec1d b)
{
    assert(a.size() == b.size());
    vec1d res;
    for (size_t i = 0; i < a.size(); i++)
        res.push_back(a[i] - b[i]);
    return res;
}

vec1d operator*(const double a, const vec1d b)
{
    vec1d res;
    for (auto it : b)
        res.push_back(a * it);
    return res;
}

vec1d operator/(const vec1d a, const double b)
{
    vec1d res;
    for (auto it : a)
    {
        res.push_back(it / b);
    }
    return res;
}

double vabs(const vec1d &a)
{
    double res = 0.;
    for (auto it : a)
    {
        res += (it * it);
    }
    return sqrt(res);
}

bool first_save = true;
void save_data(const std::string &output_file, const vstring &header, const vec1d &data)
{
    std::string filesave = "../" + output_file;

    if (first_save)
    {
        std::ofstream reset;
        reset.open(filesave, std::ofstream::out | std::ofstream::trunc);
        reset.close();
    }

    std::ofstream outfile(filesave, std::ios::out | std::ios::app);
    outfile.seekp(0, std::ios::end);

    if (first_save)
    {
        for (auto &it : header)
        {
            outfile << it << "\t";
        }
        outfile << "\n";
        first_save = false;
    }

    for (auto &it : data)
    {
        outfile << it << "\t";
    }

    outfile << "\n";

    outfile.close();
}
