#include "utils.hpp"

double generate_random(const double a, const double b)
{
    double random = rand();
    random /= (double)RAND_MAX;
    return a + (b - a) * random;
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
