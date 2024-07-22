#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <vector>

#define STR(a) #a

extern bool first_save;

typedef std::vector<double> vec1d;
typedef std::vector<vec1d> vec2d;
typedef std::vector<std::string> vstring;

double generate_random(const double a, const double b);

vec1d operator+(const vec1d a, const vec1d b);

vec1d operator-(const vec1d a, const vec1d b);

vec1d operator*(const double a, const vec1d b);

vec1d operator/(const vec1d a, const double b);

double vabs(const vec1d &a);

void save_data(const std::string &output_file, const vstring &header, const vec1d &data);
