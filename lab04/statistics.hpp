#pragma once

#include <tuple>
#include <vector>

double CalculateMean(const std::vector<double>& y);
std::tuple<double, double> CalculateMeanAndVariance(const std::vector<double>& y);
