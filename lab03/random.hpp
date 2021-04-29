#pragma once

#include <tuple>

double uniform_real(double a, double b);
double weibull_real(double k, double lambda);

std::tuple<double, double> uniform_mean_and_variance(double a, double b);
std::tuple<double, double> weibull_mean_and_variance(double k, double lambda);

std::tuple<double, double> uniform_parameters_from_mean_and_std(double mean, double std);
std::tuple<double, double> weibull_parameters_from_mean_and_std(double mean, double std);
std::tuple<double, double> weibull_parameters_from_mean(double mean);
