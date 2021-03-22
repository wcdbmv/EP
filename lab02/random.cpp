#include "random.hpp"

#include <cmath>
#include <random>

double uniform_real(double a, double b) {
	static thread_local std::mt19937 generator(std::random_device{}());
	std::uniform_real_distribution distribution(a, b);
	return distribution(generator);
}

double weibull_real(double k, double lambda) {
	static thread_local std::mt19937 generator(std::random_device{}());
	std::weibull_distribution<double> distribution(k, lambda);
	return distribution(generator);
}

std::tuple<double, double> uniform_mean_and_variance(double a, double b) {
	const double mean = (a + b) / 2.0;
	const double variance = pow(b - a, 2) / 12.0;
	return {mean, variance};
}

std::tuple<double, double> weibull_mean_and_variance(double k, double lambda) {
	const double mean = lambda * tgamma(1.0 + 1.0 / k);
	const double variance = pow(lambda, 2) * tgamma(1.0 + 2.0 / k) - pow(mean, 2);
	return {mean, variance};
}

std::tuple<double, double> uniform_parameters_from_mean_and_std(double mean, double std) {
	const double a = mean - sqrt(3) * std;
	const double b = mean + sqrt(3) * std;
	return {a, b};
}

std::tuple<double, double> weibull_parameters_from_mean_and_std(double mean, double std) {
	// https://stats.stackexchange.com/questions/159452/how-can-i-recreate-a-weibull-distribution-given-mean-and-standard-deviation-and
	const double k = pow(std / mean, -1.086);
	const double lambda = mean / tgamma(1 + 1 / k);
	return {k, lambda};
}

std::tuple<double, double> weibull_parameters_from_mean(double mean) {
	// std ~=~ mean / 2;
	const double k = 2.0;
	const double lambda = mean / tgamma(1 + 1 / k);
	return {k, lambda};
}
