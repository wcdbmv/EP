#pragma once

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
