#include "statistics.hpp"

#include <cassert>
#include <cmath>
#include <numeric>

double CalculateMean(const std::vector<double>& y) {
	return std::accumulate(y.begin(), y.end(), 0.0) / static_cast<double>(y.size());
}

std::tuple<double, double> CalculateMeanAndVariance(const std::vector<double>& y) {
	const auto mean = CalculateMean(y);
	if (y.size() <= 1) {
		return {mean, 0.0};
	}

	const auto size_minus_1 = static_cast<double>(y.size()) - 1.0;
	const auto var = std::accumulate(y.begin(), y.end(), 0.0, [mean, size_minus_1](double acc, double cur) {
		return acc + std::pow(cur - mean, 2) / size_minus_1;
	});

	return {mean, var};
}
