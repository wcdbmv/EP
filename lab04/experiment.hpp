#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <tuple>
#include <vector>

#include "combinatorics.hpp"


struct Range {
	double min;
	double max;

	constexpr double Choose(bool c) const {
		return c ? max : min;
	}

	constexpr double Choose(double norm) const {
		return (norm + 1.0) * (max - min) / 2.0 + min;
	}

	constexpr double Norm(double x) const {
		return 2.0 * (x - min) / (max - min) - 1.0;
	}
};


struct OccdParameters {
	Range lambda1;
	Range lambda2;
	Range mu1;
	Range mu2;
	Range sigma_lambda1;
	Range sigma_lambda2;
	size_t times;
};

struct OccdTableRow {
	size_t index;
	double x1;
	double x2;
	double x3;
	double x4;
	double x5;
	double x6;
	double x12mS;
	double x22mS;
	double x32mS;
	double x42mS;
	double x52mS;
	double x62mS;
	double y_mean;
	double y_var;
	double y_hat;
	double dy_hat;
};

template <size_t k>
class NonlinearCoefficients {
public:
	NonlinearCoefficients()
		: M_(2 * k + k * (k - 1) / 2 + 1)
		, a_(M_)
	{
	}

	[[nodiscard]] size_t M() const {
		return M_;
	}

	[[nodiscard]] double& operator[](size_t i) {
		return a_.at(i);
	}

	[[nodiscard]] double operator[](size_t i) const {
		return a_.at(i);
	}

	[[nodiscard]] double a(size_t i, size_t j) const {
		return i == j ? aii(i) : aij(i, j);
	}

private:
	size_t M_;
	std::vector<double> a_;

	double aij(size_t i, size_t j) const {
		const std::array<size_t, 2> indices_array{i, j};
		const size_t index = 1 + k + CalculateIndexOfCombination(k, indices_array) - 1;
		return a_.at(index);
	}

	double aii(size_t i) const {
		const size_t index = 1 + k + 2 * k + i;
		return a_.at(index);
	}
};

using OccdTable = std::vector<OccdTableRow>;

struct OccdResult {
	OccdTable table;
	NonlinearCoefficients<6> coefficients;
};

struct DotParameters {
	double lambda1;
	double lambda2;
	double mu1;
	double mu2;
	double sigma_lambda1;
	double sigma_lambda2;
};

OccdResult OrthogonalCentralCompositeDesign(const OccdParameters& params);

std::vector<double> NormalizeFactors(const OccdParameters& ffe_params, const DotParameters& dot_params);
double CalculateDotWithRegression(const NonlinearCoefficients<6>& coefficients, const std::vector<double>& factors);
double CalculateDot(const DotParameters& dot_params, size_t times);
