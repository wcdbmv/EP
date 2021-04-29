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

	constexpr double Norm(double x) const {
		return 2.0 * (x - min) / (max - min) - 1.0;
	}

	constexpr double Choose(double norm) const {
		return (norm + 1.0) * (max - min) / 2.0 + min;
	}
};


struct ExperimentParameters {
	Range lambda1;
	Range lambda2;
	Range mu1;
	Range mu2;
	Range sigma_lambda1;
	Range sigma_lambda2;
	size_t times;
};

struct ExperimentResult {
	std::vector<std::vector<double>> planning_matrix;
	std::vector<double> coefficients;
};

struct DotParameters {
	double lambda1;
	double lambda2;
	double mu1;
	double mu2;
	double sigma_lambda1;
	double sigma_lambda2;
};

ExperimentResult RunExperiment(const ExperimentParameters& params);

std::vector<double> NormalizeFactors(const ExperimentParameters& exp_params, const DotParameters& dot_params);
double CalculateDotWithRegression(const std::vector<double>& coefficients, const std::vector<double>& factors);
double CalculateDot(const DotParameters& dot_params, size_t times);
