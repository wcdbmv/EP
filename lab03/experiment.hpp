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
};


struct FfeParameters {
	Range lambda1;
	Range lambda2;
	Range mu;
	Range sigma_lambda1;
	Range sigma_lambda2;
	size_t times;
};

struct FfeTableRow {
	int index;
	double x1;
	double x2;
	double x3;
	double x4;
	double x5;
	double y_mean;
	double y_var;
	double y_hat;
	double dy_hat;
	double u_hat;
	double du_hat;
};

template <size_t k>
class PartialNonlinearCoefficients {
public:
	PartialNonlinearCoefficients()
		: N_(std::pow(2u, k))
		, a_(N_)
	{
		static_assert(1 <= k && k <= 5, "1 <= k <= 5");
	}

	size_t N() const {
		return N_;
	}

	double& operator[](size_t i) {
		return a_.at(i);
	}

	double operator[](size_t i) const {
		return a_.at(i);
	}

	template <typename... Args, typename = std::enable_if_t<std::conjunction_v<std::is_convertible<Args, size_t>...>>>
	double& a(Args... indices) {
		constexpr size_t size = sizeof...(indices);
		static_assert(size <= k, "size > k");

		const std::array<size_t, size> indices_array{{static_cast<size_t>(indices)...}};

		const size_t index = CalculateIndexStartOf(size) + CalculateIndexOfCombination(k, indices_array) - 1;

		return a_[index];
	}

private:
	size_t N_;
	std::vector<double> a_;

	size_t CalculateIndexStartOf(size_t together) {
		return together ? CalculateIndexStartOf(together - 1) + CalculateBinomialCoefficient(k, together - 1) : 0;
	}
};

using FfeTable = std::vector<FfeTableRow>;

struct FfeResult {
	FfeTable table;
	double cochran_test;
	double reproducibility_var;
	double y_hat_adequacy_var;
	double y_hat_f_test;
	double u_hat_adequacy_var;
	double u_hat_f_test;
	PartialNonlinearCoefficients<5> coefficients;
};

struct DotParameters {
	double lambda1;
	double lambda2;
	double mu;
	double sigma_lambda1;
	double sigma_lambda2;
};

struct DotResult {
	double estimated_y;
	double actual_y;
};

FfeResult FullFactorialExperiment(const FfeParameters& params);
DotResult CalculateDot(const FfeParameters& ffe_params, const PartialNonlinearCoefficients<5>& coefficients, const DotParameters& dot_params);
