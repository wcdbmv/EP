#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <tuple>
#include <vector>

#include "combinatorics.hpp"


struct FfeParameters {
	double lambda_min;
	double lambda_max;
	double sigma_lambda_min;
	double sigma_lambda_max;
	double mu_min;
	double mu_max;
	size_t times;
};

struct FfeTableRow {
	int index;
	double x1;
	double x2;
	double x3;
	double y_mean;
	double y_var;
	double partial_nonlinear;
	double dpn;
	double linear;
	double dl;
};

struct PlanningMatrix3 {
	double a0, a1, a2, a3;
	double a12, a13, a23;
	double a123;

	double& operator[](int i) {
		assert(0 <= i && i < 8);
		return *(&a0 + i);
	}
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
	double adequacy_var;
	double f_test;
	PartialNonlinearCoefficients<3> coefficients;
};

struct DotResult {
	double estimated_y;
	double actual_y;
};

FfeResult FullFactorialExperiment(const FfeParameters& params);
DotResult CalculateDot(const FfeParameters& params, const PartialNonlinearCoefficients<3>& m, double lambda, double sigma_lambda, double mu);
