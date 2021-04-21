#include "experiment.hpp"

#include <cassert>
#include <functional>
#include <numeric>

#include "random.hpp"
#include "simulate.hpp"
#include "statistics.hpp"


std::vector<std::vector<int>> CreateUnitPlanningMatrixCore(size_t k) {
	const auto N = static_cast<size_t>(std::pow(2ul, k));

	std::vector<std::vector<int>> ones;
	ones.reserve(N);
	for (size_t i = 0; i < N; ++i) {
		ones.push_back(std::vector<int>{});
		ones[i].reserve(k);
		for (size_t j = 0; j < k; ++j) {
			const size_t bit = 1u << j;
			ones[i].push_back(i & bit ? +1 : -1);
		}
	}

	return ones;
}

template <typename T>
T CalculateProductOfCombinationInRow(const std::vector<size_t>& combination, const std::vector<T>& row) {
	return std::accumulate(combination.begin(), combination.end(), T{1}, [&row](auto&& acc, auto&& cur) {
		return acc * row[static_cast<size_t>(cur)];
	});
}

void ProceedUnitPlanningMatrixCore(size_t k, std::vector<std::vector<int>>& ones) {
	const auto N = static_cast<size_t>(std::pow(2ul, k));
	const auto M = N - 1;

	for (size_t n = 0; n < ones.size(); ++n) {
		ones[n].reserve(M);

		for (size_t i = 2; i <= k; ++i) {
			for (auto&& combination : GenerateAllCombinations(k, i)) {
				ones[n].push_back(CalculateProductOfCombinationInRow(combination, ones[n]));
			}
		}
	}
}

std::vector<std::vector<int>> CreateUnitPlanningMatrix(size_t k) {
	auto ones = CreateUnitPlanningMatrixCore(k);
	ProceedUnitPlanningMatrixCore(k, ones);
	return ones;
}


std::vector<std::vector<int>> CreateUnitPlanningMatrixCore2() {
	auto ones = CreateUnitPlanningMatrixCore(5);
	for (auto&& row : ones) {
		row.push_back(row[0] * row[1] * row[2] * row[3] * row[4]);
	}
	return ones;
}

std::vector<std::vector<int>> CreateUnitPlanningMatrix2() {
	auto ones = CreateUnitPlanningMatrixCore2();
	ProceedUnitPlanningMatrixCore(6, ones);
	return ones;
}


template <size_t k>
PartialNonlinearCoefficients<k> CalculateCoefficients(const FfeTable& table) {
	PartialNonlinearCoefficients<k> coefficients{};
	const size_t N = std::min(coefficients.N(), table.size());

	const auto ones = N == coefficients.N() ? CreateUnitPlanningMatrix(k) : CreateUnitPlanningMatrix2();
	for (size_t i = 0; i < N; ++i) {
		const auto y = table[i].y_mean;
		coefficients[0] += y;
		for (size_t j = 0; j + 1 < coefficients.N(); ++j) {
			coefficients[j + 1] += ones[i][j] * y;
		}
	}
	for (size_t i = 0; i < coefficients.N(); ++i) {
		coefficients[i] /= static_cast<double>(N);
	}

	return coefficients;
}

template <size_t k>
double CalculatePartialNonlinearRegression(const PartialNonlinearCoefficients<k>& coefficients, const std::vector<double>& factors, size_t limit) {
	assert(factors.size() == k);
	assert(limit <= k);

	double y = 0.0;
	for (size_t i = 0, j = 0; i <= limit; ++i) {
		for (auto&& combination : GenerateAllCombinations(k, i)) {
			y += coefficients[j++] * CalculateProductOfCombinationInRow(combination, factors);
		}
	}

	return y;
}

template <size_t k>
using TCalculateRegression = std::function<double(const PartialNonlinearCoefficients<k>& coefficients, const std::vector<double>& factors)>;

template <size_t k>
std::vector<double> CalculateYUsingRegression(const PartialNonlinearCoefficients<k>& coefficients, TCalculateRegression<k>&& calculate_regression) {
	const auto ones = CreateUnitPlanningMatrixCore(k);

	std::vector<double> y_hat;
	for (size_t i = 0; i < coefficients.N(); ++i) {
		assert(ones[i].size() == k);

		std::vector<double> factors;
		factors.reserve(k);
		for (int one : ones[i]) {
			factors.push_back(static_cast<double>(one));
		}

		y_hat.push_back(calculate_regression(coefficients, factors));
	}
	return y_hat;
}

template <size_t k>
std::vector<double> CalculateYWithPartialNonlinearRegression(const PartialNonlinearCoefficients<k>& coefficients, size_t limit = k) {
	return CalculateYUsingRegression<6>(coefficients, [limit](const auto& c, const auto& f) { return CalculatePartialNonlinearRegression(c, f, limit); });
}

template <size_t k>
std::vector<double> CalculateYWithLinearRegression(const PartialNonlinearCoefficients<k>& coefficients) {
	return CalculateYUsingRegression<k>(coefficients, [](const auto& c, const auto& f) { return CalculatePartialNonlinearRegression(c, f, 1); });
}

std::vector<double> SimulateNTimes(const SimulateParams& params, size_t times) {
	std::vector<double> y;
	for (size_t i = 0; i < times; ++i) {
		const auto result = Simulate(params);
		y.push_back(result.average_waiting);
	}
	return y;
}

std::vector<double> CalculateY(const DotParameters& dot_params, size_t times) {
	const auto [a1, b1] = uniform_parameters_from_mean_and_std(1.0 / dot_params.lambda1, dot_params.sigma_lambda1);
	const auto [a2, b2] = uniform_parameters_from_mean_and_std(1.0 / dot_params.lambda2, dot_params.sigma_lambda2);
	const auto [k1, l1] = weibull_parameters_from_mean(1.0 / dot_params.mu1);
	const auto [k2, l2] = weibull_parameters_from_mean(1.0 / dot_params.mu2);

	const auto y = SimulateNTimes({
		.a1 = a1,
		.b1 = b1,

		.a2 = a2,
		.b2 = b2,

		.k1 = k1,
		.lambda1 = l1,

		.k2 = k2,
		.lambda2 = l2,

		.t = 2000.0,
	}, times);

	return y;
}

FfeResult FullFactorialExperiment(const FfeParameters& params) {
	FfeResult result;
	for (size_t i = 0; i < 64u; ++i) {
		const DotParameters dot_params = {
			.lambda1       = params.lambda1      .Choose(i & (0b1 << 0)),
			.lambda2       = params.lambda2      .Choose(i & (0b1 << 1)),
			.mu1           = params.mu1          .Choose(i & (0b1 << 2)),
			.mu2           = params.mu2          .Choose(i & (0b1 << 3)),
			.sigma_lambda1 = params.sigma_lambda1.Choose(i & (0b1 << 4)),
			.sigma_lambda2 = params.sigma_lambda2.Choose(i & (0b1 << 5)),
		};

		const auto y = CalculateY(dot_params, params.times);
		const auto [y_mean, y_var] = CalculateMeanAndVariance(y);

		result.table.push_back({
			.index = i + 1,
			.x1 = dot_params.lambda1,
			.x2 = dot_params.lambda2,
			.x3 = dot_params.mu1,
			.x4 = dot_params.mu2,
			.x5 = dot_params.sigma_lambda1,
			.x6 = dot_params.sigma_lambda2,

			.y_mean = y_mean,
			.y_var  = y_var,

			.y_hat  = 0.0,
			.dy_hat = 0.0,
			.u_hat  = 0.0,
			.du_hat = 0.0,
		});
	}

	result.coefficients = CalculateCoefficients<6>(result.table);
	const auto y_hat = CalculateYWithLinearRegression(result.coefficients);
	const auto u_hat = CalculateYWithPartialNonlinearRegression(result.coefficients, 6);
	for (size_t i = 0; i < result.coefficients.N(); ++i) {
		result.table[i].y_hat  = y_hat[i];
		result.table[i].dy_hat = std::abs(result.table[i].y_mean - y_hat[i]);
		result.table[i].u_hat  = u_hat[i];
		result.table[i].du_hat = std::abs(result.table[i].y_mean - u_hat[i]);
	}
	return result;
}

FfeResult FractionalFactorialExperiment(const FfeParameters& params) {
	const auto xs = CreateUnitPlanningMatrixCore2();

	FfeResult result;
	for (size_t i = 0; i < xs.size(); ++i) {
		const DotParameters dot_params = {
			.lambda1       = params.lambda1      .Choose(xs[i][0] == 1),
			.lambda2       = params.lambda2      .Choose(xs[i][1] == 1),
			.mu1           = params.mu1          .Choose(xs[i][2] == 1),
			.mu2           = params.mu1          .Choose(xs[i][3] == 1),
			.sigma_lambda1 = params.sigma_lambda1.Choose(xs[i][4] == 1),
			.sigma_lambda2 = params.sigma_lambda2.Choose(xs[i][5] == 1), // x1*x2*x3*x4*x5
		};

		const auto y = CalculateY(dot_params, params.times);
		const auto [y_mean, y_var] = CalculateMeanAndVariance(y);

		result.table.push_back({
			.index = i + 1,

			.x1 = dot_params.lambda1,
			.x2 = dot_params.lambda2,
			.x3 = dot_params.mu1,
			.x4 = dot_params.mu2,
			.x5 = dot_params.sigma_lambda1,
			.x6 = dot_params.sigma_lambda2,

			.y_mean = y_mean,
			.y_var  = y_var,

			.y_hat  = 0.0,
			.dy_hat = 0.0,
			.u_hat  = 0.0,
			.du_hat = 0.0,
		});
	}

	result.coefficients = CalculateCoefficients<6>(result.table);
	const auto y_hat = CalculateYWithLinearRegression(result.coefficients);
	const auto u_hat = CalculateYWithPartialNonlinearRegression(result.coefficients, 2);
	for (size_t i = 0; i < xs.size(); ++i) {
		result.table[i].y_hat  = y_hat[i];
		result.table[i].dy_hat = std::abs(result.table[i].y_mean - y_hat[i]);
		result.table[i].u_hat  = u_hat[i];
		result.table[i].du_hat = std::abs(result.table[i].y_mean - u_hat[i]);
	}
	return result;
}

std::vector<double> NormalizeFactors(const FfeParameters& ffe_params, const DotParameters& dot_params) {
	return {
		ffe_params.lambda1.Norm(dot_params.lambda1),
		ffe_params.lambda2.Norm(dot_params.lambda2),
		ffe_params.mu1.Norm(dot_params.mu1),
		ffe_params.mu2.Norm(dot_params.mu2),
		ffe_params.sigma_lambda1.Norm(dot_params.sigma_lambda1),
		ffe_params.sigma_lambda2.Norm(dot_params.sigma_lambda2),
	};
}

double CalculateDotWithRegression(const PartialNonlinearCoefficients<6>& coefficients, const std::vector<double>& factors, size_t limit) {
	return CalculatePartialNonlinearRegression(coefficients, factors, limit);
}

double CalculateDot(const DotParameters& dot_params, size_t times) {
	return CalculateMean(CalculateY(dot_params, times));
}
