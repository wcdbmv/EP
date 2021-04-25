#include "experiment.hpp"

#include <cassert>
#include <functional>
#include <numeric>

#include "random.hpp"
#include "simulate.hpp"
#include "statistics.hpp"


std::vector<std::vector<double>> CreateUnitPlanningMatrixCore(size_t k) {
	const auto N0 = static_cast<size_t>(std::pow(2ul, k));
	const auto N  = N0 + 2 * k + 1;

	std::vector<std::vector<double>> ones;
	ones.reserve(N);
	for (size_t i = 0; i < N0; ++i) {
		ones.push_back(std::vector<double>{});
		ones[i].reserve(k);
		for (size_t j = 0; j < k; ++j) {
			const size_t bit = 1u << j;
			ones[i].push_back(i & bit ? +1.0 : -1.0);
		}
	}

	const auto S = std::sqrt(static_cast<double>(k) / static_cast<double>(N));
	const auto alpha = std::sqrt(static_cast<double>(k) / 2.0 * (1.0 / S - 1.0));

	for (size_t i = 0; i < k; ++i) {
		const auto i_row = N0 + 2 * i;
		ones.emplace_back(k);
		ones.emplace_back(k);
		ones[i_row][i] = alpha;
		ones[i_row + 1][i] = -alpha;
	}
	ones.emplace_back(k);

	return ones;
}

template <typename T>
T CalculateProductOfCombinationInRow(const std::vector<size_t>& combination, const std::vector<T>& row) {
	return std::accumulate(combination.begin(), combination.end(), T{1}, [&row](auto&& acc, auto&& cur) {
		return acc * row[static_cast<size_t>(cur)];
	});
}

void ProceedUnitPlanningMatrixCore(size_t k, std::vector<std::vector<double>>& ones) {
	const auto N0 = static_cast<size_t>(std::pow(2ul, k));
	const auto N = N0 + 2 * k + 1;
	const auto M = 2 * k + k * (k - 1) / 2;
	const auto S = std::sqrt(static_cast<double>(k) / static_cast<double>(N));

	for (size_t n = 0; n < ones.size(); ++n) {
		ones[n].reserve(M);

		for (auto&& combination : GenerateAllCombinations(k, 2)) {
			ones[n].push_back(CalculateProductOfCombinationInRow(combination, ones[n]));
		}

		for (size_t i = 0; i < k; ++i) {
			ones[n].push_back(std::pow(ones[n][i], 2) - S);
		}

		assert(ones[n].size() == M);
	}
}

std::vector<std::vector<double>> CreateUnitPlanningMatrix(size_t k) {
	auto ones = CreateUnitPlanningMatrixCore(k);
	ProceedUnitPlanningMatrixCore(k, ones);
	return ones;
}

template <size_t k>
NonlinearCoefficients<k> CalculateCoefficients(const OccdTable& table) {
	NonlinearCoefficients<k> coefficients{};

	const auto ones = CreateUnitPlanningMatrix(k);
	assert(table.size() == ones.size());

	for (size_t i = 0; i < table.size(); ++i) {
		assert(coefficients.M() == ones[i].size() + 1);

		const auto y = table[i].y_mean;
		coefficients[0] += y;
		for (size_t j = 0; j + 1 < coefficients.M(); ++j) {
			coefficients[j + 1] += ones[i][j] * y;
		}
	}
	for (size_t i = 0; i < coefficients.M(); ++i) {
		coefficients[i] /= static_cast<double>(table.size());
	}

	return coefficients;
}

template <size_t k>
double CalculateNonlinearRegression(const NonlinearCoefficients<k>& coefficients, const std::vector<double>& factors) {
	assert(factors.size() == k);

	double y = 0.0;
	for (size_t i = 0, j = 0; i <= 2; ++i) {
		for (auto&& combination : GenerateAllCombinations(k, i)) {
			y += coefficients[j++] * CalculateProductOfCombinationInRow(combination, factors);
		}
	}
	for (size_t i = 0; i < k; ++i) {
		y += coefficients.a(i, i) * std::pow(factors[i], 2);
	}

	return y;
}

template <size_t k>
std::vector<double> CalculateYUsingRegression(const NonlinearCoefficients<k>& coefficients) {
	const auto ones = CreateUnitPlanningMatrixCore(k);

	std::vector<double> y_hat;
	for (size_t i = 0; i < ones.size(); ++i) {
		assert(ones[i].size() == k);
		y_hat.push_back(CalculateNonlinearRegression(coefficients, ones[i]));
	}
	return y_hat;
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

OccdResult OrthogonalCentralCompositeDesign(const OccdParameters& params) {
	constexpr size_t k = 6;
	const auto N0 = static_cast<size_t>(std::pow(2ul, k));
	const auto N = N0 + 2 * k + 1;
	const auto S = std::sqrt(static_cast<double>(k) / static_cast<double>(N));

	const auto ones = CreateUnitPlanningMatrix(k);
	OccdResult result;
	for (size_t i = 0; i < N; ++i) {
		const DotParameters dot_params = {
			.lambda1       = params.lambda1      .Choose(ones[i][0]),
			.lambda2       = params.lambda2      .Choose(ones[i][1]),
			.mu1           = params.mu1          .Choose(ones[i][2]),
			.mu2           = params.mu2          .Choose(ones[i][3]),
			.sigma_lambda1 = params.sigma_lambda1.Choose(ones[i][4]),
			.sigma_lambda2 = params.sigma_lambda2.Choose(ones[i][5]),
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

			.x12mS = params.lambda1.Choose(1.0 - S),
			.x22mS = params.lambda2.Choose(1.0 - S),
			.x32mS = params.mu1.Choose(1.0 - S),
			.x42mS = params.mu1.Choose(1.0 - S),
			.x52mS = params.sigma_lambda2.Choose(1.0 - S),
			.x62mS = params.sigma_lambda2.Choose(1.0 - S),

			.y_mean = y_mean,
			.y_var  = y_var,

			.y_hat  = 0.0,
			.dy_hat = 0.0,
		});
	}

	result.coefficients = CalculateCoefficients<6>(result.table);
	const auto y_hat = CalculateYUsingRegression(result.coefficients);
	for (size_t i = 0; i < N; ++i) {
		result.table[i].y_hat  = y_hat[i];
		result.table[i].dy_hat = std::abs(result.table[i].y_mean - y_hat[i]);
	}
	return result;
}

std::vector<double> NormalizeFactors(const OccdParameters& occd_params, const DotParameters& dot_params) {
	return {
		occd_params.lambda1.Norm(dot_params.lambda1),
		occd_params.lambda2.Norm(dot_params.lambda2),
		occd_params.mu1.Norm(dot_params.mu1),
		occd_params.mu2.Norm(dot_params.mu2),
		occd_params.sigma_lambda1.Norm(dot_params.sigma_lambda1),
		occd_params.sigma_lambda2.Norm(dot_params.sigma_lambda2),
	};
}

double CalculateDotWithRegression(const NonlinearCoefficients<6>& coefficients, const std::vector<double>& factors) {
	return CalculateNonlinearRegression(coefficients, factors);
}

double CalculateDot(const DotParameters& dot_params, size_t times) {
	return CalculateMean(CalculateY(dot_params, times));
}
