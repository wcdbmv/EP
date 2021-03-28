#include "experiment.hpp"

#include <functional>
#include <numeric>

#include "random.hpp"
#include "simulate.hpp"
#include "statistics.hpp"


std::vector<std::vector<int>> CalculatePlanningMatrixCore(size_t k) {
	const auto N = static_cast<size_t>(std::pow(2ul, k));

	std::vector<std::vector<int>> ones;
	ones.reserve(N);
	for (size_t i = 0; i < N; ++i) {
		ones.push_back(std::vector<int>{});
		ones[i].reserve(k);
		for (size_t j = 0; j < k; ++j) {
			const size_t bit = 1u << (k - j - 1);
			ones[i].push_back(i & bit ? +1 : -1);
		}
	}

	return ones;
}


std::vector<std::vector<int>> CalculatePlanningMatrix(size_t k) {
	constexpr size_t K = 5;
	assert(1 <= k && k <= 5);

	const auto N = static_cast<size_t>(std::pow(2ul, k));
	const auto M = N - 1;

	auto ones = CalculatePlanningMatrixCore(k);

	// hell no
	for (size_t n = 0; n < N; ++n) {
		ones[n].reserve(M);
		for (size_t i = 2; i <= K; ++i) {
			for (size_t a = 0; a + i < k + 1; ++a) {
				for (size_t b = a + 1; b + i < k + 2; ++b) {
					if (i == 2) {
						ones[n].push_back(ones[n][a] * ones[n][b]);
					} else {
						for (size_t c = b + 1; c + i < k + 3; ++c) {
							if (i == 3) {
								ones[n].push_back(ones[n][a] * ones[n][b] * ones[n][c]);
							} else {
								for (size_t d = c + 1; d + i < k + 4; ++d) {
									ones[n].push_back(ones[n][a] * ones[n][b] * ones[n][c] * ones[n][d]);
								}
							}
						}
					}
				}
			}
		}
		if (k == 5) {
			ones[n].push_back(ones[n][0] * ones[n][1] * ones[n][2] * ones[n][3] * ones[n][4]);
		}
	}

	return ones;
}

template <size_t k>
PartialNonlinearCoefficients<k> CalculateCoefficients(const FfeTable& table) {
	PartialNonlinearCoefficients<k> coefficients{};
	assert(table.size() == coefficients.N());

	const auto ones = CalculatePlanningMatrix(k);
	for (size_t i = 0; i < coefficients.N(); ++i) {
		const auto y = table[i].y_mean;
		coefficients[0] += y;
		for (size_t j = 0; j + 1 < coefficients.N(); ++j) {
			coefficients[j + 1] += ones[i][j] * y;
		}
	}
	for (size_t i = 0; i < coefficients.N(); ++i) {
		coefficients[i] /= static_cast<double>(coefficients.N());
	}

	return coefficients;
}

double CalculateWithCoefficientsNonlinear(const PartialNonlinearCoefficients<5>& coefficients, const std::vector<double>& factors) {
	constexpr size_t K = 5;
	assert(factors.size() == K);

	constexpr size_t I_SINGLES_START_FOR_5 = 1;
	constexpr size_t I_DOUBLES_START_FOR_5 = 6;
	constexpr size_t I_TRIPLES_START_FOR_5 = 16;
	constexpr size_t I_QUADRUPLES_START_FOR_5 = 26;
	constexpr size_t I_QUINTUPLES_START_FOR_5 = 31;

	double y = coefficients[0] + coefficients[I_QUINTUPLES_START_FOR_5] * std::accumulate(factors.begin(), factors.end(), 1.0, std::multiplies<>{});
	size_t i_singles = I_SINGLES_START_FOR_5;
	size_t i_doubles = I_DOUBLES_START_FOR_5;
	size_t i_triples = I_TRIPLES_START_FOR_5;
	size_t i_quadruples = I_QUADRUPLES_START_FOR_5;
	for (size_t a = 0; a < K; ++a) {
		y += coefficients[i_singles++] * factors[a];
		for (size_t b = a + 1; b < K; ++b) {
			y += coefficients[i_doubles++] * factors[a] * factors[b];
			for (size_t c = b + 1; c < K; ++c) {
				y += coefficients[i_triples++] * factors[a] * factors[b] * factors[c];
				for (size_t d = c + 1; d < K; ++d) {
					y += coefficients[i_quadruples++] * factors[a] * factors[b] * factors[c] * factors[d];
				}
			}
		}
	}
	assert(i_singles == I_DOUBLES_START_FOR_5);
	assert(i_doubles == I_TRIPLES_START_FOR_5);
	assert(i_triples == I_QUADRUPLES_START_FOR_5);
	assert(i_quadruples == I_QUINTUPLES_START_FOR_5);

	return y;
}

template <size_t k>
double CalculateWithCoefficientsLinear(const PartialNonlinearCoefficients<k>& coefficients, const std::vector<double>& factors) {
	assert(factors.size() == k);

	double y = coefficients[0];
	for (size_t i = 0; i < k; ++i) {
		y += coefficients[i + 1] * factors[i];
	}

	return y;
}

template <size_t k>
using TCalculateWith = std::function<double(const PartialNonlinearCoefficients<k>& coefficients, const std::vector<double>& factors)>;

template <size_t k>
std::vector<double> CalculateYUsingRegression(const PartialNonlinearCoefficients<k>& coefficients, TCalculateWith<k>&& calculate_with) {
	const auto ones = CalculatePlanningMatrixCore(k);

	std::vector<double> y_hat;
	for (size_t i = 0; i < coefficients.N(); ++i) {
		assert(ones[i].size() == k);

		std::vector<double> factors;
		factors.reserve(k);
		for (int one : ones[i]) {
			factors.push_back(static_cast<double>(one));
		}

		y_hat.push_back(calculate_with(coefficients, factors));
	}
	return y_hat;
}

std::vector<double> CalculateYPartialNonlinear(const PartialNonlinearCoefficients<5>& coefficients) {
	return CalculateYUsingRegression<5>(coefficients, CalculateWithCoefficientsNonlinear);
}

template <size_t k>
std::vector<double> CalculateYLinear(const PartialNonlinearCoefficients<k>& coefficients) {
	return CalculateYUsingRegression<k>(coefficients, CalculateWithCoefficientsLinear<k>);
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
	const auto [k, l] = weibull_parameters_from_mean(1.0 / dot_params.mu);
	const auto y = SimulateNTimes({
		.a1 = a1,
		.b1 = b1,

		.a2 = a2,
		.b2 = b2,

		.k = k,
		.lambda = l,

		.t = 1000.0
	}, times);
	return y;
}

#include <QtDebug>

FfeResult FullFactorialExperiment(const FfeParameters& params) {
	const auto choose = [](int i, int bit_shift, const Range& range) {
		return (i & (1 << bit_shift)) ? range.max : range.min;
	};

	FfeResult result;
	double sum_var = 0.0;
	double max_var = 0.0;
	for (int i = 0; i < 32; ++i) {
		const DotParameters dot_params = {
			.lambda1       = choose(i, 0, params.lambda1),
			.lambda2       = choose(i, 1, params.lambda2),
			.mu            = choose(i, 2, params.mu),
			.sigma_lambda1 = choose(i, 3, params.sigma_lambda1),
			.sigma_lambda2 = choose(i, 4, params.sigma_lambda2),
		};

		const auto y = CalculateY(dot_params, params.times);
		const auto [y_mean, y_var] = CalculateMeanAndVariance(y);

		sum_var += y_var;
		max_var = std::max(y_var, max_var);

		result.table.push_back({
			.index = i + 1,
			.x1 = dot_params.lambda1,
			.x2 = dot_params.lambda2,
			.x3 = dot_params.mu,
			.x4 = dot_params.sigma_lambda1,
			.x5 = dot_params.sigma_lambda2,
			.y_mean = y_mean,
			.y_var = y_var,

			.y_hat = 0.0,
			.dy_hat = 0.0,
			.u_hat = 0.0,
			.du_hat = 0.0,
		});
	}

	result.coefficients = CalculateCoefficients<5>(result.table);
	const auto y_hat = CalculateYLinear(result.coefficients);
	const auto u_hat = CalculateYPartialNonlinear(result.coefficients);
	double y_hat_diff_sum_squared = 0.0;
	double u_hat_diff_sum_squared = 0.0;
	for (size_t i = 0; i < result.coefficients.N(); ++i) {
		result.table[i].y_hat = y_hat[i];
		result.table[i].dy_hat = result.table[i].y_mean - y_hat[i];
		result.table[i].u_hat = u_hat[i];
		result.table[i].du_hat = result.table[i].y_mean - u_hat[i];
		y_hat_diff_sum_squared += std::pow(result.table[i].dy_hat, 2);
		u_hat_diff_sum_squared += std::pow(result.table[i].du_hat, 2);
	}
	result.y_hat_adequacy_var = static_cast<double>(params.times) / 8.0 * y_hat_diff_sum_squared;
	result.u_hat_adequacy_var = static_cast<double>(params.times) / 8.0 * u_hat_diff_sum_squared;
	result.y_hat_f_test = result.y_hat_adequacy_var / result.reproducibility_var;
	result.u_hat_f_test = result.u_hat_adequacy_var / result.reproducibility_var;

	return result;
}

DotResult CalculateDot(const FfeParameters& ffe_params, const PartialNonlinearCoefficients<5>& coefficients, const DotParameters& dot_params) {
	const auto norm = [](double x, const Range& range) {
		return 2.0 * (x - range.min) / (range.max - range.min) - 1.0;
	};

	const std::vector<double> factors{
		norm(dot_params.lambda1, ffe_params.lambda1),
		norm(dot_params.lambda2, ffe_params.lambda2),
		norm(dot_params.mu, ffe_params.mu),
		norm(dot_params.sigma_lambda1, ffe_params.sigma_lambda1),
		norm(dot_params.sigma_lambda2, ffe_params.sigma_lambda2),
	};

	DotResult result{
		.estimated_y = CalculateWithCoefficientsNonlinear(coefficients, factors),
		.actual_y = CalculateMean(CalculateY(dot_params, ffe_params.times)),
	};

	return result;
}
