#include "experiment.hpp"

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

double CalculateWithCoefficientsNonlinear(const PartialNonlinearCoefficients<3>& c, double x1, double x2, double x3) {
	return c[0] + c[1] * x1 + c[2] * x2 + c[3] * x3 + c[4] * x1 * x2 + c[5] * x1 * x3 + c[6] * x2 * x3 + c[7] * x1 * x2 * x3;
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

std::vector<double> CalculateYHat(const PartialNonlinearCoefficients<3>& coefficients) {
	const auto ones = CalculatePlanningMatrixCore(3);

	std::vector<double> y_hat;
	for (size_t i = 0; i < coefficients.N(); ++i) {
		y_hat.push_back(CalculateWithCoefficientsNonlinear(coefficients, ones[i][0], ones[i][1], ones[i][2]));
	}

	return y_hat;
}

template <size_t k>
std::vector<double> CalculateYHatLinear(const PartialNonlinearCoefficients<k>& coefficients) {
	const auto ones = CalculatePlanningMatrixCore(k);

	std::vector<double> y_hat_linear;
	for (size_t i = 0; i < coefficients.N(); ++i) {
		assert(ones[i].size() == k);

		std::vector<double> factors;
		factors.reserve(k);
		for (int one : ones[i]) {
			factors.push_back(static_cast<double>(one));
		}

		y_hat_linear.push_back(CalculateWithCoefficientsLinear(coefficients, factors));
	}
	return y_hat_linear;
}

std::vector<double> SimulateNTimes(const SimulateParams& params, size_t times) {
	std::vector<double> y;
	for (size_t i = 0; i < times; ++i) {
		const auto result = Simulate(params);
		y.push_back(result.average_waiting);
	}
	return y;
}

std::vector<double> CalculateY(double lambda, double sigma_lambda, double mu, size_t times) {
	const auto [a, b] = uniform_parameters_from_mean_and_std(1.0 / lambda, sigma_lambda);
	const auto [k, l] = weibull_parameters_from_mean(1.0 / mu);
	const auto y = SimulateNTimes({
		.a = a,
		.b = b,

		.k = k,
		.lambda = l,

		.t = 1000.0
	}, times);
	return y;
}

#include <QtDebug>

FfeResult FullFactorialExperiment(const FfeParameters& params) {
	FfeResult result;
	double sum_var = 0.0;
	double max_var = 0.0;
	for (int i = 0; i < 8; ++i) {
		const double lambda = (i & 0b100) == 0b100 ? params.lambda_max : params.lambda_min;
		const double sigma_lambda = (i & 0b10) == 0b10 ? params.sigma_lambda_max : params.sigma_lambda_min;
		const double mu = (i & 0b1) == 0b1 ? params.mu_max : params.mu_min;

		const auto y = CalculateY(lambda, sigma_lambda, mu, params.times);
		const auto [y_mean, y_var] = CalculateMeanAndVariance(y);

		sum_var += y_var;
		max_var = std::max(y_var, max_var);

		result.table.push_back({
			.index = i + 1,
			.x1 = lambda,
			.x2 = sigma_lambda,
			.x3 = mu,
			.y_mean = y_mean,
			.y_var = y_var,

			.partial_nonlinear = 0.0,
			.dpn = 0.0,
			.linear = 0.0,
			.dl = 0.0,
		});
	}
	result.cochran_test = max_var / sum_var;
	result.reproducibility_var = sum_var / 8.0;
	if (result.cochran_test > 0.3910) {
		qDebug() << "result.cochran_test > 0.3910";
	}

	result.coefficients = CalculateCoefficients<3>(result.table);
	const auto y_hat = CalculateYHat(result.coefficients);
	const auto y_hat_linear = CalculateYHatLinear(result.coefficients);
	double diff_sum_squared = 0.0;
	for (size_t i = 0; i < result.coefficients.N(); ++i) {
		result.table[i].partial_nonlinear = y_hat[i];
		result.table[i].dpn = result.table[i].y_mean - y_hat[i];
		result.table[i].linear = y_hat_linear[i];
		result.table[i].dl = result.table[i].y_mean - y_hat_linear[i];
		diff_sum_squared += std::pow(result.table[i].dpn, 2);
	}
	result.adequacy_var = static_cast<double>(params.times) / 8.0 * diff_sum_squared;
	result.f_test = result.adequacy_var / result.reproducibility_var;

	return result;
}

DotResult CalculateDot(const FfeParameters& params, const PartialNonlinearCoefficients<3>& coefficients, double lambda, double sigma_lambda, double mu) {
	const auto norm = [](double x, double x_min, double x_max) {
		return 2.0 * (x - x_min) / (x_max - x_min) - 1.0;
	};

	const double x1 = norm(lambda, params.lambda_min, params.lambda_max);
	const double x2 = norm(sigma_lambda, params.sigma_lambda_min, params.sigma_lambda_max);
	const double x3 = norm(mu, params.mu_min, params.mu_max);
	DotResult result;
	result.estimated_y = CalculateWithCoefficientsNonlinear(coefficients, x1, x2, x3);

	const auto y = CalculateY(lambda, sigma_lambda, mu, params.times);
	result.actual_y = CalculateMean(y);

	return result;
}
