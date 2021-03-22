#include "experiment.hpp"

#include <cmath>

#include "random.hpp"
#include "simulate.hpp"
#include "statistics.hpp"


static constexpr size_t N = 8;
static constexpr int ones[N][3] = {
	{-1, -1, -1},
	{-1, -1,  1},
	{-1,  1, -1},
	{-1,  1,  1},
	{ 1, -1, -1},
	{ 1, -1,  1},
	{ 1,  1, -1},
	{ 1,  1,  1},
};

FfeCoefficients CalculateCoefficients(const QVector<FfeTableRow>& rows) {
	Q_ASSERT(rows.size() == N);

	FfeCoefficients coefficients;
	coefficients.a.a0 = std::accumulate(rows.begin(), rows.end(), 0.0, [](double acc, const FfeTableRow& row) {
		return acc += row.y_mean;
	});
	coefficients.a.a0 /= N;
	for (int i = 0; i < 3; ++i) {
		double ai = 0.0;
		for (size_t n = 0; n < N; ++n) {
			ai += ones[n][i] * rows[static_cast<int>(n)].y_mean;
		}
		coefficients.as[i + 1] = ai / N;
	}
	for (int i = 0; i < 2; ++i) {
		for (int j = i + 1; j < 3; ++j) {
			double aij = 0.0;
			for (size_t n = 0; n < N; ++n) {
				aij += ones[n][i] * ones[n][j] * rows[static_cast<int>(n)].y_mean;
			}
			coefficients.as[3 + i + j] = aij / N;
		}
	}
	return coefficients;
}

double CalculateWithCoefficients(const FfeCoefficients& c, double x1, double x2, double x3) {
	return c.a.a0 + c.a.a1 * x1 + c.a.a2 * x2 + c.a.a3 * x3 + c.a.a12 * x1 * x2 + c.a.a13 * x1 * x3 + c.a.a23 * x2 * x3;
}

QVector<double> CalculateYHat(const FfeCoefficients& c) {
	QVector<double> y_hat;
	for (size_t i = 0; i < N; ++i) {
		y_hat.append(CalculateWithCoefficients(c, ones[i][0], ones[i][1], ones[i][2]));
	}
	return y_hat;
}

QVector<double> CalculateYHatLinear(const FfeCoefficients& c) {
	const auto calc = [&c](double x1, double x2, double x3) {
		return c.a.a0 + c.a.a1 * x1 + c.a.a2 * x2 + c.a.a3 * x3;
	};

	QVector<double> y_hat_linear;
	for (size_t i = 0; i < N; ++i) {
		y_hat_linear.append(calc(ones[i][0], ones[i][1], ones[i][2]));
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

#include <QtDebug>

FfeResult FullFactorialExperiment(const FfeParameters& params) {
	FfeResult result;
	double sum_var = 0.0;
	double max_var = 0.0;
	for (int i = 0; i < 8; ++i) {
		const double lambda = (i & 0b100) == 0b100 ? params.lambda_max : params.lambda_min;
		const double sigma_lambda = (i & 0b10) == 0b10 ? params.sigma_lambda_max : params.sigma_lambda_min;
		const double mu = (i & 0b1) == 0b1 ? params.mu_max : params.mu_min;

		const auto [a, b] = uniform_parameters_from_mean_and_std(1.0 / lambda, sigma_lambda);
		const auto [k, l] = weibull_parameters_from_mean(1.0 / mu);
		const auto y = SimulateNTimes({
			.a = a,
			.b = b,

			.k = k,
			.lambda = l,

			.t = 1000.0
		}, params.times);

		const auto [m, v] = weibull_mean_and_variance(k, l);
		qDebug() << "mu" << mu << ", k" << k << ", l" << l << ", mean" << m << ", variance " << v;

		const auto [y_mean, y_var] = CalculateMeanAndVariance(y);

		sum_var += y_var;
		max_var = std::max(y_var, max_var);

		result.rows.push_back({
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
	result.reproducibility_var = result.cochran_test / 8.0;
	if (result.cochran_test < 0.3910) {
		qDebug() << "result.cochran_test < 0.3910";
	}

	result.coefficients = CalculateCoefficients(result.rows);
	const auto y_hat = CalculateYHat(result.coefficients);
	const auto y_hat_linear = CalculateYHatLinear(result.coefficients);
	double diff_sum_squared = 0.0;
	for (int i = 0; i < 8; ++i) {
		result.rows[i].partial_nonlinear = y_hat[i];
		result.rows[i].dpn = result.rows[i].y_mean - y_hat[i];
		result.rows[i].linear = y_hat_linear[i];
		result.rows[i].dl = result.rows[i].y_mean - y_hat_linear[i];
		diff_sum_squared += std::pow(result.rows[i].dpn, 2);
	}
	result.adequacy_var = static_cast<double>(params.times) / 8.0 * diff_sum_squared;
	result.f_test = result.adequacy_var / result.reproducibility_var;

	return result;
}

DotResult CalculateDot(const FfeParameters& params, const FfeCoefficients& c, double lambda, double sigma_lambda, double mu) {
	const double x1 = 2 * (lambda - params.lambda_min) / (params.lambda_max - params.lambda_min) - 1.0;
	const double x2 = 2 * (sigma_lambda - params.sigma_lambda_min) / (params.sigma_lambda_max - params.sigma_lambda_min) - 1.0;
	const double x3 = 2 * (mu - params.mu_min) / (params.mu_max - params.mu_min) - 1.0;
	DotResult result;
	result.estimated_y = CalculateWithCoefficients(c, x1, x2, x3);

	const auto [a, b] = uniform_parameters_from_mean_and_std(1.0 / lambda, sigma_lambda);
	const auto [k, l] = weibull_parameters_from_mean(1.0 / mu);
	const auto y = SimulateNTimes({
		.a = a,
		.b = b,

		.k = k,
		.lambda = l,

		.t = 1000.0
	}, params.times);
	result.actual_y = CalculateMean(y);

	return result;
}
