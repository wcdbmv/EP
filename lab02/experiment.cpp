#include "experiment.hpp"

#include <cmath>

#include "random.hpp"
#include "simulate.hpp"
#include "statistics.hpp"


static constexpr size_t N = 8;
static constexpr int ones[N][7] = {
	// {x1, x2, x3, x1x2, x1x3, x2x3, x1x2x3}
	{-1, -1, -1, +1, +1, +1, -1},
	{-1, -1,  1, +1, -1, -1, +1},
	{-1,  1, -1, -1, +1, -1, +1},
	{-1,  1,  1, -1, -1, +1, -1},
	{ 1, -1, -1, -1, -1, +1, +1},
	{ 1, -1,  1, -1, +1, -1, -1},
	{ 1,  1, -1, +1, -1, -1, -1},
	{ 1,  1,  1, +1, +1, +1, +1},
};

PlanningMatrix3 CalculateCoefficients(const QVector<FfeTableRow>& rows) {
	Q_ASSERT(rows.size() == N);

	PlanningMatrix3 planning_matrix{};
	for (int i = 0; i < static_cast<int>(N); ++i) {
		const auto y = rows[i].y_mean;
		planning_matrix.a0 += y;
		for (int j = 0; j < 7; ++j) {
			planning_matrix[j + 1] += ones[i][j] * y;
		}
	}
	for (int i = 0; i < 8; ++i) {
		planning_matrix[i] /= N;
	}

	return planning_matrix;
}

double CalculateWithCoefficients(const PlanningMatrix3& m, double x1, double x2, double x3) {
	return m.a0 + m.a1 * x1 + m.a2 * x2 + m.a3 * x3 + m.a12 * x1 * x2 + m.a13 * x1 * x3 + m.a23 * x2 * x3 + m.a123 * x1 * x2 * x3;
}

QVector<double> CalculateYHat(const PlanningMatrix3& m) {
	QVector<double> y_hat;
	for (size_t i = 0; i < N; ++i) {
		y_hat.append(CalculateWithCoefficients(m, ones[i][0], ones[i][1], ones[i][2]));
	}
	return y_hat;
}

QVector<double> CalculateYHatLinear(const PlanningMatrix3& m) {
	const auto calc = [&m](double x1, double x2, double x3) {
		return m.a0 + m.a1 * x1 + m.a2 * x2 + m.a3 * x3;
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
	result.reproducibility_var = sum_var / 8.0;
	if (result.cochran_test > 0.3910) {
		qDebug() << "result.cochran_test > 0.3910";
	}

	result.planning_matrix = CalculateCoefficients(result.rows);
	const auto y_hat = CalculateYHat(result.planning_matrix);
	const auto y_hat_linear = CalculateYHatLinear(result.planning_matrix);
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

DotResult CalculateDot(const FfeParameters& params, const PlanningMatrix3& m, double lambda, double sigma_lambda, double mu) {
	const auto norm = [](double x, double x_min, double x_max) {
		return 2 * (x - x_min) / (x_max - x_min) - 1.0;
	};

	const double x1 = norm(lambda, params.lambda_min, params.lambda_max);
	const double x2 = norm(sigma_lambda, params.sigma_lambda_min, params.sigma_lambda_max);
	const double x3 = norm(mu, params.mu_min, params.mu_max);
	DotResult result;
	result.estimated_y = CalculateWithCoefficients(m, x1, x2, x3);

	const auto y = CalculateY(lambda, sigma_lambda, mu, params.times);
	result.actual_y = CalculateMean(y);

	return result;
}
