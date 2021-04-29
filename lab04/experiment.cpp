#include "experiment.hpp"

#include <cassert>
#include <functional>
#include <numeric>

#include "random.hpp"
#include "simulate.hpp"
#include "statistics.hpp"


template <typename T>
T CalculateProductOfCombinationInRow(const std::vector<size_t>& combination, const std::vector<T>& row, size_t offset = 0) {
	return std::accumulate(combination.begin(), combination.end(), T{1}, [&row, offset](auto&& acc, auto&& cur) {
		return acc * row[offset + static_cast<size_t>(cur)];
	});
}


template <size_t n>
class OrthogonalCentralCompositeDesign;


template <size_t n>  // количество факторов
class NonlinearCoefficients {
public:
	NonlinearCoefficients()
		: M_(CalculateM())
		, a_(M_)
	{}

	NonlinearCoefficients(const std::vector<double>& coefficients)
		: M_(CalculateM())
		, a_(coefficients)
	{
		assert(a_.size() == M_);
	}

	double& operator[](size_t i) {
		return a_[i];
	}

	double operator[](size_t i) const {
		return a_[i];
	}

	double a(size_t i, size_t j) const {
		return i == j ? aii(i) : aij(i, j);
	}

	double CalculateRegression(const std::vector<double>& factors) const {
		assert(factors.size() >= n);

		double y = 0.0;

		for (size_t i = 0, j = 0; i <= 2; ++i) {
			for (auto&& combination : GenerateAllCombinations(n, i)) {
				y += a_[j++] * CalculateProductOfCombinationInRow(combination, factors);
			}
		}

		const auto Nkernel = static_cast<size_t>(std::pow(2ul, n));
		const auto Nalpha = 2ul * n;
		const auto N0 = 1ul;
		const auto N = Nkernel + Nalpha + N0;
		const auto S = std::sqrt(static_cast<double>(Nkernel) / static_cast<double>(N));

		for (size_t i = 0; i < n; ++i) {
			y += aii(i) * (std::pow(factors[i], 2) - S);
		}

		return y;
	}

private:
	const size_t M_;  // количество коэффициентов полинома нелинейной регрессии
	std::vector<double> a_;  // коэффициенты полинома нелинейной регрессии


	static size_t CalculateM() {
		return CalculateBinomialCoefficient(n, 2) + 2 * n + 1;
	}

	double aij(size_t i, size_t j) const {
		const std::array<size_t, 2> indices_array{i, j};
		const size_t index = 1 + n + CalculateIndexOfCombination(n, indices_array) - 1;
		return a_[index];
	}

	double aii(size_t i) const {
		const size_t index = 1 + n + 2 * n + i;
		return a_[index];
	}

	friend class OrthogonalCentralCompositeDesign<n>;
	friend class ExperimentPlanning;
	friend ExperimentResult RunExperiment(const ExperimentParameters& params);
};


template <size_t n>  // количество факторов
class OrthogonalCentralCompositeDesign {
public:
	OrthogonalCentralCompositeDesign()
		: Nkernel_{static_cast<size_t>(std::pow(2ul, n))}
		, Nalpha_{2 * n}
		, N0_{1}
		, N_{Nkernel_ + Nalpha_ + N0_}

		, S_{std::sqrt(static_cast<double>(Nkernel_) / static_cast<double>(N_))}
		, alpha_{std::sqrt(static_cast<double>(Nkernel_) / 2.0 * (1.0 / S_ - 1.0))}

		, all_2_combinations_of_n_{GenerateAllCombinations(n, 2)}

		, unit_planning_matrix_{}
		, coefficients_{}
	{
		InitializeUnitPlanningMatrixCore();
		ProceedUnitPlanningMatrixCore();
	}

	void CalculateCoefficients(const std::vector<double>& y) {
		assert(y.size() == N_);

		for (size_t i = 0; i < N_; ++i) {
			coefficients_[0] += y[i];
			for (size_t j = 0; j < coefficients_.M_ - 1; ++j) {
				coefficients_[j + 1] += unit_planning_matrix_[i][j] * y[i];
			}
		}

		for (auto&& coefficient : coefficients_.a_) {
			coefficient /= static_cast<double>(N_);
		}
	}

	std::vector<double> CalculateYUsingRegression() const {
		std::vector<double> y_hat;
		std::transform(unit_planning_matrix_.begin(), unit_planning_matrix_.end(), std::back_inserter(y_hat),
			[&](auto&& row) { return coefficients_.CalculateRegression(row); });
		return y_hat;
	}

private:
	const size_t Nkernel_;  // количество опытов в ядре плана (при ПФЭ)
	const size_t Nalpha_;  // количество звёздных точек
	const size_t N0_;  // количество точек в центре плана
	const size_t N_;  // количество опытов в ОЦКП

	const double S_;  // коэффициент, обеспечивающий ортогональность ЦКП
	const double alpha_;  // длина звёздного плеча

	const std::vector<std::vector<size_t>> all_2_combinations_of_n_;  // множество всех сочетаний из n по 2

	std::vector<std::vector<double>> unit_planning_matrix_;  // безразмерная матрица планирования без столбцов "№", "x0" и "y"...
	NonlinearCoefficients<n> coefficients_;  // коэфициенты нелинейной регрессии


	void InitializeUnitPlanningMatrixCore() {
		unit_planning_matrix_.reserve(N_);

		for (size_t i = 0; i < Nkernel_; ++i) {
			unit_planning_matrix_.push_back(std::vector<double>{});
			unit_planning_matrix_.back().reserve(n);
			for (size_t j = 0; j < n; ++j) {
				const size_t bit = 1ul << j;
				unit_planning_matrix_.back().push_back(i & bit ? +1.0 : -1.0);
			}
		}

		for (size_t i = 0; i < n; ++i) {
			unit_planning_matrix_.emplace_back(n);
			unit_planning_matrix_.back()[i] = alpha_;
			unit_planning_matrix_.emplace_back(n);
			unit_planning_matrix_.back()[i] = -alpha_;
		}

		for (size_t i = 0; i < N0_; ++i) {
			unit_planning_matrix_.emplace_back(n);
		}
	}

	void ProceedUnitPlanningMatrixCore() {
		for (auto&& row : unit_planning_matrix_) {
			row.reserve(coefficients_.M_ - 1);

			for (auto&& combination : all_2_combinations_of_n_) {
				row.push_back(CalculateProductOfCombinationInRow(combination, row));  // xi * xj
			}

			for (size_t i = 0; i < n; ++i) {
				row.push_back(std::pow(row[i], 2) - S_);  // xi^2 - S
			}

			assert(row.size() == coefficients_.M_ - 1);
		}
	}

	friend class ExperimentPlanning;
	friend ExperimentResult RunExperiment(const ExperimentParameters& params);
};


class ExperimentPlanning {
public:
	explicit ExperimentPlanning(const ExperimentParameters& params)
		: occd_{}
		, planning_matrix_{}
	{
		InitializePlanningMatrix(params);
	}

	static std::vector<double> CalculateY(const DotParameters& dot_params, size_t times) {
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

private:
	OrthogonalCentralCompositeDesign<6> occd_;

	std::vector<std::vector<double>> planning_matrix_;


	void InitializePlanningMatrix(const ExperimentParameters& params) {
		std::vector<double> ys;
		ys.reserve(occd_.N_);

		for (size_t i = 0; i < occd_.N_; ++i) {
			const DotParameters dot_params = {
				.lambda1       = params.lambda1      .Choose(occd_.unit_planning_matrix_[i][0]),
				.lambda2       = params.lambda2      .Choose(occd_.unit_planning_matrix_[i][1]),
				.mu1           = params.mu1          .Choose(occd_.unit_planning_matrix_[i][2]),
				.mu2           = params.mu2          .Choose(occd_.unit_planning_matrix_[i][3]),
				.sigma_lambda1 = params.sigma_lambda1.Choose(occd_.unit_planning_matrix_[i][4]),
				.sigma_lambda2 = params.sigma_lambda2.Choose(occd_.unit_planning_matrix_[i][5]),
			};

			const auto y = CalculateMean(CalculateY(dot_params, params.times));

			planning_matrix_.push_back(std::vector<double>{});
			planning_matrix_.back().push_back(static_cast<double>(i + 1));  // №
			planning_matrix_.back().push_back(1.0);  // x0
			std::copy(occd_.unit_planning_matrix_[i].begin(), occd_.unit_planning_matrix_[i].end(), std::back_inserter(planning_matrix_.back()));
			planning_matrix_.back().push_back(y);
			ys.push_back(y);
		}

		occd_.CalculateCoefficients(ys);

		const auto y_hat = occd_.CalculateYUsingRegression();
		for (size_t i = 0; i < occd_.N_; ++i) {
			planning_matrix_[i].push_back(y_hat[i]);
			planning_matrix_[i].push_back(std::abs(planning_matrix_[i][1 + occd_.coefficients_.M_] - y_hat[i]));
		}
	}

	static std::vector<double> SimulateNTimes(const SimulateParams& params, size_t times) {
		std::vector<double> y;
		for (size_t i = 0; i < times; ++i) {
			const auto result = Simulate(params);
			y.push_back(result.average_waiting);
		}
		return y;
	}

	friend ExperimentResult RunExperiment(const ExperimentParameters& params);
};

ExperimentResult RunExperiment(const ExperimentParameters& params) {
	ExperimentPlanning experiment(params);
	return {
		.planning_matrix = experiment.planning_matrix_,
		.coefficients = experiment.occd_.coefficients_.a_,
	};
}

std::vector<double> NormalizeFactors(const ExperimentParameters& exp_params, const DotParameters& dot_params) {
	return {
		exp_params.lambda1.Norm(dot_params.lambda1),
		exp_params.lambda2.Norm(dot_params.lambda2),
		exp_params.mu1.Norm(dot_params.mu1),
		exp_params.mu2.Norm(dot_params.mu2),
		exp_params.sigma_lambda1.Norm(dot_params.sigma_lambda1),
		exp_params.sigma_lambda2.Norm(dot_params.sigma_lambda2),
	};
}

double CalculateDotWithRegression(const std::vector<double>& coefficients, const std::vector<double>& factors) {
	return NonlinearCoefficients<6>{coefficients}.CalculateRegression(factors);
}

double CalculateDot(const DotParameters& dot_params, size_t times) {
	return CalculateMean(ExperimentPlanning::CalculateY(dot_params, times));
}
