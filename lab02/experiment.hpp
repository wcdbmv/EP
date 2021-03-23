#pragma once

#include <tuple>
#include <QVector>

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

struct FfeResult {
	QVector<FfeTableRow> rows;
	double cochran_test;
	double reproducibility_var;
	double adequacy_var;
	double f_test;
	PlanningMatrix3 planning_matrix;
};

struct DotResult {
	double estimated_y;
	double actual_y;
};

FfeResult FullFactorialExperiment(const FfeParameters& params);
DotResult CalculateDot(const FfeParameters& params, const PlanningMatrix3& m, double lambda, double sigma_lambda, double mu);
