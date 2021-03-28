#pragma once

#include <vector>

struct SimulateParams {
	double a1, b1;
	double a2, b2;
	double k, lambda;

	double t;
};

struct SimulateResult {
	double average_waiting;
};

SimulateResult Simulate(const SimulateParams& params);
