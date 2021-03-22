#pragma once

#include <vector>

struct SimulateParams {
	double a, b;
	double k, lambda;

	double t;
};

struct SimulateResult {
	double average_waiting;
};

SimulateResult Simulate(const SimulateParams& params);
