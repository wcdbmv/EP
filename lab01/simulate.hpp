#pragma once

struct SimulateParams {
	double a, b;
	double k, lambda;

	double t;
};

struct SimulateResult {
	double load;
	double average_waiting;
};

SimulateResult Simulate(const SimulateParams& params);
