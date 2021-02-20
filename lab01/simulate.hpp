#pragma once

struct SimulateParams {
	double a, b;
	double k, lambda;

	double t;
};

struct SimulateResult {
	// i don't know
};

SimulateResult Simulate(const SimulateParams& params);
