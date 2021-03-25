#pragma once

#include <array>


size_t CalculateBinomialCoefficient(size_t n, size_t k);

template <size_t k>
size_t CalculateIndexOfCombination(size_t n, const std::array<size_t, k>& combination) {
	// https://math.stackexchange.com/questions/1363239/fast-way-to-get-a-position-of-combination-without-repetitions

	// combination starts with 1 up to n
	size_t index = CalculateBinomialCoefficient(n, k);
	for (size_t i = 0; i < k; ++i) {
		index -= CalculateBinomialCoefficient(n - combination[i], k - i);
	}

	return index;
}
