#pragma once

#include <algorithm>
#include <array>


constexpr size_t CalculateBinomialCoefficient(size_t n, size_t k) {
	if (n < k) {
		return 0;
	}

	k = std::min(k, n - k);

	size_t result = 1;
	for (size_t i = 0; i < k; ++i) {
		result *= n - i;
		result /= i + 1;
	}

	return result;
}

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

std::vector<std::vector<size_t>> GenerateAllCombinations(size_t n, size_t k);
