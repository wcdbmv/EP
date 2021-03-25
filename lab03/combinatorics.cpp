#include "combinatorics.hpp"

#include <algorithm>


size_t CalculateBinomialCoefficient(size_t n, size_t k) {
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
