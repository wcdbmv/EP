#include "combinatorics.hpp"

#include <cassert>


std::vector<std::vector<size_t>> GenerateAllCombinations(size_t n, size_t k) {
	// http://rosettacode.org/wiki/Combinations#C.2B.2B

	const auto size = CalculateBinomialCoefficient(n, k);
	std::vector<std::vector<size_t>> result;
	result.reserve(size);

	std::string bitmask(k, 1);  // K leading 1's
	bitmask.resize(n, 0);       // N-K trailing 0's

	// save integers and permute bitmask
	do {
		result.push_back({});
		result.back().reserve(k);

		for (size_t i = 0; i < n; ++i) {  // [0..N-1] integers
			if (bitmask[i]) {
				result.back().push_back(i);
			}
		}

		assert(result.back().size() == k);
	} while (std::prev_permutation(bitmask.begin(), bitmask.end()));

	assert(result.size() == size);
	return result;
}
