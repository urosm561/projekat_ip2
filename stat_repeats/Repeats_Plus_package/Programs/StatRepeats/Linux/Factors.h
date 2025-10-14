#ifndef FACTORS_H
#define FACTORS_H

#include <vector>
#include "Primes.h"

class Factors{
	mutable std::vector<size_t> _smallestDivisor;
	const Primes& _primes;

	size_t GetSmallestDivisor(size_t num) const;
public:
	Factors(const Primes& primes);
	void FillFactors(const size_t num, std::vector<size_t>& factors) const;

	void operator = (const Factors&) = delete;
	Factors(const Factors&) = delete;
};

#endif