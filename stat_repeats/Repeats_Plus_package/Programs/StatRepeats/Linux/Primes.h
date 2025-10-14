#ifndef PRIMES_H
#define PRIMES_H

#include <cstdlib>
#include <vector>

class Primes
{
	mutable std::vector<unsigned> _primes;

public:
	Primes();
	unsigned operator [] (std::size_t index) const;
};

#endif