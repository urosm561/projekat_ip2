#include "Primes.h"

using namespace std;

Primes::Primes ()
{
	_primes.push_back(1);
}

unsigned Primes::operator [] (size_t index) const
{
	for (size_t x = _primes.size(); x <= index; x++){
		unsigned candidate = _primes[x - 1] + 1;
		for (;;) {
			bool isPrime = true;
			for (size_t i = 1; i < x; i++){
				unsigned div = _primes[i];
				if (div*div > candidate) break;
				if ((candidate % div) == 0) {
					isPrime = false;
					break;
				}
			}

			if (isPrime) break;
			candidate++;
		}

		_primes.push_back(candidate);
	}

	return _primes[index];
}
