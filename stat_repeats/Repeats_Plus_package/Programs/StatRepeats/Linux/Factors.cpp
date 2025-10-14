#include <utility>
#include <tuple>
#include "Factors.h"

using namespace std;

Factors::Factors(const Primes& primes) :_primes(primes)
{
	_smallestDivisor.push_back(0);
	_smallestDivisor.push_back(1);
	_smallestDivisor.push_back(2);
}

size_t Factors::GetSmallestDivisor(size_t num) const
{
	for(size_t i = _smallestDivisor.size(); i <= num ;i++) {
		size_t primeNo = 1;
		size_t divisor;
		for(;;) {
			size_t prime = _primes[primeNo];
			if (prime*prime > i) {
				divisor = i;
				break;
			}
			else if ((i%prime) == 0) {
				divisor = prime;
				break;
			}
			else {
				++primeNo;
			}
		}
		_smallestDivisor.push_back(divisor);
	}

	return _smallestDivisor[num];
}

void Factors::FillFactors(const size_t number, vector<size_t>& factors) const {
	factors.clear ();

	pair<size_t, size_t> primeFactorsStorage[100]; // (prime, count), 100 is more than needed as the most primes we can have is 32 for 32-bit unsigneds
	pair<size_t, size_t>* const primeFactorsStart = primeFactorsStorage;
	pair<size_t, size_t>* primeFactorsEnd = primeFactorsStorage;

	size_t num = number;
	while (num > 1) {
		auto smallestDivisor = GetSmallestDivisor(num);
		if (primeFactorsStart == primeFactorsEnd || (primeFactorsEnd - 1)->first != smallestDivisor) {
			*primeFactorsEnd++ = make_pair(smallestDivisor, 1u);
		}
		else {
			(primeFactorsEnd - 1)->second++;
		}
		num /= smallestDivisor;
	}

	auto combinations = size_t(1);
	for (auto x = primeFactorsStart; x != primeFactorsEnd; ++x) {
		combinations *= x->second + 1u;
	}

	for (size_t combination = 0; combination < combinations; ++combination) {
		auto factor = size_t(1);
		auto comb = combination;
		for (auto x = primeFactorsStart; x != primeFactorsEnd; ++x) {
			if (comb == 0) break; 

			auto count = comb % (x->second + 1);
			for (auto y = 0u; y < count; ++y) {
				factor *= x->first;
			}

			if (factor > number / factor) break;

			comb /= x->second + 1;
		}

		if (factor <= number / factor) {
			factors.push_back (factor);
		}
	}
}

