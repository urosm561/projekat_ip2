#ifndef BIOLOGICAL_PALINDROME_ESTIMATOR_H
#define BIOLOGICAL_PALINDROME_ESTIMATOR

#include <vector>
#include <cmath>
#include <math.h>
#include "Factors.h"
#include "MathematicalPalindromeEstimator.h"
#include "BiologicalRepeatEstimator.h"

class BiologicalPalindromeEstimator {
	MathematicalPalindromeEstimator matPalEst_;
	BiologicalRepeatEstimator bioRepEst_;
	size_t fragLen_;

public:
	BiologicalPalindromeEstimator(size_t seqLen, size_t fragLen, size_t alphabetSize, const Factors& factors) :
		matPalEst_(seqLen, fragLen, alphabetSize, factors), bioRepEst_(seqLen, fragLen, alphabetSize, factors)
	{
		fragLen_ = fragLen;
	}

	double Compute(size_t pairs) {
		return
			(fragLen_ & 1) == 0 ?
			matPalEst_.Compute(pairs) :
			bioRepEst_.Compute(pairs)
			;
	}

	size_t getFragLen() const {
		return fragLen_;
	}
};

#endif