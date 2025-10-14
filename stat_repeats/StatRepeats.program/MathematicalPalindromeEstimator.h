#ifndef MATHEMATICAL_PALINDROME_ESTIMATOR_H
#define MATHEMATICAL_PALINDROME_ESTIMATOR_H

#include <vector>
#include <cmath>
#include <math.h>
#include "Factors.h"


class MathematicalPalindromeEstimator {
	size_t _seqLen;
	size_t _fragLen;
	size_t _fragCount;

	double _lnFragCombo;
	double _lnFragComboNonPal;
	double _lnFragComboPal;

	double _minimalFragmentProbability;
	std::vector<double> _lnfaktoriel;
	double _lambda;
	double _lnLambda;

	double _lnMinimalFragmentProbability;
	double _lnOneMinusMfP;

	double _lnSelfMinimalFragmentProbability;
	double _lnSelfOneMinusMfP;

	std::vector<double> _lnDistribution;
	size_t alphabetSize_;

	void EnsureFactoriel(size_t k);
	double LnBinomialDistribution(size_t k, size_t n);
	double LnBinomialSelfDistribution(size_t k, size_t n);
	double LnMatRepBinomialDistribution(size_t k, size_t words);
	double GetLnBinomRepeats(size_t k, size_t n, size_t comb, double p);
	double LnPoissonDistribution(size_t k);
	double CalculateVarianceShift(const size_t k);
	double NormalDistribution(const size_t k, const double mean, const double standardDeviation);
	double LnDistribution(const size_t k);

	double ForNonRepeats(size_t searchStart, size_t last, int step);
	double ForRepeats(size_t searchStart);

	const Factors* _factors;

	double ComputeForNonPalindromes(size_t searchStart);
	double ComputeForPalindromes(size_t num);

public:
	MathematicalPalindromeEstimator(size_t seqLen, size_t fragLen, size_t alphabetSize, const Factors& factors);

	double Compute(size_t pairs) {
		return ComputeForNonPalindromes(pairs) + ComputeForPalindromes(pairs);
	}

	size_t getFragLen() const {
		return _fragLen;
	}
};

#endif