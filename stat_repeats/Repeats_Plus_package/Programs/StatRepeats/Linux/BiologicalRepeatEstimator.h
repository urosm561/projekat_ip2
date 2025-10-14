#ifndef BIOLOGICAL_REPEAT_ESTIMATOR_H
#define BIOLOGICAL_REPEAT_ESTIMATOR_H

#include <vector>
#include <cmath>
#include <math.h>
#include "Factors.h"

class BiologicalRepeatEstimator
{   
	size_t _seqLen;
	size_t _fragLen;
	size_t _fragCount;
	double _lnFragCombo;
	double _minimalFragmentProbability;
	std::vector<double> _lnfaktoriel;
	double _lambda;
	double _lnLambda;
	double _lnMinimalFragmentProbability;
	double _lnOneMinusMfP;
	std::vector<double> _lnDistribution;
	size_t alphabetSize_;

	double GetLnFactoriel(size_t k);
	double LnBinomialDistribution(size_t k, size_t n);
	double GetLnBinomRepeats(size_t k, size_t n, size_t comb, double p);
	double LnPoissonDistribution(size_t k);
	double CalculateVarianceShift(const size_t k);
	double NormalDistribution(const size_t k, const double mean, const double standardDeviation);
	double LnDistribution(const size_t k);

	double ForNonRepeats(size_t searchStart, size_t last, int step);

	const Factors* _factors;

public:
	BiologicalRepeatEstimator(size_t seqLen, size_t fragLen, size_t alphabetSize, const Factors& factors);
	double Compute(size_t searchStart);
};

#endif