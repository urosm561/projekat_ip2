#define _USE_MATH_DEFINES

#include <vector>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include "BiologicalRepeatEstimator.h"
#include "RepeatCoefficients.h"

using namespace std;

BiologicalRepeatEstimator::BiologicalRepeatEstimator(size_t seqLen, size_t fragLen, size_t alphabetSize, const Factors& factors)
{
	alphabetSize_ = alphabetSize;
	_factors = &factors;
	_seqLen = seqLen;
	_fragLen = fragLen;
	_fragCount = _seqLen - _fragLen + 1;
	_lnFragCombo = log(alphabetSize_) * _fragLen;
	_minimalFragmentProbability = double(alphabetSize_ - 1) * (alphabetSize_ - 1) / alphabetSize_ / alphabetSize_;
	_lnfaktoriel.push_back(0);
	_lnLambda = log(_fragCount) - _lnFragCombo;
	_lambda = exp(_lnLambda);
	_lnMinimalFragmentProbability = log(_minimalFragmentProbability);
	_lnOneMinusMfP = log(1 - _minimalFragmentProbability);
	_lnDistribution.reserve(_fragCount);

	for (size_t x = 0; x < _fragCount; ++x)
	{
		double pois = LnDistribution(x);
		if (x > _lambda && exp(pois + _lnFragCombo) == 0) break;
		_lnDistribution.push_back(pois);
	}
}

double BiologicalRepeatEstimator::GetLnFactoriel(size_t k)
{
	if (k > 10000) {
		return k * log(k) - k + 0.5*(log(2) + log(k) + log(M_PI));
	}
	else {
		if (k >= _lnfaktoriel.size())
		{
			for (size_t i = _lnfaktoriel.size(); i <= k; ++i)
			{
				_lnfaktoriel.push_back(_lnfaktoriel[i - 1] + log(i));
			}
		}

		return _lnfaktoriel[k];
	}
}

double BiologicalRepeatEstimator::LnBinomialDistribution(size_t k, size_t n)
{
    double result =
        + GetLnFactoriel(n)
        - GetLnFactoriel(n - k)
        - GetLnFactoriel(k)
        + k * _lnMinimalFragmentProbability
        + (n - k) * _lnOneMinusMfP
        ;

    return result;
}

double BiologicalRepeatEstimator::GetLnBinomRepeats(size_t k, size_t n, size_t comb, double p)
{
	if (alphabetSize_ < GenericRepeatCoefficients.size() && n < GenericMatPalCoefficients[alphabetSize_].size()) {
		auto arr = GenericRepeatCoefficients[alphabetSize_][n];

		if (arr != nullptr) {
			return log(arr[k]);
		}
	}

	return LnBinomialDistribution(k, comb);
}

double BiologicalRepeatEstimator::CalculateVarianceShift(const size_t k){
	return (_seqLen - _fragLen + 1) / (exp(_lnFragCombo)) - ((2*_fragLen-1)*_seqLen - 3*_fragLen*_fragLen + 4*_fragLen - 1) / (exp(_lnFragCombo)*exp(_lnFragCombo));
}


double BiologicalRepeatEstimator::NormalDistribution(const size_t k, const double mean, const double standardDeviation){
	const double PI = 3.141592653589793238463;
	double nd = exp(-0.5*((k - mean) / standardDeviation)*((k - mean) / standardDeviation)) / (standardDeviation * sqrt(2 * PI));
	return nd;
}


double BiologicalRepeatEstimator::LnPoissonDistribution(size_t k)
{
    return -_lambda + k * _lnLambda - GetLnFactoriel(k);
}

double BiologicalRepeatEstimator::LnDistribution(const size_t k)
{
	auto vs = CalculateVarianceShift(k);
	if (_lambda >= 100) return log(NormalDistribution(k, _lambda, sqrt(_lambda + vs)));
	else return LnPoissonDistribution(k);
}


double BiologicalRepeatEstimator::ForNonRepeats(size_t k, size_t last, int step)
{
	vector<size_t> factors;

    double totalExpected = 0.0;

	auto loopStart = (size_t)(k / _minimalFragmentProbability);
	if(step > 0) loopStart+=step;
	for (size_t num = loopStart; step > 0 ? num <= last : num >= last; num += step)
	{
        double lnBinom = LnBinomialDistribution(k, num);
		if (exp(lnBinom + _lnFragCombo) == 0) break;

		_factors->FillFactors (num, factors);

		for (auto factor1 : factors){
			size_t factor2 = num / factor1;

			if (factor1 < _lnDistribution.size() && factor2 < _lnDistribution.size())
            {
				auto lnbin = (alphabetSize_ == 4 && factor1 <= 16 && factor2 <= 16 && DnaNonRepeatCoefficients[factor1][factor2] != nullptr) ? log(DnaNonRepeatCoefficients[factor1][factor2][k]) : lnBinom;
				double expected = exp(_lnDistribution[factor1] + _lnDistribution[factor2] + lnbin + _lnFragCombo);
				auto help1 = exp(lnbin);
				auto help2 = exp(_lnDistribution[factor1]);
				auto help3 = exp(_lnDistribution[factor2]);
				if (factor1 != factor2) expected *= 2;

				totalExpected += expected;
            }
		}
    }

	return totalExpected;
}

double BiologicalRepeatEstimator::Compute(size_t k)
{
	size_t searchEnd = (_lnDistribution.size() - 1) * (_lnDistribution.size() - 1);

    double expected = ForNonRepeats(k, searchEnd, 1) + ForNonRepeats(k, k,-1);
	if (std::isnan(expected)) throw logic_error("Numeric computation failed.");
	return expected / 2;
}

