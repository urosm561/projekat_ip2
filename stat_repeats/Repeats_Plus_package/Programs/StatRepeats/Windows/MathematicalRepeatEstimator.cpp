#define _USE_MATH_DEFINES

#include <vector>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include "MathematicalRepeatEstimator.h"
#include "RepeatCoefficients.h"

using namespace std;

MathematicalRepeatEstimator::MathematicalRepeatEstimator(size_t seqLen, size_t fragLen, size_t alphabetSize, const Factors& factors)
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

double MathematicalRepeatEstimator::GetLnFactoriel(size_t k)
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

double MathematicalRepeatEstimator::LnBinomialDistribution(size_t k, size_t n)
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

double MathematicalRepeatEstimator::GetLnBinomRepeats(size_t k, size_t n, size_t comb, double p)
{
	if (alphabetSize_ < GenericRepeatCoefficients.size() && n < GenericRepeatCoefficients[alphabetSize_].size()) {
		auto arr = GenericRepeatCoefficients[alphabetSize_][n];

		if (arr != nullptr) {
			return log(arr[k]);
		}
	}

	return LnBinomialDistribution(k, comb);
}

double MathematicalRepeatEstimator::CalculateVarianceShift(const size_t k){
	return (_seqLen - _fragLen + 1) / (exp(_lnFragCombo)) - ((2*_fragLen-1)*_seqLen - 3*_fragLen*_fragLen + 4*_fragLen - 1) / (exp(_lnFragCombo)*exp(_lnFragCombo));
}


double MathematicalRepeatEstimator::NormalDistribution(const size_t k, const double mean, const double standardDeviation){
	const double PI = 3.141592653589793238463;
	double nd = exp(-0.5*((k - mean) / standardDeviation)*((k - mean) / standardDeviation)) / (standardDeviation * sqrt(2 * PI));
	return nd;
}


double MathematicalRepeatEstimator::LnPoissonDistribution(size_t k)
{
    return -_lambda + k * _lnLambda - GetLnFactoriel(k);
}

double MathematicalRepeatEstimator::LnDistribution(const size_t k)
{
	auto vs = CalculateVarianceShift(k);
	if (_lambda >= 100) return log(NormalDistribution(k, _lambda, sqrt(_lambda + vs)));
	else return LnPoissonDistribution(k);
}


double MathematicalRepeatEstimator::ForRepeats(size_t k, size_t last, int step)
{
	size_t limit = 30;
	if (alphabetSize_ == 20) limit = 5;
    double totalExpected = 0.0;

	auto begin = size_t(floor((1 + sqrt(1 + 8 * k)) / 2 + 0.5));
	if (begin*(begin-1)/2 < k) begin++;

	auto lnDistSize = _lnDistribution.size();

	for (size_t num = begin; num != last + step && num < lnDistSize; num += step)
    {
		auto comb = num*(num-1)/2;
		if (comb == 0) {
			continue;
		}

		double lnBinom = GetLnBinomRepeats(k, num, comb, _minimalFragmentProbability);
		if (std::isinf(lnBinom + _lnFragCombo) && num > limit) {
			break;
		}

		auto x = lnBinom + _lnDistribution[num] + _lnFragCombo;
		if (!std::isinf(x)) {
			totalExpected += exp(x);
		}
    }

	return totalExpected;
}

double MathematicalRepeatEstimator::Compute(size_t k)
{
	size_t searchEnd = _lnDistribution.size() - 1;
	
	double expected = ForRepeats(k, searchEnd, 1);
	if (std::isnan(expected)) throw logic_error("Numeric computation failed.");
	return expected;
}

