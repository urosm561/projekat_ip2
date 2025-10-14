#include <vector>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <cassert>
#include "MathematicalPalindromeEstimator.h"
#include "RepeatCoefficients.h"

using namespace std;

MathematicalPalindromeEstimator::MathematicalPalindromeEstimator(size_t seqLen, size_t fragLen, size_t alphabetSize, const Factors& factors)
{
	alphabetSize_ = alphabetSize;
	_factors = &factors;
	_seqLen = seqLen;
	_fragLen = fragLen;
	_fragCount = _seqLen - _fragLen + 1;

	auto palLen = fragLen / 2;

	_lnFragCombo = log(alphabetSize_) * _fragLen;

	auto checkedPow = pow(alphabetSize, palLen);
	auto fraction = std::isinf(checkedPow) ? 1 : (checkedPow - 1) / checkedPow;
	_lnFragComboNonPal = fragLen * log(alphabetSize) + log(fraction);

	_lnFragComboPal = (fragLen - palLen) * log(alphabetSize);

	_lnfaktoriel.push_back(0);
	_lnLambda = log(_fragCount) - _lnFragCombo;
	_lambda = exp(_lnLambda);
	
	_minimalFragmentProbability = double(alphabetSize_ - 1) * (alphabetSize_ - 1) / alphabetSize_ / alphabetSize_;
	_lnMinimalFragmentProbability = log(_minimalFragmentProbability);
	_lnOneMinusMfP = log(1 - _minimalFragmentProbability);

	_lnSelfMinimalFragmentProbability = log(double(alphabetSize - 1) / alphabetSize);
	_lnSelfOneMinusMfP = log(1 - double(alphabetSize - 1) / alphabetSize);

	_lnDistribution.reserve(_fragCount);

	for (size_t x = 0; x < _fragCount; ++x)
	{
		double pois = LnDistribution(x);
		if (x > _lambda && exp(pois + _lnFragCombo) == 0) break;
		_lnDistribution.push_back(pois);
	}
}

void MathematicalPalindromeEstimator::EnsureFactoriel(size_t k)
{
	if (k >= _lnfaktoriel.size())
	{
		for (size_t i = _lnfaktoriel.size(); i <= k; ++i)
		{
			_lnfaktoriel.push_back(_lnfaktoriel[i - 1] + log(i));
		}
	}
}

double MathematicalPalindromeEstimator::LnBinomialDistribution(size_t k, size_t n)
{
	assert(0 <= k && k <= n);

	EnsureFactoriel(n);

	double result =
		+_lnfaktoriel[n]
		- _lnfaktoriel[n - k]
		- _lnfaktoriel[k]
		+ k * _lnMinimalFragmentProbability
		+ (n - k) * _lnOneMinusMfP
		;

	return result;
}

double MathematicalPalindromeEstimator::LnBinomialSelfDistribution(size_t k, size_t n) {
	assert(0 <= k && k <= n);

	EnsureFactoriel(n);

	double result =
		+_lnfaktoriel[n]
		- _lnfaktoriel[n - k]
		- _lnfaktoriel[k]
		+ k * _lnSelfMinimalFragmentProbability
		+ (n - k) * _lnSelfOneMinusMfP
		;

	return result;
}

double MathematicalPalindromeEstimator::LnMatRepBinomialDistribution(size_t k, size_t words) {
	auto pairs = words * (words - 1) / 2;

	assert(0 <= k && k <= pairs + words);

	EnsureFactoriel(pairs);
	EnsureFactoriel(words);

	auto base = k > words ? k - words : 0;
	if (base > pairs) base = pairs;
	auto maxBase = k;
	if (maxBase > pairs) maxBase = pairs;

	auto ret = 0.0;

	for (auto x = base; x <= maxBase; ++x) {
		auto lnBinom1 = LnBinomialDistribution(base, pairs);
		auto lnBinom2 = LnBinomialSelfDistribution(k - base, words);
		ret += exp(lnBinom1 + lnBinom2);
	}

	return log(ret);
}

double MathematicalPalindromeEstimator::GetLnBinomRepeats(size_t k, size_t n, size_t comb, double p)
{
	if (alphabetSize_ < GenericMatPalCoefficients.size() && n < GenericMatPalCoefficients[alphabetSize_].size()) {
		auto arr = GenericMatPalCoefficients[alphabetSize_][n];

		if (arr != nullptr) {
			return log(arr[k]);
		}
	}

	return LnMatRepBinomialDistribution(k, n);
}

double MathematicalPalindromeEstimator::CalculateVarianceShift(const size_t k) {
	return (_seqLen - _fragLen + 1) / (exp(_lnFragCombo)) - ((2 * _fragLen - 1)*_seqLen - 3 * _fragLen*_fragLen + 4 * _fragLen - 1) / (exp(_lnFragCombo)*exp(_lnFragCombo));
}


double MathematicalPalindromeEstimator::NormalDistribution(const size_t k, const double mean, const double standardDeviation) {
	const double PI = 3.141592653589793238463;
	double nd = exp(-0.5*((k - mean) / standardDeviation)*((k - mean) / standardDeviation)) / (standardDeviation * sqrt(2 * PI));
	return nd;
}


double MathematicalPalindromeEstimator::LnPoissonDistribution(size_t k)
{
	EnsureFactoriel(k);
	return -_lambda + k * _lnLambda - _lnfaktoriel[k];
}

double MathematicalPalindromeEstimator::LnDistribution(const size_t k)
{
	auto vs = CalculateVarianceShift(k);
	if (_lambda >= 1000) return log(NormalDistribution(k, _lambda, sqrt(_lambda + vs)));
	else return LnPoissonDistribution(k);
}


double MathematicalPalindromeEstimator::ForNonRepeats(size_t k, size_t last, int step)
{
	vector<size_t> factors;

	double totalExpected = 0.0;

	auto loopStart = (size_t)(k / _minimalFragmentProbability);
	if (step > 0) loopStart += step;
	for (size_t num = loopStart; step > 0 ? num <= last : num >= last; num += step)
	{
		double lnBinom = LnBinomialDistribution(k, num);
		if (exp(lnBinom + _lnFragComboNonPal) == 0) break;

		_factors->FillFactors(num, factors);

		for (auto factor1 : factors) {
			size_t factor2 = num / factor1;

			if (factor1 < _lnDistribution.size() && factor2 < _lnDistribution.size())
			{
				auto lnbin = (alphabetSize_ == 4 && factor1 <= 16 && factor2 <= 16 && DnaNonRepeatCoefficients[factor1][factor2] != nullptr) ? log(DnaNonRepeatCoefficients[factor1][factor2][k]) : lnBinom;
				double expected = exp(_lnDistribution[factor1] + _lnDistribution[factor2] + lnbin + _lnFragComboNonPal);
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

double MathematicalPalindromeEstimator::ForRepeats(size_t k)
{
	size_t limit = 32;
	if (alphabetSize_ == 20) limit = 5;
	double totalExpected = 0.0;

	auto begin = size_t(floor((sqrt(1 + 8 * k) - 1) / 2 + 0.5));
	if (begin*(begin + 1) / 2 < k) begin++;

	auto max = _lnDistribution.size(); // max poisson distribution value
	for (size_t numWords = begin; numWords < max; ++numWords) {
		auto comb = numWords * (numWords + 1) / 2;
		if (comb == 0) {
			continue;
		}

		assert(k <= comb);

		double lnBinom = GetLnBinomRepeats(k, numWords, comb, _minimalFragmentProbability);
		if (std::isinf(lnBinom + _lnFragComboPal) && numWords > limit) {
			break;
		}

		auto x = lnBinom + _lnDistribution[numWords] + _lnFragComboPal;
		if (!std::isinf(x)) {
			totalExpected += exp(x);
		}
	}

	return totalExpected;
}

double MathematicalPalindromeEstimator::ComputeForNonPalindromes(size_t k)
{
	size_t searchEnd = (_lnDistribution.size() - 1) * (_lnDistribution.size() - 1);

	double expected = ForNonRepeats(k, searchEnd, 1) + ForNonRepeats(k, k, -1);
	if (std::isnan(expected)) throw logic_error("Numeric computation failed.");
	return expected / 2;
}

double MathematicalPalindromeEstimator::ComputeForPalindromes(size_t k)
{
	double expected = ForRepeats(k);
	if (std::isnan(expected)) throw logic_error("Numeric computation failed.");
	return expected;
}

