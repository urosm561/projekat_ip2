#include <utility>
#include <stack>
#include <vector>
#include <ostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <list>
#include <memory>
#include <array>
#include <unordered_set>
#include <unordered_map>
#include <locale>
#include <algorithm>
#include <math.h>
#include <functional>
#include <ctime>
#include <string>
#include <tuple>
#include <cstdio>
#include <climits>
#include <stdexcept>
//#include "Db2EmbeddedSqlRepeatPrinter.h"
#include <regex>
#include <limits>
#include "divsufsort.h"
#include "CommandLineParser.h"
#include "FastaReader.h"
#include "CommandLine.h"
#include "Primes.h"
#include "Factors.h"
#ifdef ODBCOUTPUT
#include "CppOdbc.h"
#pragma comment(lib, "odbc32.lib")
#endif
#include "ConfidenceBound.h"
#include "RepeatCoefficients.h"
#include "SharedLifetimeAllocator.h"
#include "StringSlice.h"
#include "RepeatPrinter.h"
#include "MathematicalPalindromeEstimator.h"
#include "BiologicalRepeatEstimator.h"
#include "MathematicalRepeatEstimator.h"
#include "BiologicalPalindromeEstimator.h"
#include "RepeatPrinter.h"


using namespace std::rel_ops;
using namespace std;

struct NakedLcpInterval
{
	size_t lcp;
	size_t lb;
	size_t rb;
};

struct LcpInterval : public NakedLcpInterval
{
	bool defined = false;
	vector<NakedLcpInterval> children;
};

LcpInterval CreateLcpInterval(int lcp, int lb, int rb)
{
	LcpInterval x;
	x.defined = true;
	x.lcp = lcp;
	x.lb = lb;
	x.rb = rb;
	return x;
}

struct SearchType
{
	virtual void UpdateRepeatLocations(size_t& x, size_t& y, size_t length) const = 0;
	virtual bool IsResult(size_t x, size_t y, size_t length, const regex& motif) const = 0;
	virtual StringSlice GetComplement(const StringSlice& s) const = 0;
	virtual bool IncludePalindromes() const = 0;
};

void Union(array<vector<size_t>, 256> &results1, array<vector<size_t>, 256> &results2)
{
	for (size_t i = 0; i < results1.size(); i++) {
		results2[i].insert(results2[i].begin(), results1[i].begin(), results1[i].end());
		results1[i].clear();
	}
}

static size_t absDif(size_t x, size_t y) {
	return x > y ? x - y : y - x;
}

template <typename Func>
void PrintRepeats(const StringSlice& repeat, size_t x, size_t y, size_t h, SearchType* searchType,
	Func forEachRepeat, size_t maxGap, const regex& motif)
{
	if (!searchType->IsResult(x, y, h, motif)) return;


	searchType->UpdateRepeatLocations(x, y, h);

	if (x != y) {
		if (searchType->GetComplement(repeat) < repeat) {
			swap(x, y);
		}
		if (absDif(h, absDif(x, y)) > maxGap) return;
		forEachRepeat(x, y, h);
	}
	else {
		if (searchType->IncludePalindromes()) {
			if (absDif(h, absDif(x, y)) > maxGap) return;
			forEachRepeat(x, y, h);
		}
	}
}

void ToUpper(string& s)
{
	for (auto it = s.begin(); it != s.end(); ++it) {
		*it = toupper(*it);
	}
}

template <typename Func>
void Combine(const StringSlice& repeat,
	const array<vector<size_t>, 256> &results,
	size_t h,
	SearchType* searchType,
	Func& forEachRepeat,
	size_t maxGap,
	const regex& motif)
{
	for (size_t i = 0; i < results.size(); i++) {
		if (results[i].empty()) continue;

		for (size_t k = i + 1; k < results.size(); k++) {
			if (results[k].empty()) continue;

			for (auto x : results[i]) {
				for (auto y : results[k]) {
					PrintRepeats(repeat, x, y, h, searchType, forEachRepeat, maxGap, motif);
				}
			}
		}
	}
}

bool Empty(const array<vector<size_t>, 256>& results) {
	for (auto& x : results) {
		if (!x.empty()) {
			return false;
		}
	}
	return true;
}

template <typename Func>
inline void CombineTwoResults(const StringSlice& repeat,
	const array<vector<size_t>, 256> &results1,
	const array<vector<size_t>, 256> &results2,
	size_t h,
	SearchType* searchType,
	Func& forEachRepeat,
	size_t maxGap,
	const regex& motif)
{
	for (size_t i = 0; i < results1.size(); i++) {
		if (results1[i].empty()) continue;

		for (size_t k = 0; k < results2.size(); k++) {
			if (results2[k].empty() || k == i) continue;

			for (auto x : results1[i]) {
				for (auto y : results2[k]) {
					PrintRepeats(repeat, x, y, h, searchType, forEachRepeat, maxGap, motif);
				}
			}
		}
	}
}

static size_t last_processed_size_ = 100;

template <typename Func>
void ComputeResults(
	const LcpInterval &l,
	const vector<size_t> &suftab,
	const string &text,
	size_t fragmentLength,
	array<vector<size_t>, 256>& results,
	array<vector<size_t>, 256>& lastResults,
	SearchType* searchType,
	Func& forEachRepeat,
	size_t maxGap,
	const regex& motif)
{
	vector<size_t> processed;
	processed.reserve(last_processed_size_);

	for (size_t i = 0; i < l.children.size(); i++) {
		if (l.lcp >= fragmentLength) {
			for (auto k = l.children[i].lb; k <= l.children[i].rb; k++) {
				processed.push_back(k);
				if (suftab[k] == 0) continue;
				results[static_cast<unsigned char> (text[suftab[k] - 1])].push_back(suftab[k]);
			}
			if (!Empty(lastResults)) {
				CombineTwoResults(StringSlice(text, suftab[l.lb], l.lcp), results, lastResults, l.lcp, searchType, forEachRepeat, maxGap, motif);
			}

			Union(results, lastResults);
		}
	}

	for (size_t x = 1; x < processed.size(); ++x) {
		if (processed[x - 1] >= processed[x]) throw logic_error("Processed not in ascending order.");
	}

	if (l.lcp >= fragmentLength) {
		for (auto k = l.lb; k <= l.rb; k++) {
			if (binary_search(processed.begin(), processed.end(), k)) continue;
			if (suftab[k] == 0) continue;
			results[static_cast<unsigned char> (text[suftab[k] - 1])].push_back(suftab[k]);
		}
		Combine(StringSlice(text, suftab[l.lb], l.lcp), results, l.lcp, searchType, forEachRepeat, maxGap, motif);
		CombineTwoResults(StringSlice(text, suftab[l.lb], l.lcp), results, lastResults, l.lcp, searchType, forEachRepeat, maxGap, motif);
	}

	for (auto& x : results) x.clear();
	for (auto& x : lastResults) x.clear();

	last_processed_size_ = processed.size();
}

template <typename Func>
void ComputeLcpIntervals(
	const vector<size_t> &lcp,
	const vector<size_t> &suftab,
	const string &text,
	const size_t fragmentLength,
	SearchType* searchType,
	Func forEachRepeat,
	unsigned maxGap,
	string motif)
{
	regex motifRegex{ motif, regex_constants::ECMAScript | regex_constants::optimize };

	array<vector<size_t>, 256> results, lastResults;

	LcpInterval lcpInterval = CreateLcpInterval(0, 0, -1);
	stack<LcpInterval, vector<LcpInterval>> lcpIntervals;
	lcpIntervals.push(lcpInterval);
	LcpInterval lastInterval;
	lastInterval.defined = false;
	size_t lb = 0;
	for (size_t i = 1; i < lcp.size(); i++) {
		lb = i - 1;
		while (lcp[i] < lcpIntervals.top().lcp) {
			lcpIntervals.top().rb = i - 1;
			lastInterval = lcpIntervals.top();
			lcpIntervals.pop();

			ComputeResults(lastInterval, suftab, text, fragmentLength, results, lastResults, searchType, forEachRepeat, maxGap, motifRegex);

			lb = lastInterval.lb;
			if (lcp[i] <= lcpIntervals.top().lcp) {
				lcpIntervals.top().children.push_back(lastInterval);
				lastInterval.defined = false;
			}
		}

		if (lcp[i] > lcpIntervals.top().lcp) {
			if (lastInterval.defined) {
				lcpInterval.lcp = lcp[i];
				lcpInterval.lb = lb;
				lcpInterval.rb = -1;
				lcpIntervals.push(lcpInterval);
				lcpIntervals.top().children.push_back(lastInterval);
				lastInterval.defined = false;
			}
			else {
				lcpInterval.lcp = lcp[i];
				lcpInterval.lb = lb;
				lcpInterval.rb = -1;
				lcpIntervals.push(lcpInterval);
			}
		}
	}

	lcpIntervals.top().rb = (int)lcp.size() - 1;
	lastInterval = lcpIntervals.top();
	lcpIntervals.pop();
}

vector<size_t> MakeSubInv(const vector<size_t>& suftab)
{
	vector<size_t> sufinv(suftab.size());
	for (size_t i = 0; i < suftab.size(); i++) {
		sufinv[suftab[i]] = i;
	}
	return sufinv;
}

bool LettersEqual(const string& text, size_t pos1, size_t pos2) {
	auto c1 = text[pos1];
	auto c2 = text[pos2];

	if (c1 || c2) {
		return c1 == c2;
	}
	else {
		return pos1 == pos2;
	}
}

vector<size_t> MakeLCP(const string &text, vector<size_t>& suftab, vector<size_t>& sufinv)
{
	vector<size_t> lcp(suftab.size());
	int l = 0;

	for (size_t i = 0; i <= suftab.size() - 2; i++) {
		auto k = sufinv[i];
		if (k == 0) continue;
		auto j = suftab[k - 1];
		while (LettersEqual(text, i + l, j + l)) l++;
		lcp[k] = l;
		if (l > 0) l--;
	}

	return lcp;
}


class MathematicalRepeat : public SearchType
{
private:
	const string& _text;

public:
	MathematicalRepeat(const string& text) : _text(text)
	{
	};

	void UpdateRepeatLocations(size_t& x, size_t& y, size_t length) const
	{
		if (y < x) {
			swap(x, y);
		}
	}

	bool IsResult(size_t x, size_t, size_t length, const regex& r) const
	{
		return regex_search(_text.begin() + x, _text.begin() + x + length, r);
	}

	StringSlice GetComplement(const StringSlice& slice) const
	{
		return slice;
	}

	bool IncludePalindromes() const
	{
		return false;
	}
};

class MathematicalPalindrome : public SearchType
{
private:
	size_t _median;
	const string& _text;

public:
	MathematicalPalindrome(size_t median, const string &text) : _text(text)
	{
		_median = median;
	};

	void UpdateRepeatLocations(size_t& x, size_t& y, size_t length) const
	{
		if (y > _median) {
			y = _median - (y - _median) - length + 1;
		}
		else {
			x = _median - (x - _median) - length + 1;
		}

		if (x > y) {
			swap(x, y);
		}
	}

	bool IsResult(size_t x, size_t y, size_t length, const regex& r) const
	{
		if ((y < _median) && (x > _median) && (_median + 1 >= y + length + (x - _median))) {
			return
				regex_search(_text.begin() + x, _text.begin() + x + length, r) ||
				regex_search(_text.begin() + y, _text.begin() + y + length, r)
				;
		}
		else if (((x < _median) && (y > _median) && (_median + 1 >= x + length + (y - _median))) || ((y < _median) && (x > _median) && (_median + 1 >= y + length + (x - _median)))) {
			return
				regex_search(_text.begin() + x, _text.begin() + x + length, r) ||
				regex_search(_text.begin() + y, _text.begin() + y + length, r)
				;
		}
		else {
			return false;
		}

	}

	StringSlice GetComplement(const StringSlice& slice) const
	{
		auto offset = slice.begin() - _text.begin();
		auto length = slice.size();
		auto first = _text.end() - offset - length;
		StringSlice s(first, first + length);
		return s;
	}

	bool IncludePalindromes() const
	{
		return true;
	}
};

class BiologicalPalindrome : public SearchType
{
private:
	size_t _median;
	const string& _text;
public:
	BiologicalPalindrome(size_t median, const string &text) : _text(text)
	{
		_median = median;
	}

	void UpdateRepeatLocations(size_t& x, size_t& y, size_t length) const
	{
		if (y > _median) {
			y = _median - (y - _median) - length + 1;
		}
		else {
			x = _median - (x - _median) - length + 1;
		}

		if (x > y) {
			swap(x, y);
		}
	}

	bool IsResult(size_t x, size_t y, size_t length, const regex& r) const
	{
		if (((x < _median) && (y > _median) && (_median + 1 >= x + length + (y - _median))) || ((y < _median) && (x > _median) && (_median + 1 >= y + length + (x - _median)))) {
			return
				regex_search(_text.begin() + x, _text.begin() + x + length, r) ||
				regex_search(_text.begin() + y, _text.begin() + y + length, r)
				;
		}
		else {
			return false;
		}

	}

	StringSlice GetComplement(const StringSlice& slice) const
	{
		auto offset = slice.begin() - _text.begin();
		auto length = slice.size();
		auto first = _text.end() - offset - length;
		StringSlice s(first, first + length);
		return s;
	}

	bool IncludePalindromes() const
	{
		return true;
	}
};

class BiologicalRepeat : public SearchType
{
private:
	size_t _median;
	const string& _text;
public:
	BiologicalRepeat(size_t median, const string &text) : _text(text)
	{
		_median = median;
	}

	void UpdateRepeatLocations(size_t& x, size_t& y, size_t length) const
	{
		if (y > x) {
			y = y - _median;
		}
		else {
			x = x - _median;
		}

		if (x > y) {
			swap(x, y);
		}
	}

	bool IsResult(size_t x, size_t y, size_t length, const regex& r) const
	{
		if (((x < _median) && (y > _median) && (y - _median > x)) || ((y < _median) && (x > _median) && (x - _median > y))) {
			return
				regex_search(_text.begin() + x, _text.begin() + x + length, r) ||
				regex_search(_text.begin() + y, _text.begin() + y + length, r)
				;
		}
		else return false;
	}

	StringSlice GetComplement(const StringSlice& slice) const
	{
		if (slice.begin() > _text.begin() + _median) {
			return StringSlice(slice.begin() - _median, slice.end() - _median);
		}
		else {
			return StringSlice(slice.begin() + _median, slice.end() + _median);
		}
	}

	bool IncludePalindromes() const
	{
		return false;
	}
};

static const char* Instructions() {
	return 		"Usage:\n"
		"StatRepeats <input file name> <fragment length> [<options>]\n"
		"\n"
		" \n"

		"Mandatory options are:\n"
		"inputFileName specifies the input file name, does not have a default value, must be specified.\n"
		"fragment length is a number that specifies the minimal fragment length, where fragments are words that are results. Does not have default value, must be specified. If zero is specified StatRepeats automatically chooses the lowest possible value for witch the probability estimation is valid.\n"
		" \n"
		" \n"
		"Other options are:\n"
		"-dn for direct non-complementary repeats, -dc for direct complementary repeats, -in for inverse non-complementary reapeats and -ic for inverse complementary repeats. Default is dn. \n"
		"-rna, for working with rna sequence, default is dna sequence with alphabet containing a,c,t,g. Alphabet size in rna is four and contains a,c,u,g.\n"
		"-protein, for working with protein sequence, default is dna sequence with alphabet containing a,c,t,g. With option proteins alphabet size is between 20 and 22, depending on the input sequence\n"
		"When working with proteins user can specify mappings for specific amino acids into predefined groups. Available groups are: aliphatic, sulphur, tiny, aromatic, hydrophobic, charged, positive, polar, acidic, small and hydroxylic\n"
		"-maxgap number, where this is maximal gap between two results, default is finding all results, can be abbreviated as -max[gap] number. \n"
		"-noprobability, for working without probability estimation. Default is working with probability estimation.Can be abbreviated as -n[otprobability]. \n"
		"-pvalue number, specifies p value, default is 0.05. Can be abbreviated as -pv[alue] number. \n"
		"-msr. Finds results in all input files at once. In this case probability estimation is unabled. Can be abbreviated as -ms[r]. \n"
		"-complement filename. Reads complement rules from specified file. Can be abbreviated as -compl[ement] filename.\n"
		"-exclude number, groups of N consecutive letters that are larger or equal then given number are excluded from input sequence if option letter is not specified. If option -letter with given letter is specified then consecutive letters that are larger or equal then given number are excluded from input sequence. Can be abbreviated as -ex[clude]\n"
		"-motif string. The string specifies the motif that must be contained in all results.\n"
		" \n"
		" \n"
		"Output options:\n"
		"Default for output is writing results on console.\n"
		"-output outputFileName, can be abbreviated as -o[utput] outputFileName.\n"
		"-split, option for making three output files. Sequence name is written in file named sequencename.id, statistics are written in file nemed sequencename.stat, and results are written in file named sequencename.load. Can be abbreviated as -sp[lit].\n"
		"-load name,option for making three output files. Sequence name is written in file named name.id, statistics are written in file nemed name.stat, and results are written in file named name.load,can be abbreviated as -loa[d] name. \n"
#ifdef ODBCOUTPUT
		"-database string, specifies database connection string, writes results to database, can be abbreviated as -dat[abase] string.\n"
		"-db2output string, specifies database connection string for db2, writes results to tables Sequence, Fragment and Match. Can be abbreviated as -db2[output] string.\n"
#endif
		"Only one output option can be chosen. Output options are outputfile,split,load,database."
		" \n"
		" \n"
#ifdef ODBCOUTPUT
		"-logging true/false. Default is true. If false alters table fragment and match activate not logged initially,can be abbreviated as -log[ging] true/false. Used only with DB2.\n"
		"-commitcount number, number of inserts in one transaction, default is 0 and then there will be only one commit at the end of the program, can be abbreviated as -com[mitcount] number.\n"
#endif
		"-printinstances true/false, false for not printing instances, default is true. Can be abbreviated as -pri[ntinstances] true/false. \n"
		" \n"

		"Some examples of usage:\n"
		"SmartRepeats sars.fasta 6 -in -printinstances false -exclude 33 -letter A \n"
		"SmartRepeats sars.fasta 6 -in -motif T.C[CT]AG[^A]G \n"
		"SmartRepeats sars.fasta 5 -dc -rna -noprobability -complement rules.txt -output rez.txt\n"
		"SmartRepeats sars.fasta 6 -in -conf 0.99 -load true\n"
		"SmartRepeats sars.fasta 6 -in -msr \n"
#ifdef ODBCOUTPUT
		"SmartRepeats sars.fasta 6 -ic -db2output connectionstring \n"
#endif
		"SmartRepeats input.fasta 5 -in -proteins -charged -acidic \n"
		"SmartRepeats sars.fasta 6 -in -split -proteins -noprobability \n"
#ifdef ODBCOUTPUT
#endif
		"\n"
		"\n";
}

vector<tuple<size_t, size_t, size_t>> CategorizeStrings(const shared_lifetime_unordered_map<StringSlice, size_t>& strings) {
	shared_lifetime_unordered_map<pair<size_t, size_t>, size_t> occurence(1000, hash<pair<size_t, size_t>>(), equal_to<pair<size_t, size_t>>(), shared_lifetime_allocator<pair<const pair<size_t, size_t>, size_t>>(true));
	for (const auto& x : strings) {
		auto len = x.first.size();
		auto repeatCount = x.second;

		auto key = make_pair(len, repeatCount);
		auto it = occurence.find(key);
		if (it == occurence.end()) {
			occurence.insert(make_pair(key, 1));
		}
		else {
			++(it->second);
		}
	}

	vector<tuple<size_t, size_t, size_t>> sorted;
	sorted.reserve(occurence.size());
	for (auto & x : occurence) sorted.push_back(make_tuple(x.first.first, x.first.second, x.second));

	sort(sorted.begin(), sorted.end(),
		[](const tuple<size_t, size_t, size_t>& x, const tuple<size_t, size_t, size_t>& y) {
		if (get<0>(x) < get<0>(y)) return true;
		if (get<0>(x) == get<0>(y)) return get<1>(x) < get<1>(y);
		return false;
	});

	occurence.get_allocator().fast_free_on();

	return sorted;
}

template<typename Estimator>
void ComputeStatistics(size_t sequenceLength, const shared_lifetime_unordered_map<StringSlice, size_t>& strings, size_t alphabetSize, double pval, shared_lifetime_unordered_set<StringSlice>& results, const size_t multiplier, const size_t motifSize) {
	if (strings.empty()) return;

	auto sorted = CategorizeStrings(strings);

	shared_lifetime_unordered_set<pair<size_t, size_t>> passedStatistics(1000, hash<pair<size_t, size_t>>(), equal_to<pair<size_t, size_t>>(), shared_lifetime_allocator<pair<size_t, size_t>>(true));

	Primes primes;
	Factors factors(primes);
	size_t estimatorFragmentLen = get<0>(sorted[0]);

	Estimator estimator(sequenceLength, estimatorFragmentLen, alphabetSize, factors);

	for (const auto& x : sorted) {
		if (get<0>(x) != estimatorFragmentLen) {
			estimatorFragmentLen = get<0>(x);
			estimator = Estimator(sequenceLength, estimatorFragmentLen, alphabetSize, factors);
		}

		auto a = estimator.Compute(get<1>(x));
		if (motifSize > 0) a = exp(log(a)+log(multiplier)-motifSize*log(alphabetSize));
		auto confidenceBound = GetConfidenceBound(a, 1 - pval);
		
		//cout << get<0>(x) << ", " << get<1>(x) << ", " << get<2>(x) << ", ocekivano, " << a << ", confidence, " << confidenceBound;

		if (get<2>(x) >= confidenceBound) {
			passedStatistics.insert(make_pair(get<0>(x), get<1>(x)));
			//cout << ", prosao";
		}
		//cout << endl;
	}

	for (const auto & x : strings) {
		if (passedStatistics.find(make_pair(x.first.size(), x.second)) != passedStatistics.end()) {
			results.insert(x.first);
		}
	}

	passedStatistics.get_allocator().fast_free_on();
}

RepeatPrinter* SetOutput(CommandLine& c, bool& writeInfo, size_t& outputCounter) {
	RepeatPrinter* repeatPrinter = nullptr;

	if (c.GetOutputFileName() != "") {
		repeatPrinter = new FileOutput(c.GetOutputFileName());
		outputCounter++;
		writeInfo = false;
	}
	else if (c.GetSplit()) {
		repeatPrinter = new SplitOutput();
		outputCounter++;
		writeInfo = false;
	}
	else if (c.GetLoad() != "") {
		repeatPrinter = new LoadOutput(c.GetLoad());
		outputCounter++;
		writeInfo = false;
	}
#ifdef ODBCOUTPUT
	else if (c.GetDatabase() != "") {
		repeatPrinter = CreateDataBaseOutput(c);
		outputCounter++;
		writeInfo = false;
	}
	/*else if (c.GetDb2Output() != "") {
		repeatPrinter = CreateDB2Output(c);
		outputCounter++;
		writeInfo = false;
	}*/
#endif // !NOODBC
	else if (c.GetComputeStatistics()) {
		repeatPrinter = new StatisticsOutput();
		outputCounter++;
	}
#ifdef DB2EMBEDDEDOUTPUT
	else if (c.GetStaticSql() != "") {
		repeatPrinter = CreateDb2EmbeddedSqlRepeatPrinter(c);
		outputCounter++;
		writeInfo = false;

	}
#endif
	else {
		repeatPrinter = new ConsoleOutput();
	}

	if (outputCounter > 1) {
		throw invalid_argument("There can only be one output option. Output options are outputfile,split,load,database,db2output,staticsql.");
	}
	return repeatPrinter;
}

pair<string, int> ParseSequenceName(const string& firstLine) {
	int version;
	string sequenceName;
	string afterDot = "";
	auto startswith = firstLine.find_first_of("|");
	auto firststring = firstLine.substr(0, startswith);
	if (firststring == "gi") {
		auto pos2 = firstLine.find_last_of("|");
		auto pos1 = firstLine.find_last_of("|", pos2 - 1);
		sequenceName = firstLine.substr(pos1 + 1, pos2 - pos1 - 1);
		auto position = sequenceName.find(".");
		if (position != string::npos) {
			afterDot = sequenceName.substr(position + 1, sequenceName.length() - position - 1);
			sequenceName = sequenceName.substr(0, position);
		}
	}
	else if (firststring == "DisProt") {
		sequenceName = firstLine.substr(startswith + 1, firstLine.find_first_of("|", startswith + 1) - startswith - 1);
	}
	else if (firststring == "sp" || firststring == "tr" || firststring == "gb") {
		auto first = firstLine.find_first_of("|");
		auto second = firstLine.find_first_of("|", first + 1);
		sequenceName = firstLine.substr(first + 1, second - first - 1);
	}
	else if (firststring == "pir") {
		sequenceName = firstLine.substr(firstLine.find_first_of("||") + 2, firstLine.find_first_of(" ", 1) - firstLine.find_first_of("||") - 2);
	}
	else if (firststring == "prf") {
		sequenceName = firstLine.substr(firstLine.find_first_of("||") + 2, firstLine.find_first_of(" ", 1) - firstLine.find_first_of("||") - 2);
	}
	else {
		//sequenceName = firstLine;
		sequenceName = firstLine.substr(0, firstLine.find_first_of(" ", 1));
	}

	version = atoi(afterDot.c_str());
	if (afterDot != "" && version == 0) {
		sequenceName += "." + afterDot;
	}

	if (afterDot == "") {
		version = 0;
	}

	return make_pair(sequenceName, version);
}


pair<size_t, string> CountNonAlphabetMakeComplement(string& x, const char(&complement)[256], bool(&nonAlphabet)[256])
{
	size_t count = 0;

	string sequence = x;
	for (size_t i = 0; i < sequence.length() - 1; i++) {
		auto index = static_cast<unsigned char> (x[i]);
		if (complement[index] == 0 && index != '\0') {
			nonAlphabet[index] = true;
			count++;
		}
		else sequence[i] = complement[index];
	}
	return make_pair(count, sequence);
}

void MakeCompRules(const CommandLine& c, char(&complement)[256], size_t& alphabetSize) {
	string fileName = c.GetComplementFileName();
	if (fileName != "") {
		ifstream file(fileName);
		if (!file) {
			string errorMsg = strerror(errno);
			throw runtime_error(errorMsg + " while opening the file " + fileName);
		}

		while (true) {
			char letter;
			do
			{
				file.get(letter);
			} while (isspace(letter) && !file.eof());
			if (file.eof()) break;

			char compLetter;
			do
			{
				file.get(compLetter);
			} while (isspace(compLetter) && !file.eof());
			if (file.eof()) throw invalid_argument("There must be an even number of letters in the complement file.");;

			letter = toupper(letter);
			compLetter = toupper(compLetter);

			complement[static_cast<unsigned char> (letter)] = compLetter;
			alphabetSize++;
		}
	}
	else if (c.GetProteins()) {
		alphabetSize = 20;
		complement['A'] = 'A';
		complement['C'] = 'C';
		complement['T'] = 'T';
		complement['G'] = 'G';
		complement['D'] = 'D';
		complement['E'] = 'E';
		complement['F'] = 'F';
		complement['H'] = 'H';
		complement['I'] = 'I';
		complement['K'] = 'K';
		complement['L'] = 'L';
		complement['M'] = 'M';
		complement['N'] = 'N';
		complement['P'] = 'P';
		complement['Q'] = 'Q';
		complement['R'] = 'R';
		complement['S'] = 'S';
		complement['V'] = 'V';
		complement['W'] = 'W';
		complement['Y'] = 'Y';
	}
	else if (c.GetRna()) {
		alphabetSize = 4;
		complement['A'] = 'U';
		complement['U'] = 'A';
		complement['C'] = 'G';
		complement['G'] = 'C';
	}
	else {
		alphabetSize = 4;
		complement['A'] = 'T';
		complement['T'] = 'A';
		complement['C'] = 'G';
		complement['G'] = 'C';
	}
}

size_t SuggestLength(size_t sequenceLength, const size_t alphabetSize, const size_t t) {
	size_t length = 0;
	if (alphabetSize >= 4 && alphabetSize <= 6) {
		if (t != 1) {
			if (sequenceLength <= 500) length = 3;
			if (sequenceLength > 500) length = 4;
			if (sequenceLength > 3000) length = 5;
			if (sequenceLength > 30000) length = 6;
			if (sequenceLength > 231000) length = 7;
			if (sequenceLength > 400000) length = 8;
			if (sequenceLength > 1500000) length = 9;
			if (sequenceLength > 5000000) length = 10;
			if (sequenceLength > 10000000) length = 11;
		}
		else {
			if (sequenceLength <= 500) length = 3;
			if (sequenceLength > 500) length = 4;
			if (sequenceLength > 3000) length = 5;
			if (sequenceLength > 27000) length = 6;
			if (sequenceLength > 231000) length = 7;
			if (sequenceLength > 400000) length = 8;
			if (sequenceLength > 1500000) length = 9;
			if (sequenceLength > 5000000) length = 10;
			if (sequenceLength > 10000000) length = 11;
		}
		return length;
	}
	if (sequenceLength <= 1000) {
		if (alphabetSize == 7 || alphabetSize == 8) length = 3;
		else length = 2;
		return length;
	}

	if (sequenceLength <= 2000) {
		if (alphabetSize >= 7 && alphabetSize <= 11) length = 3;
		else length = 2;
		return length;
	}

	if (sequenceLength <= 9000) {
		length = 3;
		return length;
	}

	if (sequenceLength >= 9000 && sequenceLength <= 30000) {
		if (alphabetSize >= 7 && alphabetSize <= 10) length = 4;
		else length = 3;
		return length;
	}

	if (sequenceLength > 30000) length = 6;
	if (sequenceLength > 231000) length = 7;
	if (sequenceLength > 400000) length = 8;
	if (sequenceLength > 1500000) length = 9;
	if (sequenceLength > 5000000) length = 10;
	if (sequenceLength > 10000000) length = 11;

	return length;
}

bool IsInAlphabet(const size_t letter) {
	return true;
}

size_t ExcludeN(string& sequence, const CommandLine& c) {
	size_t excludedCount = 0;
	char letter = c.GetExcludeLetter()[0];
	if (c.GetExcludeLetter().size() > 1) {
		throw invalid_argument("You may specify only one letter to exclude from input sequence.");
	}
	auto nstart = sequence.end();
	if (c.GetExclude() != -1) {
		for (auto it = sequence.begin(); it != sequence.end(); ++it) {
			if (*it == letter) {
				if (nstart == sequence.end()) {
					nstart = it;
				}
			}
			else {
				if (nstart != sequence.end()) {
					auto distance = it - nstart;
					if (distance >= c.GetExclude()) {
						fill(nstart, it, '\0');
						excludedCount = excludedCount + distance;
					}
					nstart = sequence.end();
				}
			}
		}
		if (nstart != sequence.end()) {
			auto distance = sequence.end() - nstart;
			if (distance >= c.GetExclude()) {
				fill(nstart, sequence.end(), '\0');
				excludedCount = excludedCount + distance;
			}
		}
	}
	return excludedCount;
}

void SetForProteins(string& sequence, CommandLine& c, size_t& alphabetSize, char(&complement)[256]) {
	alphabetSize = alphabetSize + c.GetAlphabetReduction();
	auto it = sequence.find_first_of('U');
	if (it != string::npos) {
		alphabetSize++;
		complement['U'] = 'U';
	}
	it = sequence.find_first_of('O');
	if (it != string::npos) {
		alphabetSize++;
		complement['O'] = 'O';
	}
	if (c.GetAlphabetReduction() != 0) {
		const auto& mapping = c.GetMapping();
		for (size_t i = 0; i != sequence.size(); i++) {
			auto a = static_cast<unsigned char>(sequence[i]);
			if (a != mapping[a]) sequence[i] = mapping[a];
			complement[mapping[a]] = mapping[a];
		}
	}
}

string NonAlphabetLetters(bool(&nonAlphabet)[256], const size_t excludedCount) {
	{
		string letters;
		for (size_t i = 0; i < 256; i++) {
			if (nonAlphabet[i]) {
				letters += char(i);
			}
		}
		if (excludedCount) {
			letters += "N";
		}
		return letters;
	}

}

bool FindSuggestLength(CommandLine& c, size_t seqLength, const size_t alphabetSize, const size_t t, string& sequenceName) {

	auto suggestedlength = SuggestLength(seqLength, alphabetSize, t);

	if (c.GetFragmentLength() == 0) {
		if (suggestedlength == 0) {
			cerr << "Minimal fragment length not specified or invalid. Cannot suggest a good minimal fragment length for alphabets that don't have cardinality 4 or 20." << endl;
		}
		else {
			c.SetFragmentLength(suggestedlength);
		}
	}
	else if (suggestedlength > c.GetFragmentLength()) {
		cout << "For sequence " << sequenceName << " with length " << seqLength << " and alphabet cardinality " << alphabetSize << " we suggest using results starting from minimal length " << suggestedlength << "." << endl;
		return false;
	}
	return true;
}


string MotifToRegex(const string& motif, size_t& multiplier, const size_t alphabetSize, size_t& motifSize, char(&complement)[256]) {
	multiplier = 1;
	size_t i = 0;
	string regex = "";
	auto m = motif;
	ToUpper(m);
	while (i < motif.size()) {
		if (motif[i] == '.') {
			regex += '.';
			motifSize++;;
		}
		else if (isalpha(motif[i])) {
			regex += m[i];
			motifSize++;;
		}
		else if (motif[i] == '[' && motif[i + 1] != '^') {
			regex += m[i];
			size_t j = 1;
			size_t k = 0;
			while (i + j < motif.size() && motif[i + j] != ']') {
				if (!isalpha(motif[i + j])) throw invalid_argument("Motif mask is not correct.");
				regex += m[i + j];
				if (complement[m[i+j]] != 0) k++;
				j++;
			}
			if (j == 1) throw invalid_argument("Motif mask is not correct.");
			i = i + j;
			regex += m[i];
			if (k!=0) multiplier = multiplier * k; 
			motifSize++;
		}
		else if (motif[i] == '[' && motif[i + 1] == '^') {
			regex += "[^";
			size_t j = 2;
			size_t k = 0;
			while (i + j < motif.size() && motif[i + j] != ']') {
				if (!isalpha(motif[i + j])) throw invalid_argument("Motif mask is not correct.");
				regex += m[i + j];
				if (complement[m[i+j]] != 0) { 
					k++;
				}
				j++;
			}
			if (j == 1) throw invalid_argument("Motif mask is not correct.");
			i = i + j;
			regex += "]";
			if (k!=0) multiplier = multiplier * (alphabetSize - k);
			motifSize++;
		}
		else throw invalid_argument("Motif mask is not correct.");
		i++;
	}
	return regex;
}


void MainImplementation(int argc, char** argv)
{
	cout << "StatRepeats v1.r6 Copyright ( C ) 2020 Ana Jelovic. All Rights Reserved." << endl;

	CommandLine c = CommandLine(argc, argv);

	size_t outputCounter = 0;
	bool writeInfo = true;
	auto repeatPrinter = unique_ptr<RepeatPrinter>(SetOutput(c, writeInfo, outputCounter));

	char complement[256] = {};
	size_t alphabetSize = 0;

	MakeCompRules(c, complement, alphabetSize);

	size_t multiplier = 1;
	size_t motifSize = 0;
	string reg = "";
	if (c.GetMotif() != "") reg = MotifToRegex(c.GetMotif(), multiplier, alphabetSize, motifSize, complement);

	repeatPrinter->Initialize(c, alphabetSize);
	size_t firstAlphabetSize = alphabetSize;
	bool first = true;

	if (writeInfo) cout << "Reading file..." << endl;
	auto start = clock();
	vector<string> sequenceStrings;
	for (auto &seq : FastaFileView(c.GetInputFileName(), !c.CombinedSequenceName().empty())) {
		auto firstLine = ParseSequenceName(c.CombinedSequenceName().empty() ? static_cast<const string&> (seq.NamesAndOffsets[0].first) : c.CombinedSequenceName());
		string sequenceName = firstLine.first;
		int version = firstLine.second;

		vector<string> sequences;
		for (auto& x : seq.NamesAndOffsets) {
			x.first = ParseSequenceName(x.first).first;
			sequences.push_back(x.first);
		}
		sequenceStrings.push_back(seq.Sequences);
		string& sequence = sequenceStrings.back();
		ToUpper(sequence);
		size_t excludedCount = ExcludeN(sequence, c);
		if (c.GetMaxGap() > sequence.length() && c.GetMaxGap() != numeric_limits<unsigned>::max()) throw invalid_argument("MaxGap can not be greater then sequence length.");

		size_t median = sequence.length() + 1;
		SearchType* searchType;
		size_t nonAlphabetLetters = 0;
		size_t r = 0;
		size_t t = 0;
		bool nonAlphabet[256] = {};

		if (c.GetProteins()) {
			alphabetSize = firstAlphabetSize;
			SetForProteins(sequence, c, alphabetSize, complement);
		}


		if (c.GetIsRepeat()) {
			if (c.GetIsMathematical()) {
				searchType = new MathematicalRepeat(sequence);
				t = 1;
				nonAlphabetLetters = CountNonAlphabetMakeComplement(sequence, complement, nonAlphabet).first;
				sequence = '\0' + sequence + '\0';
			}
			else {
				if (c.GetProteins()) {
					cout << "Complementary types can not be used with option proteins.\n " << endl;
					return;
				}
				auto a = CountNonAlphabetMakeComplement(sequence, complement, nonAlphabet);
				sequence = '\0' + sequence + '\0' + a.second + '\0';
				searchType = new BiologicalRepeat(median, sequence);
				nonAlphabetLetters = a.first;
				r = 1;
			}
		}
		else {
			if (c.GetIsMathematical()) {
				string s = sequence;
				nonAlphabetLetters = CountNonAlphabetMakeComplement(sequence, complement, nonAlphabet).first;
				reverse(s.begin(), s.end());
				sequence = '\0' + sequence + '\0' + s + '\0';
				searchType = new MathematicalPalindrome(median, sequence);
			}
			else {
				if (c.GetProteins()) {
					cout << "Complementary types can not be used with option proteins.\n " << endl;
					return;
				}
				string s = sequence;
				reverse(s.begin(), s.end());
				auto a = CountNonAlphabetMakeComplement(s, complement, nonAlphabet);
				nonAlphabetLetters = a.first;
				sequence = '\0' + sequence + '\0' + a.second + '\0';
				searchType = new BiologicalPalindrome(median, sequence);
			}
		};

		string letters;
		if (nonAlphabetLetters != 0)  letters = NonAlphabetLetters(nonAlphabet, excludedCount);
		if (repeatPrinter->InitializeForEverySequence(sequenceName, sequences, version, c, alphabetSize, nonAlphabetLetters, letters, excludedCount)) continue;

		vector<size_t> suftab(sequence.size());
		divsufsort(reinterpret_cast<const unsigned char*>(sequence.c_str()), reinterpret_cast<saidx_t*> (suftab.data()), sequence.size());
		auto subinv = MakeSubInv(suftab);
		auto lcp = MakeLCP(sequence, suftab, subinv);

		auto finder = [](size_t pos, const pair<string, size_t>& item) {
			return pos < item.second;
		};

	   if (!c.GetProbability()) {
			vector<size_t> repeatCount;

			auto forEachRepeatNoStat = [&](size_t x, size_t y, size_t len) {
				auto repeat = StringSlice(sequence, x, len);

				if (c.GetPrintInstances()) {
					auto seqX = upper_bound(seq.NamesAndOffsets.begin(), seq.NamesAndOffsets.end(), x, finder) - 1;
					auto seqY = upper_bound(seq.NamesAndOffsets.begin(), seq.NamesAndOffsets.end(), y, finder) - 1;

					repeatPrinter->OutputPairs(x - seqX->second, y - seqY->second, seqX - seq.NamesAndOffsets.begin(), seqY - seq.NamesAndOffsets.begin(), repeat, searchType->GetComplement(repeat));
				}

				if (repeatCount.size() <= len) {
					repeatCount.resize(len + 1, 0);
				}

				repeatCount[len]++;
			};

			ComputeLcpIntervals(lcp, suftab, sequence, c.GetFragmentLength(), searchType, forEachRepeatNoStat, c.GetMaxGap(), c.GetMotif());

			if (writeInfo) cout << '\n';
			repeatPrinter->AfterEverySequence(repeatCount);
		}
		else {
		   if (r != 1 && (alphabetSize > 22 || alphabetSize < 4)) {
			   cout << "Because of alphabet cardinality probability estimation is not performed." << endl;
			   return;
		   }

		   size_t seqLength;
		   if (c.GetIsMathematical() == true && c.GetIsRepeat() == true) seqLength = sequence.size() - 3;
		   else seqLength = (sequence.size() - 3) / 2 - 1;

		   FindSuggestLength(c, seqLength, alphabetSize, t, sequenceName);

		   shared_lifetime_unordered_map<StringSlice, size_t> strings(1000, hash<StringSlice>(), equal_to<StringSlice>(), shared_lifetime_allocator<pair<StringSlice, size_t>>(true));

		   auto forEachRepeat = [&](size_t x, size_t y, size_t len) {
			   auto repeat = StringSlice(sequence, x, len);
			   auto it = strings.find(repeat);
			   if (it == strings.end()) {
				   strings.insert(make_pair(repeat, 1));
			   }
			   else {
				   ++(it->second);
			   }
		   };
		   if (writeInfo && first) {
			   cout << "Working on repeats..." << endl;
			   first = false;
		   }

			ComputeLcpIntervals(lcp, suftab, sequence, c.GetFragmentLength(), searchType, forEachRepeat, c.GetMaxGap(), c.GetMotif());

			shared_lifetime_unordered_set<StringSlice> passedStatistics(1000, hash<StringSlice>(), equal_to<StringSlice>(), shared_lifetime_allocator<StringSlice>(true));

			if (c.GetIsMathematical()) {
				if (c.GetIsRepeat()) {
					ComputeStatistics<MathematicalRepeatEstimator>(median - 1 - excludedCount, strings, alphabetSize, c.GetPVal(), passedStatistics, multiplier, motifSize);
				}
				else {
					ComputeStatistics<MathematicalPalindromeEstimator>(median - 1 - excludedCount, strings, alphabetSize, c.GetPVal(), passedStatistics, multiplier, motifSize);
				}
			}
			else {
				if (c.GetIsRepeat()) {
					ComputeStatistics<BiologicalRepeatEstimator>(median - 1 - excludedCount, strings, alphabetSize, c.GetPVal(), passedStatistics, multiplier, motifSize);
				}
				else {
					ComputeStatistics<BiologicalPalindromeEstimator>(median - 1 - excludedCount, strings, alphabetSize, c.GetPVal(), passedStatistics, multiplier, motifSize);
				}
			}

			vector<size_t> repeatCount;
			auto forEachRepeat1 = [&](size_t x, size_t y, size_t len) {
				auto rep = StringSlice(sequence, x, len);

				if (passedStatistics.find(rep) != passedStatistics.end()) {
					if (c.GetPrintInstances()) {
						auto seqX = upper_bound(seq.NamesAndOffsets.begin(), seq.NamesAndOffsets.end(), x, finder) - 1;
						auto seqY = upper_bound(seq.NamesAndOffsets.begin(), seq.NamesAndOffsets.end(), y, finder) - 1;

						StringSlice repeat = StringSlice(sequence, x, len);

						repeatPrinter->OutputPairs(x - seqX->second, y - seqY->second, seqX - seq.NamesAndOffsets.begin(), seqY - seq.NamesAndOffsets.begin(), repeat, searchType->GetComplement(repeat));
					}

					if (repeatCount.size() <= len) {
						repeatCount.resize(len + 1, 0);
					}

					repeatCount[len]++;
				}
			};

			ComputeLcpIntervals(lcp, suftab, sequence, c.GetFragmentLength(), searchType, forEachRepeat1, c.GetMaxGap(), c.GetMotif());

			repeatPrinter->AfterEverySequence(repeatCount);

			strings.get_allocator().fast_free_on();
			passedStatistics.get_allocator().fast_free_on();
		}
	}
}

struct RepeatCount {
	int WordLen, Pairs, Instances;
};

vector<RepeatCount> ParseCsv(const string& file) {
	vector<RepeatCount> ret;

	ifstream matPalCsv(file);
	for (;;) {
		string line;
		getline(matPalCsv, line);
		if (line.empty()) break;

		istringstream seg(line);
		RepeatCount rc;
		char comma;
		seg >> rc.WordLen >> comma >> rc.Pairs >> comma >> rc.Instances;
		ret.push_back(rc);
	}

	return ret;
}

int main(int argc, char* argv[]) {
	try {
		ios_base::sync_with_stdio(false);
		setlocale(LC_ALL, "C");
		auto start = clock();
		MainImplementation(argc, argv);
		auto end = clock();
		cout << '\n';
		cout << "Elapsed time: " << (end - start) / (double)CLOCKS_PER_SEC << " seconds." << '\n';
		return 0;
	}
	catch (const invalid_argument& ex) {
		cout << endl;
		cerr << ex.what() << '\n';
		cout << endl;
		cout << Instructions() << '\n';
		cout << endl;
		return 1;

	}
}

