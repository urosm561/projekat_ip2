#include <iostream>
#include <bitset>
#include <iomanip>
#include <locale>
#include <utility>
#include <ctime>
#include <algorithm>
#include <ctype.h>
#include <regex>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <stdexcept>
#include <array>
#include "FastaReader.h"
#include "CommandLineParser.h"
#include "CommandLine.h"
#include "RepeatPrinter.h"

using namespace std;

array<char, 256> MakeCompRules(const CommandLine& c) {
	array<char, 256> complement;
	for (size_t x = 0; x < complement.size(); ++x) complement[x] = static_cast<char>(x);

	string fileName = c.GetComplementFileName();
	if (fileName != "") {
		ifstream file(fileName);
		if (!file) {
			throw invalid_argument("Can't open file for complement rules.");
		}

		for (;;) {
			char letter = '\0';
			do
			{
				file.get(letter);
			} while (isspace(letter) && !file.eof());
			if (file.eof()) break;

			char compLetter = '\0';
			do
			{
				file.get(compLetter);
			} while (isspace(compLetter) && !file.eof());
			if (file.eof()) throw invalid_argument("There must be an even number of letters in the complement file.");;

			letter = toupper(letter);
			compLetter = toupper(compLetter);

			complement[static_cast<unsigned char> (letter)] = compLetter;
		}
	}
	else if (c.GetRna()) {
		complement['A'] = 'U';
		complement['U'] = 'A';
		complement['C'] = 'G';
		complement['G'] = 'C';
		complement['R'] = 'Y';
		complement['Y'] = 'R';
		complement['S'] = 'W';
		complement['W'] = 'S';
		complement['K'] = 'M';
		complement['M'] = 'K';
		complement['B'] = 'V';
		complement['V'] = 'B';
		complement['D'] = 'H';
		complement['H'] = 'D';
	}
	else {
		complement['A'] = 'T';
		complement['T'] = 'A';
		complement['C'] = 'G';
		complement['G'] = 'C';
		complement['R'] = 'Y';
		complement['Y'] = 'R';
		complement['S'] = 'W';
		complement['W'] = 'S';
		complement['K'] = 'M';
		complement['M'] = 'K';
		complement['B'] = 'V';
		complement['V'] = 'B';
		complement['D'] = 'H';
		complement['H'] = 'D';
	}

	return complement;
}

static const char* Instructions() {
	return
		"Usage:\n"
		"RepeatsPlus <input file name> <fragment length> [<options>]\n"
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
		"-maxgap number, where this is maximal gap between two results, default is finding all results, can be abbreviated as -max[gap] number. \n"
		"-motif string, where string is motif mask. Only repeats that match motif mask are found.\n"
		" \n"
		" \n"
		"Output options:\n"
		"Default for output is writing results on console.\n"
		"-output outputFileName, can be abbreviated as -o[utput] outputFileName.\n"
		"-split, option for making three output files. Sequence name is written in file named sequencename.id, statistics are written in file nemed sequencename.stat, and results are written in file named sequencename.load. Can be abbreviated as -sp[lit].\n"
		"-load name,option for making three output files. Sequence name is written in file named name.id, statistics are written in file nemed name.stat, and results are written in file named name.load,can be abbreviated as -loa[d] name. \n"
		"Only one output option can be chosen. Output options are outputfile,split,load."
		" \n"
		" \n"
		"-printinstances true/false, false for not printing instances, default is true. Can be abbreviated as -pri[ntinstances] true/false. \n"
		" \n"
		"\n"
		"\n";
}

bool EqualNucleotide(char c1, char c2) {
	static uint32_t table[26] = {
		0b00000000010001000011000100011000,
		0b00000000100010000001000001000000,
		0b00000001000001000011000010010010,
		0b00000010000010000001000001000000,
		0b00000000000000000001000000000000,
		0b00000000000000000001000000000000,
		0b00000001010000001001000110010000,
		0b00000010100000000001000001000000,
		0b00000000000000000001000000000000,
		0b00000000000000000001000000000000,
		0b00000000000010000001000001000000,
		0b00000000000000000001000000000000,
		0b00000010100000000001000000000000,
		0b11111111111111111111111111111111,
		0b00000000000000000001000000000000,
		0b00000000000000000001000000000000,
		0b00000000000000000001000000000000,
		0b00000010000010000001000000000000,
		0b00000000100010000001000000000000,
		0b00000001010001001001000000001010,
		0b00000000000000000001000000000000,
		0b00000010100010000001000000000000,
		0b00000010000000000001000001000000,
		0b00000000000000000001000000000000,
		0b00000000100000000001000001000000,
		0b00000000000000000001000000000000
	};

	if (c1 == c2) return true;

	if (c1 >= 'A' && c1 <= 'Z' && c2 >= 'A' && c2 <= 'Z') {
		size_t i = c1 - 'A';
		size_t j = 25 - (c2 - 'A');

		// vidi da li je j ti bit u help jednak 1? Ako jeste vrati true i mismatc = true;
		auto letter = table[i];
		if (letter & (1 << j)) {
			return true;
		};
	}

	return false;
}


bool EqualProtein(char c1, char c2) {
	static uint32_t table[26] = {
		0b00000000000000000000000000000100,
		0b00000000010000000001000000000100,
		0b00000000000000000000000000000100,
		0b00000001000000000000000000000100,
		0b00000000000000000000000000000101,
		0b00000000000000000000000000000100,
		0b00000000000000000000000000000100,
		0b00000000000000000000000000000100,
		0b00000000000000010000000000000100,
		0b00000000000000100100000000000100,
		0b00000000000000000000000000000100,
		0b00000000000000010000000000000100,
		0b00000000000000000000000000000100,
		0b00000001000000000000000000000100,
		0b00000000000000000000000000000100,
		0b00000000000000000000000000000100,
		0b00000000000000000000000000000101,
		0b00000000000000000000000000000100,
		0b00000000000000000000000000000100,
		0b00000000000000000000000000000100,
		0b00000000000000000000000000000100,
		0b00000000000000000000000000000100,
		0b00000000000000000000000000000100,
		0b11111111111111111111111111111111,
		0b00000000000000000000000000000100,
		0b00000000001000000000001000000100
	};
	if (c1 == c2) return true;

	if (c1 >= 'A' && c1 <= 'Z' && c2 >= 'A' && c2 <= 'Z') {
		size_t i = c1 - 'A';
		size_t j = 25 - (c2 - 'A');

		// vidi da li je j ti bit u help jednak 1? Ako jeste vrati true i mismatc = true;
		auto letter = table[i];
		if (letter & (1 << j)) {
			return true;
		};
	}
	return false;
}

bool EqualProtein(char c1, char c2, bool& mismatch) {
	if (c1 == c2) return true;
	
	if (c1 == 'B' && (c2 == 'D' || c2 == 'N')) {
		mismatch = true;
		return true;
	}
	if (c2 == 'B' && (c1 == 'D' || c1 == 'N')) {
		mismatch = true;
		return true;
	}
	if (c1 == 'Z' && (c2 == 'E' || c2 == 'Q')) {
		mismatch = true;
		return true;
	}
	if (c2 == 'Z' && (c1 == 'E' || c1 == 'Q')) {
		mismatch = true;
		return true;
	}
	if (c1 == 'J' && (c2 == 'I' || c2 == 'L')) {
		mismatch = true;
		return true;
	}
	if (c2 == 'J' && (c1 == 'I' || c1 == 'L')) {
		mismatch = true;
		return true;
	}
	if (c1 == 'X' || c2 == 'X') {
		mismatch = true;
		return true;
	}

	return false;
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

void ToUpper(string& s)
{
	for (auto it = s.begin(); it != s.end(); ++it) {
		*it = toupper(*it);
	}
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
	else {
		repeatPrinter = new ConsoleOutput();
	}

	if (outputCounter > 1) {
		throw invalid_argument("There can only be one output option. Output options are outputfile,split,load,database,db2output,staticsql.");
	}
	return repeatPrinter;
}

template<typename T>
class circular_buffer {
	vector<T> vec_;
	size_t capacity_, pos_;

public:
	circular_buffer(size_t capacity) {
		pos_ = 0;
		capacity_ = capacity;
		vec_.reserve(capacity);
	}

	void push_back(const T& value) {
		if (vec_.size() == capacity_) {
			vec_[pos_++] = value;
			pos_ = pos_ == capacity_ ? 0 : pos_;
		}
		else {
			vec_.push_back(value);
			pos_ = pos_ + 1 == capacity_ ? 0 : pos_ + 1;
		}
	}

	size_t size() const {
		return vec_.size();
	}

	T& front() {
		return vec_.size() < capacity_ ? vec_[0] : vec_[pos_];
	}
};

template<typename It1, typename It2, typename BinaryPredicate, typename OutputFunc>
void scan_repeats(It1 begin1, It1 end1, It2 begin2, size_t minLength, size_t numMismatches, BinaryPredicate equal, const OutputFunc& outputFunc) {
	circular_buffer<It1> startOfRuns(numMismatches + 1);
	startOfRuns.push_back(begin1);

	while (begin1 != end1) {
		if (!equal(*begin1, *begin2)) {
			It1 startOfRun = startOfRuns.front();

			auto len = static_cast<size_t> (begin1 - startOfRun);
			if (len >= minLength) {
				outputFunc(startOfRun, begin2 - len, len);
			}

			startOfRuns.push_back(begin1 + 1);
		}

		++begin1;
		++begin2;
	}

	{
		auto startOfRun = startOfRuns.front();
		auto len = static_cast<size_t> (begin1 - startOfRun);
		if (len >= minLength) {
			outputFunc(startOfRun, begin2 - len, len);
		}
	}
}

template<typename It1, typename It2, typename BinaryPredicate, typename OutputFunc>
void slide_one(It1 begin1, It1 end1, It2 begin2, size_t minLength, size_t numMismatches, BinaryPredicate equal, const OutputFunc& outputFunc) {
	for (auto start1 = begin1; start1 != end1; ++start1) {
		scan_repeats(start1, end1, begin2, minLength, numMismatches, equal, outputFunc);
	}
}

template<typename It1, typename It2, typename BinaryPredicate, typename OutputFunc>
void slide_both(It1 begin1, It1 end1, It2 begin2, size_t minLength, size_t numMismatches, BinaryPredicate equal, const OutputFunc& outputFunc) {
	for (auto start1 = begin1; start1 != end1; ++start1) {
		scan_repeats(start1, end1, begin2, minLength, numMismatches, equal, outputFunc);
	}

	auto end2 = begin2 + (end1 - begin1);

	for (auto start2 = begin2 + 1; start2 != end2; ++start2) {
		auto new_end1 = begin1 + (end2 - start2);
		scan_repeats(begin1, new_end1, start2, minLength, numMismatches, equal, outputFunc);
	}
}

string map_letters(const string& input, const array<char, 256>& complements) {
	string ret;
	ret.reserve(input.size());

	for (auto x : input) {
		ret.push_back(complements[x]);
	}

	return ret;
}

bool regex_search(const string_view& sv, const regex& r) {
	return regex_search(sv.begin(), sv.end(), r);
}

template<typename It1, typename It2>
bool regex_search_two(It1 begin1, It2 begin2, size_t len, const regex& regex) {
	It1 end1 = begin1 + len;

	while (begin1 != end1) {
		match_results<It1> match;
		bool found = std::regex_search(begin1, end1, match, regex);
		if (!found) break;
		auto pos = match.position(0);
		auto len = match.length(0);

		auto from2 = begin2 + pos;
		auto to2 = from2 + len;

		if (regex_match(from2, to2, regex) && equal(begin1, begin1 + pos, begin2) && equal(begin1 + pos + len, end1, to2)) {
			return true;
		}

		begin1 += pos + 1;
		begin2 += pos + 1;
	}

	return false;
}

template<typename BinaryPredicate, typename OutputFunc>
void find_repeats(const string& sequence, bool reverseDirection, const optional<array<char, 256>>& complements, size_t minLength, const optional<regex>& maybeMotif, size_t numMismatches, BinaryPredicate equal, const OutputFunc& outputFunc) {
	if (reverseDirection && complements) {
		auto reversedBothWays = map_letters(sequence, *complements);
		reverse(reversedBothWays.begin(), reversedBothWays.end());

		slide_both(sequence.cbegin(), sequence.cend(), reversedBothWays.cbegin(), minLength, numMismatches, equal, [&sequence, &reversedBothWays, &maybeMotif, &outputFunc](auto pos1, auto pos2, auto len) {
			auto altPos2 = sequence.cend() - (pos2 - reversedBothWays.cbegin()) - len;

			if (altPos2 >= pos1) {
				auto sv1 = string_view(&*pos1, len);
				auto sv2 = string_view(&*altPos2, len);

				auto altPos1 = reversedBothWays.cend() - (pos1 - sequence.cbegin()) - len;

				if (!maybeMotif || regex_search_two(pos1, pos2, len, *maybeMotif) || regex_search_two(altPos1, altPos2, len, *maybeMotif)) {
					outputFunc(pos1 - sequence.cbegin(), altPos2 - sequence.cbegin(), len);
				}
			}
		});
	}
	else if (reverseDirection) {
		auto reversedSequence = sequence;
		reverse(reversedSequence.begin(), reversedSequence.end());

		slide_both(sequence.cbegin(), sequence.cend(), reversedSequence.cbegin(), minLength, numMismatches, equal, [&sequence, &reversedSequence, &maybeMotif, &outputFunc](auto pos1, auto pos2, auto len) {
			auto seqPos2 = sequence.cend() - (pos2 - reversedSequence.cbegin()) - len;

			if (seqPos2 >= pos1) {
				string_view sv1(&*pos1, len);
				string_view sv2(&*seqPos2, len);

				if (!maybeMotif || regex_search_two(sv1.begin(), sv2.rbegin(), sv1.size(), *maybeMotif) || regex_search_two(sv1.rbegin(), sv2.begin(), sv1.size(), *maybeMotif)) {
					outputFunc(pos1 - sequence.cbegin(), seqPos2 - sequence.cbegin(), len);
				}
			}
		});
	}
	else if (complements) {
		auto reversedLetters = map_letters(sequence, *complements);

		slide_both(sequence.cbegin(), sequence.cend(), reversedLetters.cbegin(), minLength, numMismatches, equal, [&sequence, &reversedLetters, &maybeMotif, &outputFunc](auto pos1, auto pos2, auto len) {
			auto altPos2 = sequence.cbegin() + (pos2 - reversedLetters.cbegin());

			if (altPos2 >= pos1) {
				string_view sv1(&*pos1, len);
				string_view sv2(&*altPos2, len);

				auto altPos1 = reversedLetters.cbegin() + (pos1 - sequence.cbegin());

				if (!maybeMotif || regex_search_two(pos1, pos2, len, *maybeMotif) || regex_search_two(altPos1, altPos2, len, *maybeMotif)) {
					outputFunc(pos1 - sequence.cbegin(), altPos2 - sequence.cbegin(), len);
				}
			}
		});
	}
	else {
		slide_one(sequence.cbegin() + 1, sequence.cend(), sequence.cbegin(), minLength, numMismatches, equal, [&sequence, &maybeMotif, &outputFunc](auto pos1, auto pos2, auto len) {
			string_view sv1(&*pos1, len);
			string_view sv2(&*pos2, len);

			if (!maybeMotif || regex_search_two(sv1.begin(), sv2.begin(), sv1.size(), *maybeMotif)) {
				outputFunc(pos1 - sequence.cbegin(), pos2 - sequence.cbegin(), len);
			}
		});
	}
}

int main(int argc, char** argv)
{
	try {
		cout << "RepeatPlus v1.1 Copyright ( C ) 2018 Ana Jelovic. All Rights Reserved." << endl;
		auto start = clock();

		CommandLine cmdLine = CommandLine(argc, argv);

		size_t minLength = cmdLine.GetFragmentLength();

		string motifString = cmdLine.GetMotif();
		ToUpper(motifString);
		
		optional<regex> maybeMotif;
		size_t numMismatches = 0;
		if (cmdLine.GetMotif() != "") {
			maybeMotif = regex{ motifString, regex_constants::ECMAScript | regex_constants::optimize };
			numMismatches = count_if(motifString.begin(), motifString.end(), [](char c) { return c == '.' || c == '['; });
		}

		bool reverseDirection = !cmdLine.GetIsRepeat();
		optional<array<char, 256>> complements;
		if (!cmdLine.GetIsMathematical()) {
			if (cmdLine.GetProteins()) throw invalid_argument("Complementary repeats are not provided for protein sequences.");
			complements = MakeCompRules(cmdLine);
		}

		size_t outputCounter = 0;
		bool writeInfo = true;
		auto repeatPrinter = unique_ptr<RepeatPrinter>(SetOutput(cmdLine, writeInfo, outputCounter));

		repeatPrinter->Initialize(cmdLine);
		bool first = true;
		cout << "Reading file..." << endl;
		vector<string> sequenceStrings;
		for (auto& seq : FastaFileView(cmdLine.GetInputFileName(), !cmdLine.CombinedSequenceName().empty())) {
			auto firstLine = ParseSequenceName(cmdLine.CombinedSequenceName().empty() ? static_cast<const string&> (seq.NamesAndOffsets[0].first) : cmdLine.CombinedSequenceName());
			string sequenceName = firstLine.first;
			int version = firstLine.second;
			vector<string> sequences;
			for (auto& x : seq.NamesAndOffsets) {
				x.first = ParseSequenceName(x.first).first;
				sequences.push_back(x.first);
			}
			sequenceStrings.push_back(seq.Sequences);
			string& sequence1 = sequenceStrings.back();
			ToUpper(sequence1);

			vector<size_t> repeatCount;

			size_t numRepeats = 0;
			if (first) {
				cout << "Finding repeats..." << endl;
				first = false;
			}
			unordered_map<string, size_t> strings;

			auto forEachRepeat = [&](size_t x, size_t y, size_t len) {
				auto repeat = StringSlice(sequence1, x, len);
				auto complement = StringSlice(sequence1, y, len);

				repeatPrinter->OutputPairs(sequenceName, x + 1, y + 1, x + 1 + len, y + 1 + len, repeat, complement);

				if (repeatCount.size() <= len) {
					repeatCount.resize(len + 1, 0);
				}

				repeatCount[len]++;
			};

			if (repeatPrinter->InitializeForEverySequence(sequenceName, sequences, version, cmdLine)) {
				continue;
			}

			if (maybeMotif) {
				find_repeats(sequence1, reverseDirection, complements, minLength, maybeMotif, numMismatches, std::equal_to(), forEachRepeat);
			}
			else if (cmdLine.GetProteins()) {
				find_repeats(sequence1, reverseDirection, complements, minLength, maybeMotif, numMismatches, [](char c1, char c2) { return EqualProtein(c1, c2); }, forEachRepeat);
			}
			else {
				find_repeats(sequence1, reverseDirection, complements, minLength, maybeMotif, numMismatches, [](char c1, char c2) { return EqualNucleotide(c1, c2); }, forEachRepeat);
			}

			repeatPrinter->AfterEverySequence(repeatCount);
		}

		auto end = clock();
		cout << '\n';
		cout << "Elapsed time: " << (end - start) / (double)CLOCKS_PER_SEC << " seconds." << '\n';

		return 0;
	}
	catch (const regex_error e) {
		cout << endl;
		cout << e.what() << '\n';
		cout << endl;
		cout << Instructions() << '\n';
		cout << endl;
		return 1;
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

