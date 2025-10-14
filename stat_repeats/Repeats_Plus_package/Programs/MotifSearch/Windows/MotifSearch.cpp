
#include <iostream>
#include <fstream>
#include <bitset>
#include <iomanip>
#include "FastaReader.h"
#include "CommandLineParser.h"
#include "CommandLine.h"
#include <locale>
#include <utility>
#include <algorithm>
#include <ctype.h>
#include <regex>
#include <limits>
#include <cmath>
#include<stdio.h>

#pragma warning (disable:4996)


using namespace std;

void ToUpper(string& s)
{
	for (auto it = s.begin(); it != s.end(); ++it) {
		*it = toupper(*it);
	}
}

string MotifToRegex(const string& motif, vector<vector<char>>& v, const vector<char>& alphabet, const size_t alphabetSize, size_t& motifSize) {
	size_t i = 0;
	string regex = "";
	auto m = motif;
	ToUpper(m);
	while (i < motif.size()) {
		if (motif[i] == '.') {
			regex += '.';
			motifSize++;
			v.push_back(alphabet);
		}
		else if (isalpha(motif[i])) {
			regex += m[i];
			motifSize++;
			v.push_back({ motif[i] });
		}
		else if (motif[i] == '[' && motif[i + 1] != '^') {
			vector<char> w = {};
			regex += m[i];
			size_t j = 1;
			size_t k = 0;
			while (i + j < motif.size() && motif[i + j] != ']') {
				if (!isalpha(motif[i + j])) throw invalid_argument("Motif mask is not correct.");
				regex += m[i + j];
				w.push_back(motif[i + j]);
				j++;
			}
			if (j == 1) throw invalid_argument("Motif mask is not correct.");
			i = i + j;
			regex += m[i];
			v.push_back(w);
			motifSize++;
		}
		else if (motif[i] == '[' && motif[i + 1] == '^') {
			v.push_back(alphabet);
			regex += "[^";
			size_t j = 2;
			size_t k = 0;
			while (i + j < motif.size() && motif[i + j] != ']') {
				if (!isalpha(motif[i + j])) throw invalid_argument("Motif mask is not correct.");
				regex += m[i + j];
				auto position = find(v[i].begin(), v[i].end(), m[i + j]);
				if (position != v[i].end()) v[i].erase(position);
				j++;
			}
			if (j == 1) throw invalid_argument("Motif mask is not correct.");
			i = i + j;
			regex += "]";
			motifSize++;
		}
		else throw invalid_argument("Motif mask is not correct.");
		i++;
	}
	return regex;
}


static const char* Instructions() {
	return 		"Usage:\n"
		"MotifSearch <input file name> -motif motifString  [<options>]\n"
		"\n"
		" \n"

		"Mandatory options are:\n"
		"inputFileName specifies the input file name, does not have a default value, must be specified.\n"
		" \n"
		" \n"
		"Output options:\n"
		"Default for output is writing results on console.\n"
		"-output outputFileName, can be abbreviated as -o[utput] outputFileName.\n"
		"-printinstances true/false, false for not printing instances, default is true. Can be abbreviated as -pri[ntinstances] true/false. \n"
		" \n"
		"\n"
		"\n";
}

size_t WordOverlap(const vector<char>& word) {

	for (size_t i = 1; i <= word.size() - 1; i++) {
		bool same = true;
		for (size_t j = 0; j < word.size() - i; j++) {
			if (word[j] != word[j + i]) {
				same = false;
				break;
			}
		}
		if (same == true) return word.size() - i;
	}
	return 0;

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

void MotifCoefficients(const vector<vector<char>>& motif, vector<char>& vec, unordered_map<size_t, size_t>& motifOverlapping) {
	if (motif.size() == vec.size()) {
		auto x = WordOverlap(vec);
		if (motifOverlapping.find(x) != motifOverlapping.end()) motifOverlapping.find(x)->second++;
		else motifOverlapping[x] = 1;
	}
	else {
		auto& characters = motif[vec.size()];
		for (const auto a : characters) {
			vec.push_back(a);
			MotifCoefficients(motif, vec, motifOverlapping);
			vec.pop_back();
		}
	}
}
void MakeCompRules(const CommandLine& c, size_t& alphabetSize, vector<char>& alphabet) {
	if (c.GetProteins()) {
		alphabetSize = 20;
		alphabet = { 'A','C','T','G','D','E','F','H','I','K','L','M','N','P','Q','R','S','V','W','Y' };
	}
	else if (c.GetRna()) {
		alphabetSize = 4;
		alphabet = { 'A', 'C', 'U', 'G' };

	}
	else {
		alphabetSize = 4;
		alphabet = { 'A', 'C', 'T', 'G' };
	}
}


int main(int argc, char** argv)
{
    try {
		CommandLine c = CommandLine(argc, argv);
		FILE* writer;
		writer = fopen(c.GetOutputFileName().c_str(), "w");
			
	

		string motifString = c.GetMotif();
		string reg;
		vector<char> alphabet;
		size_t alphabetSize = 0;

		MakeCompRules(c, alphabetSize, alphabet);

		size_t motifSize = 0;

		vector<vector<char>> v;
		if (c.GetMotif() != "") reg = MotifToRegex(c.GetMotif(), v, alphabet, alphabetSize, motifSize);
		unordered_map <size_t, size_t> multiplier;
		vector<char> w;
		if (c.GetMotif().size() != 0) MotifCoefficients(v, w, multiplier);

		size_t outputCounter = 0;
		bool writeInfo = true;
		bool first = true;
		//cout << "Reading file..." << endl;
		//auto start = clock();
		vector<string> sequenceStrings;
		for (auto& seq : FastaFileView(c.GetInputFileName(), !c.CombinedSequenceName().empty())) {
			auto firstLine = ParseSequenceName(c.CombinedSequenceName().empty() ? static_cast<const string&> (seq.NamesAndOffsets[0].first) : c.CombinedSequenceName());
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

			regex r{ reg };

			double lambda = (sequence1.size()-motifSize+1)/pow(alphabetSize, motifSize);
			double expected = 0;

			for (auto& a : multiplier) {
				double newLambda = 0;
				if (a.first != 0) {
					newLambda = lambda * (1 - 1/pow(alphabetSize, motifSize - a.first));
				}
				else newLambda = lambda;
				expected += newLambda*(a.second);
			}
			if (c.GetOutputFileName() != "") {

				fprintf(writer, "%s%s%s\n", "Working with  ", sequenceName.c_str(), " file.");
				fprintf(writer, "%s%s%s\n\n", "Motif mask is ", c.GetMotif().c_str(), ".");
			}
			else {
				cout << "Working with " << sequenceName << " file." << endl;
				cout << "Motif mask is " << c.GetMotif() << "." << endl;
				cout << endl;
			}

			int count = 0;
			for (sregex_iterator it = sregex_iterator(sequence1.begin(), sequence1.end(), r);  it != sregex_iterator(); it++) {
				count++;
				smatch match;
				match = *it;
				if (c.GetOutputFileName() != "") {

					fprintf(writer, "%s%s%zu%s%s\n", sequenceName.c_str(), ", ", match.position(0),", ", match.str(0).c_str());
				}
				else {
					cout << sequenceName << ", " << match.position(0) << ", " << match.str(0) << endl;
				}
			}
			if (c.GetOutputFileName() != "") {

				fprintf(writer, "\n%s%d%s\n", "Found ", count, " matches.");
				fprintf(writer, "%s%f%s\n","Expected ", expected, " in random sequence with the same length.");
			}
			else {
				cout << endl;
				cout << "\nFound " << count << " matches." << endl;
				cout << "Expected " << expected << " matches in a random sequence with the same length." << endl;
			}
		
		}

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
	catch (const invalid_argument & ex) {
		cout << endl;
		cerr << ex.what() << '\n';
		cout << endl;
		cout << Instructions() << '\n';
		cout << endl;
		return 1;
	}

}

