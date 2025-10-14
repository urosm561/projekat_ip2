#include "CommandLine.h"
#include "CommandLineParser.h"
#include "StringSlice.h"
#include <limits>

using namespace std;


template<typename T, typename X>
vector<T>& operator << (vector<T>& vec, const X& t) {
	vec.push_back(t);
	return vec;
}

template<typename T, typename F>
T ReadCommandLineArgument(const unordered_map<string, string>& args, const string& name, T defaultVal, F reader) {
	auto it = args.find(name);
	if (it == args.end()) return defaultVal;

	try {
		return reader(it->second);
	}
	catch (const exception& ex) {
		auto msg = ex.what();
		if (msg != nullptr && *msg != '\0') {
			throw invalid_argument("Invalid value for argument " + name + ": " + msg);
		}
		else {
			throw invalid_argument("Invalid value for argument " + name + ".");
		}
	}
}

std::array<char, 256> CreateMapping() {
	std::array<char, 256> mapping;

	for (size_t i = 0; i < mapping.size(); i++) {
		mapping[i] = static_cast<char> (i);
	}

	return std::move(mapping);
}

void ProteinGroupMap(std::array<char, 256>& mapping, std::string letters, char mapTo, int& x) {
	for (auto ltr : letters) {
		auto a = static_cast<unsigned char>(ltr);
		if (mapping[a] != a) {
			throw invalid_argument(string("Letter ") + ltr + " is mapped more then once.");
		}
		else mapping[a] = mapTo;
	}

	x = -static_cast<int> (letters.length() - 1);
}

void Aliphatic(std::array<char, 256>& mapping, int &x) {
	ProteinGroupMap(mapping, "IVL", '0', x);
}

void Sulphur(std::array<char, 256>& mapping, int &x) {
	ProteinGroupMap(mapping, "MC", '1', x);
}

void Tiny(std::array<char, 256>& mapping, int &x) {
	ProteinGroupMap(mapping, "AGCS", '2', x);
}

void Aromatic(std::array<char, 256>& mapping, int &x) {
	ProteinGroupMap(mapping, "FYWH", '3', x);
}

void Hydrophobic(std::array<char, 256>& mapping, int &x) {
	ProteinGroupMap(mapping, "ACTKHWYFMILV", '4', x);
}

void Charged(std::array<char, 256>& mapping, int &x) {
	ProteinGroupMap(mapping, "DEHKR", '5', x);
}

void Positiv(std::array<char, 256>& mapping, int &x) {
	ProteinGroupMap(mapping, "HKR", '6', x);
}

void Polar(std::array<char, 256>& mapping, int &x) {
	ProteinGroupMap(mapping, "NQSTCDEHKRYW", '7', x);
}

void Acidic(std::array<char, 256>& mapping, int &x) {
	ProteinGroupMap(mapping, "NQ", '8', x);
}

void Small(std::array<char, 256>& mapping, int &x) {
	ProteinGroupMap(mapping, "VPAGCTSDN", '9', x);
}

void Hydroxylic(std::array<char, 256>& mapping, int &x) {
	ProteinGroupMap(mapping, "ST", '@', x);
}


CommandLine::CommandLine(int argsLength, char** args) {
	initializer_list<tuple<string, string, string>> shortcuts = {
		make_tuple("ic", "type", "ic"),
		make_tuple("dc", "type", "dc"),
		make_tuple("in", "type", "in"),
		make_tuple("dn", "type", "dn"),
		make_tuple("split", "split", "true"),
		make_tuple("proteins", "proteins", "true"),
		make_tuple("dna", "dna", "true"),
		make_tuple("rna", "rna", "true"),
		make_tuple("logging", "logging", "false"),
		make_tuple("noprobability", "noprobability", "true"),
		make_tuple("aliphatic", "aliphatic", "true"),
		make_tuple("sulphur", "sulphur", "true"),
		make_tuple("tiny", "tiny", "true"),
		make_tuple("aromatic", "aromatic", "true"),
		make_tuple("hydrophobic", "hydrophobic", "true"),
		make_tuple("charged", "charged", "true"),
		make_tuple("positive", "positive", "true"),
		make_tuple("polar", "polar", "true"),
		make_tuple("acidic", "acidic", "true"),
		make_tuple("small", "small", "true"),
		make_tuple("hydroxylic", "hydroxylic", "true"),
		make_tuple("msr", "msr", "gi|*|*|Combined FASTA Sequence|")
	};

	auto parsed = ParseCommandLine(
		argsLength, args,
		{ "input", "length" },
		{ "printinstances", "maxgap", "type", "embeddedsql", "pvalue", "exclude", "load", "commitcount", 
		"database", "db2output", "complement", "output", "aliphatic", "sulphur", "tiny","aromatic", "hydrophobic","computestatistics",
		"charged", "positive", "polar","acidic","small","hydroxylic", "msr" , "motif", "searchmotif"},
		shortcuts
		);

	if (parsed["help"] != "true") {
		inputFileName_ = parsed["input"];

		fragmentLength_ = atoi(parsed["length"].c_str());

		auto stringToUnsigned = [](const string& val) { return stoul(val); };
		auto stringToDouble = [](const string& val) { return stod(val); };
		auto stringToInt = [](const string& val) { return stoi(val); };
		auto stringToString = [](const string& val) { return val; };
		auto stringToBool = [](const string& val) {
			if (val == "true" || val == "1") return true;
			else if (val == "false" || val == "0") return false;
			else throw invalid_argument(""); };
		auto stringToType = [](const string& val) {
			if (CaseInsensitiveStartsWith("dn", val)) return 0;
			else if (CaseInsensitiveStartsWith("dc", val)) return 1;
			else if (CaseInsensitiveStartsWith("in", val)) return 2;
			else if (CaseInsensitiveStartsWith("ic", val))return 3;
			else throw invalid_argument("Type must be dn, dc, in or ic"); };

		maxGap_ = ReadCommandLineArgument(parsed, "maxgap", numeric_limits<unsigned>::max(), stringToUnsigned);

		outputFileName_ = ReadCommandLineArgument(parsed, "output", string(), stringToString);
		motif_ = ReadCommandLineArgument(parsed, "motif", string(), stringToString);
		motif1_ = ReadCommandLineArgument(parsed, "searchmotif", string(), stringToString);
		complementFileName_ = ReadCommandLineArgument(parsed, "complement", string(), stringToString);
		inputDb2FileName_ = ReadCommandLineArgument(parsed, "db2fileinput", string(), stringToString);
		database_ = ReadCommandLineArgument(parsed, "database", string(), stringToString);
		db2output_ = ReadCommandLineArgument(parsed, "db2output", string(), stringToString);
#ifdef NOODBC
		database_ = "";
		db2output_ = "";
#endif

		printInstances_ = ReadCommandLineArgument(parsed, "printinstances", true, stringToBool);
		probability_ = !(ReadCommandLineArgument(parsed, "noprobability", false, stringToBool));
		computeStatistics_ = ReadCommandLineArgument(parsed, "computestatistics", false, stringToBool);
		pval_ = ReadCommandLineArgument(parsed, "pvalue", 0.05, stringToDouble);
		if (pval_ > 1 || pval_ <= 0) throw invalid_argument("Pvalue must be a number between 0 and 1.");
		split_ = ReadCommandLineArgument(parsed, "split", false, stringToBool);
		proteins_ = ReadCommandLineArgument(parsed, "proteins", false, stringToBool);
		rna_ = ReadCommandLineArgument(parsed, "rna", false, stringToBool);
		dna_ = ReadCommandLineArgument(parsed, "dna", false, stringToBool);
		staticsql_ = ReadCommandLineArgument(parsed, "embeddedsql", string(), stringToString);
		exclude_ = ReadCommandLineArgument(parsed, "exclude", -1, stringToInt);
		load_ = ReadCommandLineArgument(parsed, "load", string(), stringToString);
		commitCount_ = ReadCommandLineArgument(parsed, "commitcount", 0, stringToInt);
		logging_ = ReadCommandLineArgument(parsed, "logging", true, stringToBool);
		auto a = ReadCommandLineArgument(parsed, "type", 0, stringToType);
		if (a == 0 || a == 2) isMathematical_ = true;
		else isMathematical_ = false;
		if (a == 0 || a == 1) isRepeat_ = true;
		else isRepeat_ = false;
		combinedSequenceName_ = ReadCommandLineArgument(parsed, "msr", string(), stringToString);

		mapping_ = CreateMapping();
		int x = 0;
		alphabetReduction_ = 0;
		newLetters_ = "";
		if (ReadCommandLineArgument(parsed, "aliphatic", false, stringToBool)) {
			Aliphatic(mapping_, x);
			alphabetReduction_ = x;
			newLetters_ = newLetters_ + '0';
		}
		if (ReadCommandLineArgument(parsed, "sulphur", false, stringToBool)) {
			Sulphur(mapping_, x);
			alphabetReduction_ = alphabetReduction_ + x;
			newLetters_ = newLetters_ + '1';
		}
		if (ReadCommandLineArgument(parsed, "tiny", false, stringToBool)) {
			Tiny(mapping_, x);
			alphabetReduction_ = alphabetReduction_ + x;
			newLetters_ = newLetters_ + '2';
		}
		if (ReadCommandLineArgument(parsed, "aromatic", false, stringToBool)) {
			Aromatic(mapping_, x);
			alphabetReduction_ = alphabetReduction_ + x;
			newLetters_ = newLetters_ + '3';
		}
		if (ReadCommandLineArgument(parsed, "hydrophobic", false, stringToBool)) {
			Hydrophobic(mapping_, x);
			alphabetReduction_ = alphabetReduction_ + x;
			newLetters_ = newLetters_ + '4';
		}
		if (ReadCommandLineArgument(parsed, "charged", false, stringToBool)) {
			Charged(mapping_, x);
			alphabetReduction_ = alphabetReduction_ + x;
			newLetters_ = newLetters_ + '5';
		}
		if (ReadCommandLineArgument(parsed, "positive", false, stringToBool)) {
			Positiv(mapping_, x);
			alphabetReduction_ = alphabetReduction_ + x;
			newLetters_ = newLetters_ + '6';
		}
		if (ReadCommandLineArgument(parsed, "polar", false, stringToBool)) {
			Polar(mapping_, x);
			alphabetReduction_ = alphabetReduction_ + x;
			newLetters_ = newLetters_ + '7';
		}
		if (ReadCommandLineArgument(parsed, "acidic", false, stringToBool)) {
			Acidic(mapping_, x);
			alphabetReduction_ = alphabetReduction_ + x;
			newLetters_ = newLetters_ + '8';
		}
		if (ReadCommandLineArgument(parsed, "small", false, stringToBool)) {
			Small(mapping_, x);
			alphabetReduction_ = alphabetReduction_ + x;
			newLetters_ = newLetters_ + '9';
		}
		if (ReadCommandLineArgument(parsed, "hydroxylic", false, stringToBool)) {
			Hydroxylic(mapping_, x);
			alphabetReduction_ = alphabetReduction_ + x;
			newLetters_ = newLetters_ + '@';
		}
	}
}
