#include "CommandLine.h"
#include "CommandLineParser.h"
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



CommandLine::CommandLine(int argsLength, char** args) {
	printInstances_ = true;
	exclude_ = -1;
	initializer_list<tuple<string, string, string>> shortcuts = {
		make_tuple("msr", "msr", "gi|*|*|Combined FASTA Sequence|"),
		make_tuple("proteins", "proteins", "true"),
		make_tuple("dna", "dna", "true"),
		make_tuple("rna", "rna", "true"),

	};

	auto parsed = ParseCommandLine(
		argsLength, args,
		{ "input", "motif"},
		{ "printinstances",  "exclude", "output", "msr" , "motif", "searchmotif"},
		shortcuts
		);

	if (parsed["help"] != "true") {
		inputFileName_ = parsed["input"];


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


		outputFileName_ = ReadCommandLineArgument(parsed, "output", string(), stringToString);
		motif_ = ReadCommandLineArgument(parsed, "motif", string(), stringToString);
		proteins_ = ReadCommandLineArgument(parsed, "proteins", false, stringToBool);
		rna_ = ReadCommandLineArgument(parsed, "rna", false, stringToBool);
		dna_ = ReadCommandLineArgument(parsed, "dna", false, stringToBool);

#ifdef NOODBC
		database_ = "";
		db2output_ = "";
#endif

		printInstances_ = ReadCommandLineArgument(parsed, "printinstances", true, stringToBool);
		exclude_ = ReadCommandLineArgument(parsed, "exclude", -1, stringToInt);
		combinedSequenceName_ = ReadCommandLineArgument(parsed, "msr", string(), stringToString);

	}
}
