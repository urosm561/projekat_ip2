#include "CommandLineParser.h"
#include <cstring>
#include <algorithm>
#include <initializer_list>

using namespace std;

bool CaseInsensitiveEqual(const string& str1, const string& str2, string::size_type len) {
	for (size_t x = 0; x < len; ++x) {
		if (toupper(str1[x]) != toupper(str2[x])) {
			return false;
		}
	}

	return true;
}

bool CaseInsensitiveStartsWith(const string& str1, const string& str2) { 
	if(str1.size() < str2.size()) { 
		return false; 
	} 

	return CaseInsensitiveEqual(str1, str2, str2.size());
}

bool CaseInsensitiveEqual(const string& str1, const string& str2) {
	if (str1.size() != str2.size()) return false;

	return CaseInsensitiveEqual(str1, str2, str1.size());
}

static string NormalizeSwitch(const string& x, const vector<string>& allSwitches){
	for(auto& y : allSwitches) {
		if(CaseInsensitiveEqual(x, y)) {
			return y;
		}
	}

	string longer = "";

	for(auto& y : allSwitches) {
		if(CaseInsensitiveStartsWith(y, x)) {
			if(longer != "") throw invalid_argument ("Command-line argument " + x + " matches both " + longer + " and " + y + ".");  			
			longer = y;
		}
	}

	if(longer != "") {
		return longer;
	}

	throw invalid_argument ("Command-line argument " + x + " not recognized.");  			
}

static void AddOption(unordered_map<string,string>& options, const string& option, const string& optionName, const string& optionValue) {
	if(options.find(optionName) != options.end()) 
		throw invalid_argument ("Command-line argument " + option + " has already been set: " + optionName + '=' + options[optionName] + ".");  			
	options.insert (make_pair(optionName,optionValue));
}

unordered_map<string, string> ParseCommandLine(int argsLength, char** args, const initializer_list<string>& mandatoryArguments, const initializer_list<string>& otherArguments, initializer_list<tuple<string, string, string>>& shortcuts) {
	vector<string> allSwitches(mandatoryArguments.begin(), mandatoryArguments.end());
	allSwitches.insert(allSwitches.end(),otherArguments.begin(), otherArguments.end());
	for(auto& x : shortcuts) allSwitches.push_back (get<0> (x));

	unordered_map<string,string> ret;
	vector <string> secondPass;


	for(int i = 1; i < argsLength; ++i) {
		string arg = args[i];
		if(arg.size() != 0 && (arg[0] == '-' || arg[0] == '/')) {
			string normal = NormalizeSwitch(arg.substr(1), allSwitches);
			string second, third;

			for (auto& x : shortcuts) {
				if(get<0>(x) == normal) {
					second = get<1> (x);
					third = get<2> (x);
				};
			}
			
			if(second != "") {
				AddOption(ret, normal, second, third);
			}
			else {
				if(i == argsLength - 1) throw invalid_argument("Command-line options end unexpectedly, value for option + " + normal + " is missing.");
				++i;
				AddOption(ret, normal, normal, args[i]);
			}
		}
		else {
			secondPass.push_back(arg);
		}
	}
	
	reverse(secondPass.begin(), secondPass.end());

	for (auto& x : mandatoryArguments) {
		if(ret.find(x) == ret.end()) {
			if(secondPass.size() == 0) throw invalid_argument("Argument input file and length must be specified.");
			ret.insert (make_pair(x, secondPass.back()));
			secondPass.pop_back();
		}
	}

	if(!secondPass.empty()) throw invalid_argument ("Unrecognized option: " + secondPass.back());

	return ret;
}
