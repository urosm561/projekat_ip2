#ifndef COMMAND_LINE_H 
#define COMMAND_LINE_H 

#include <string>
#include <vector>
#include <tuple>
#include <stdexcept>
#include <array>   

class CommandLine {
	std::string inputFileName_;
	bool printInstances_;
	std::string outputFileName_;
	int exclude_;
	std::string combinedSequenceName_;
	std::string motif_;
	bool proteins_;
	bool rna_;
	bool dna_;

public:
	CommandLine(int argsLength, char** args);


	int GetExclude() const {
		return exclude_;
	}



	std::string GetInputFileName() const {
		return inputFileName_;
	}

	bool GetPrintInstances() const {
		return printInstances_;
	}


	std::string GetOutputFileName() const {
		return outputFileName_;
	}

	const std::string& CombinedSequenceName() const {
		return combinedSequenceName_;
	}


	std::string GetMotif() const {
		return motif_;
	}
	bool GetProteins() const {
		return proteins_;
	}

	bool GetDna() const {
		return dna_;
	}

	bool GetRna() const {
		return rna_;
	}



};

#endif