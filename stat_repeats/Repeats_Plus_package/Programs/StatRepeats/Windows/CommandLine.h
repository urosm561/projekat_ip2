#ifndef COMMAND_LINE_H 
#define COMMAND_LINE_H 

#include <string>
#include <vector>
#include <tuple>
#include <stdexcept>
#include <array>   

class CommandLine {
	std::string inputFileName_;
	size_t fragmentLength_;
	bool printInstances_;
	bool probability_;
	bool isMathematical_;
	bool isRepeat_;
	std::string complementFileName_;
	std::string outputFileName_;
	std::string inputDb2FileName_;
	bool split_;
	bool logging_;
	bool proteins_;
	bool rna_;
	bool dna_;
	int commitCount_;
	std::string load_;
	unsigned maxGap_;
	double pval_;
	std::string database_;
	std::string db2output_;
	std::string staticsql_;
	bool computeStatistics_;
	int exclude_;
	std::string excludeLetter_;
	std::string combinedSequenceName_;
	std::array<char, 256> mapping_;
	int alphabetReduction_;
	std::string newLetters_;
	std::string motif_;
	std::string motif1_;
public:
	CommandLine(int argsLength, char** args);

	void SetFragmentLength(const size_t fragmentLength) {
		fragmentLength_ = fragmentLength;
	}

	size_t GetFragmentLength() const {
		return fragmentLength_;
	}

	bool GetIsMathematical() const {
		return isMathematical_;
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

	int GetExclude() const {
		return exclude_;
	}


	bool GetIsRepeat() const {
		return isRepeat_;
	}

	std::string GetInputFileName() const {
		return inputFileName_;
	}

	std::string GetComplementFileName() const {
		return complementFileName_;
	}

	bool GetPrintInstances() const {
		return printInstances_;
	}

	bool GetSplit() const {
		return split_;
	}

	std::string GetStaticSql() const {
		return staticsql_;
	}

	std::string GetLoad() const {
		return load_;
	}

	bool GetLogging() const {
		return logging_;
	}

	int GetCommitCount() const {
		return commitCount_;
	}

	bool GetProbability() const {
		return probability_;
	}

	double GetPVal() const {
		return pval_;
	}

	std::string GetOutputFileName() const {
		return outputFileName_;
	}

	std::string GetDatabase() const {
		return database_;
	}

	std::string GetDb2Output() const {
		return db2output_;
	}

	bool GetComputeStatistics() const {
		return computeStatistics_;
	}

	unsigned GetMaxGap() const {
		return maxGap_;
	}

	const std::string& CombinedSequenceName() const {
		return combinedSequenceName_;
	}

	int GetAlphabetReduction() const {
		return alphabetReduction_;
	}

	const std::array<char, 256>& GetMapping() const {
		return mapping_;
	}

	std::string GetNewLetters() const {
		return newLetters_;
	}

	std::string GetMotif() const {
		return motif_;
	}

	std::string GetMotif1() const {
		return motif1_;
	}

	std::string GetExcludeLetter() const {
		return excludeLetter_;
	}
};

#endif