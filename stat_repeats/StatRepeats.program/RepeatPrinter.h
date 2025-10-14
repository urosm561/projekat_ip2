#ifndef REPEATPRINTER_H
#define REPEATPRINTER_H

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <utility>
#include <fstream>
#include "StringSlice.h"
#include "CommandLine.h"
#include "WriteHeader.h"
#include "SharedLifetimeAllocator.h"

template<class T>
using shared_lifetime_unordered_set = std::unordered_set<T, std::hash<T>, std::equal_to<T>, shared_lifetime_allocator<T>>;

template<class K, class V>
using shared_lifetime_unordered_map = std::unordered_map<K, V, std::hash<K>, std::equal_to<K>, shared_lifetime_allocator<std::pair<K, V>>>;

namespace std
{
	template<typename X, typename Y>
	struct hash<pair<X, Y>>
	{
		size_t operator() (const pair<X, Y>& pair) const
		{
			auto h1 = hash<X>() (pair.first);
			auto h2 = hash<Y>() (pair.second);

			return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h2 >> 2));
		}
	};
}


struct RepeatPrinter
{
	virtual void OutputPairs(size_t x, size_t y, size_t xSeq, size_t ySeq, StringSlice repeat, StringSlice complement) = 0;
	virtual void AfterEverySequence(std::vector<size_t> repeatCount) = 0;
	virtual void Initialize(const CommandLine&, const size_t alphabetSize) = 0;
	virtual bool InitializeForEverySequence(const std::string& sequenceName, const std::vector<std::string>& containedSequences, const int version, const CommandLine& c, const size_t alphabetSize, const size_t nonAlphabetLetters, const std::string& letters, const size_t excludedCount) = 0;

	virtual ~RepeatPrinter() {}
};

class ConsoleOutput : public RepeatPrinter
{
	std::vector<std::string> containedSequences_;
	bool printSequenceNames_;
public:
	ConsoleOutput();

	void OutputPairs(size_t x, size_t y, size_t xSeq, size_t ySeq, StringSlice repeat, StringSlice complement) override;

	void AfterEverySequence(std::vector<size_t> repeatCount);

	virtual void Initialize(const CommandLine& c, const size_t alphabetSize);

	bool InitializeForEverySequence(const std::string& sequenceName, const std::vector<std::string>& containedSequences, const int version, const CommandLine& c, const size_t alphabetSize, const size_t nonAlphabetLetters, const std::string& letters, const size_t excludedCount) override;
};

template<typename K>
void IncrementCount(std::unordered_map<K, size_t>& map, const K& key);

class StatisticsOutput : public RepeatPrinter
{
	std::unordered_map<StringSlice, size_t> counts_;
	std::unordered_map<std::pair<size_t, size_t>, size_t> totals_;

public:
	virtual ~StatisticsOutput();

	void OutputPairs(size_t x, size_t y, size_t xSeq, size_t ySeq, StringSlice repeat, StringSlice complement) override;

	void AfterEverySequence(std::vector<size_t> repeatCount);

	virtual void Initialize(const CommandLine& c, const size_t alphabetSize);

	bool InitializeForEverySequence(const std::string& sequenceName, const std::vector<std::string>& containedSequences, const int version, const CommandLine& c, const size_t alphabetSize, const size_t nonAlphabetLetters, const std::string& letters, const size_t excludedCount) override;
};

class FileOutput : public RepeatPrinter
{
	std::string _sequenceName;
	int _version;
	FILE* _writer;
	std::vector<std::string> sequences_;
	bool printSequenceNames_;

public:
	FileOutput(std::string fileName);

	~FileOutput(); 

	void OutputPairs(size_t x, size_t y, size_t xSeq, size_t ySeq, StringSlice repeat, StringSlice complement) override;

	void AfterEverySequence(std::vector<size_t> repeatCount);

	void Initialize(const CommandLine& c, const size_t alphabetSize);

	bool InitializeForEverySequence(const std::string& sequenceName, const std::vector<std::string>& containedSequences, const int version, const CommandLine& c, const size_t alphabetSize, const size_t nonAlphabetLetters, const std::string& letters, const size_t excludedCount) override;
};

class LoadOutput : public RepeatPrinter
{
	std::string _sequenceName;
	int _version;
	FILE* _writer;
	FILE* _stat;
	FILE* _id;

	std::vector<std::string> sequences_;
	bool printSequenceNames_;

public:
	LoadOutput(std:: string name);

	~LoadOutput();

	void OutputPairs(size_t x, size_t y, size_t xSeq, size_t ySeq, StringSlice repeat, StringSlice complement) override;

	void AfterEverySequence(std::vector<size_t> repeatCount);

	void Initialize(const CommandLine& c, const size_t alphabetSize);

	bool InitializeForEverySequence(const std::string& sequenceName, const std::vector<std::string>& containedSequences, const int version, const CommandLine& c, const size_t alphabetSize, const size_t nonAlphabetLetters, const std::string& letters, const size_t excludedCount) override;
};

class SplitOutput : public RepeatPrinter
{
	int _version;
	FILE* _writer;
	FILE* _stat;
	FILE* _id;
	std::string _sequenceName;
	bool _firstSequence;

	std::vector<std::string> sequences_;
	bool printSequenceNames_;

	~SplitOutput(); 

public:
	SplitOutput();

	void OutputPairs(size_t x, size_t y, size_t xSeq, size_t ySeq, StringSlice repeat, StringSlice complement) override;

	void AfterEverySequence(std::vector<size_t> repeatCount);

	virtual void Initialize(const CommandLine& c, const size_t alphabetSize);

	bool InitializeForEverySequence(const std::string& sequenceName, const std::vector<std::string>& containedSequences, const int version, const CommandLine& c, const size_t alphabetSize, const size_t nonAlphabetLetters, const std::string& letters, const size_t excludedCount) override;
};


#ifdef ODBCOUTPUT

RepeatPrinter* CreateDataBaseOutput(const CommandLine& cl);
RepeatPrinter* CreateDb2Output(const CommandLine& cl);

#endif

#endif