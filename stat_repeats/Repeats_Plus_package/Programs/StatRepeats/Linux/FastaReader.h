#ifndef FASTA_READER_H 
#define FASTA_READER_H 

#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <memory>
#include <tuple>

struct FastaSequence {
	std::vector<std::pair<std::string, size_t>> NamesAndOffsets;
	std::string Sequences;
};

class FastaFileIterator {
public:
	typedef FastaSequence ElementType;

private:
	std::unique_ptr<std::ifstream> inputFile_;
	bool readAll_;
	ElementType current_;
	bool hasCurrent_;

	std::tuple<bool, std::string, std::string> ReadNextSequence();
	bool ReadNext();

public:
	FastaFileIterator(const std::string& inputFileName, bool readAll);
	FastaFileIterator();

	FastaFileIterator& operator++();

	ElementType& operator * ();
	ElementType* operator -> ();

	friend bool operator == (const FastaFileIterator& x, const FastaFileIterator& y) {
		if (&x == &y) return true;
		return !x.hasCurrent_ && !y.hasCurrent_;
	}

	friend bool operator != (const FastaFileIterator& x, const FastaFileIterator& y) {
		return !(x == y);
	}
};

class FastaFileView {
	std::string fileName_;
	bool readAll_;

public:
	FastaFileView(const std::string& fileName, bool readAll) : fileName_ (fileName), readAll_(readAll) {
	}

	FastaFileIterator begin() {
		return FastaFileIterator(fileName_, readAll_);
	}

	FastaFileIterator end() {
		return FastaFileIterator();
	}
};

#endif

