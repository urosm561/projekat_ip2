#include "FastaReader.h"
#include <stdexcept>
#include <utility>

using namespace std;

istream& Readline( istream& istr, string& s ) {
	ostringstream readbuffer;
	while(istr) {
		int c = istr.get();
		if(!istr)
			break;
		else if( c=='\n' )
			break;
		else if( c=='\r' ){
			if( istr.peek() == '\n' )
				istr.get();
			break;
		}
		else
			readbuffer.put(c);
	}
	s = readbuffer.str();
	if(s.length() > 0 && !istr)
		istr.clear();
	return istr;
}

FastaFileIterator::FastaFileIterator(const std::string& inputFileName, bool readAll) : inputFile_(new ifstream (inputFileName)), readAll_(readAll) {
	if (!*inputFile_) {
		throw invalid_argument("Input file does not exist.");
	}

	hasCurrent_ = ReadNext();
}

FastaFileIterator::FastaFileIterator() {
	hasCurrent_ = false;
}

FastaFileIterator& FastaFileIterator::operator++() {
	hasCurrent_ = ReadNext();
	return *this;
}

FastaFileIterator::ElementType& FastaFileIterator::operator * () {
	return current_;
}

FastaFileIterator::ElementType* FastaFileIterator::operator -> () {
	return &current_;
}

static void ThrowOnBadFile(const istream& stream) {
	if (stream.bad()) throw runtime_error("Error reading file.");
}

tuple<bool, string, string> FastaFileIterator::ReadNextSequence() {
	ThrowOnBadFile(*inputFile_ >> ws);

	if (inputFile_->eof()) return make_tuple(false, string(), string());

	string name;

	if (inputFile_->peek() == '>') {
		inputFile_->get();
		Readline(*inputFile_, name);
	}

	string sequence;

	for (;;) {
		ThrowOnBadFile(*inputFile_ >> ws);

		if (inputFile_->eof() || inputFile_->peek() == '>') break;

		string s;
		Readline(*inputFile_, s);

		sequence += s;
	}

	return make_tuple(true, name, sequence);
}

bool FastaFileIterator::ReadNext() {
	auto readSeq = ReadNextSequence();

	if (get<0>(readSeq)) {
		auto seq = move (get<2>(readSeq));

		seq += '\0';

		vector<pair<string, size_t>> namesAndOffsets;
		namesAndOffsets.emplace_back(make_pair(move(get<1>(readSeq)), size_t(0)));

		if (readAll_) {
			for (;;) {
				readSeq = ReadNextSequence();
				if (!get<0>(readSeq)) break;
				
				auto offset = seq.size();
				namesAndOffsets.emplace_back(make_pair(move(get<1>(readSeq)), offset));

				seq += move(get<2>(readSeq));
				seq += '\0';
			}
		}

		current_ = FastaSequence{ move(namesAndOffsets), move(seq) };
		return true;
	}
	else {
		return false;
	}
}
