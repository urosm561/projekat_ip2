#include "RepeatPrinter.h"

#pragma warning (disable:4996)

using namespace std;

ConsoleOutput::ConsoleOutput()
{
	printSequenceNames_ = false;
}

void ConsoleOutput::OutputPairs(string sequenceName, size_t x, size_t y, size_t xSeq, size_t ySeq, StringSlice repeat, StringSlice complement)
{
	if (printSequenceNames_) {
		cout << containedSequences_[xSeq] << ',' << x << ',' << (x + repeat.size() - 1) << ',' << containedSequences_[ySeq] << ',' << y << ',' << (y + repeat.size() - 1) << "," << repeat.size() << "," << repeat << "," << complement << '\n';
	}
	else {
		cout << x << "," << (x + repeat.size() - 1) << "," << y << "," << (y + repeat.size() - 1) << "," << repeat.size() << "," << repeat << "," << complement << '\n';
	}
}

void ConsoleOutput::OutputSingles(std::string sequenceName, size_t x, size_t xSeq, StringSlice repeat) {
	if (printSequenceNames_) {
		cout << containedSequences_[xSeq] << ',' << x << ',' << (x + repeat.size() - 1) << ',' << repeat.size() << "," << repeat << "," <<  '\n';
	}
	else {
		cout << x << "," << (x + repeat.size() - 1) << "," << repeat.size() << "," << repeat << '\n';
	}

}


void ConsoleOutput::BeforeAmbiguousRepeats() {
	cout << "Repeats with ambiguous letters:" << endl;
}


void ConsoleOutput::AfterEverySequence(std::vector<size_t> repeatCount)
{
	for (size_t i = 0; i < repeatCount.size(); i++) {
		if (repeatCount[i] != 0) {
			cout << "Total for length " << i << " is " << repeatCount[i] << ". " << '\n';
		}
	}
}

void ConsoleOutput::Initialize(const CommandLine& c)
{
	::WriteHeader(c, cout);

	printSequenceNames_ = !c.CombinedSequenceName().empty();
}

bool ConsoleOutput::InitializeForEverySequence(const string sequenceName, const vector<string>& containedSequences, const int version, const CommandLine& c) {
	containedSequences_ = containedSequences;

	cout << endl;
		cout << "Processing sequence " << sequenceName << "." << '\n';
	if (version != 0) {
		cout << "Version " << version << "." << '\n';
		cout << "Complete sequence name is " << sequenceName + "." << version << "." << '\n';
	}
	return false;
}


FileOutput::FileOutput(string fileName)
	{
		_writer = fopen(fileName.c_str(), "w");
		printSequenceNames_ = false;
	}

FileOutput::~FileOutput() {
		fclose(_writer);
	}

void FileOutput::OutputSingles(std::string sequenceName, size_t x, size_t xSeq, StringSlice repeat) {
	fprintf(_writer, "%s,%zu,%zu,%zu,", _sequenceName.c_str(), x, x + repeat.size() - 1, repeat.size());
	fwrite(&*repeat.begin(), sizeof(char), repeat.size(), _writer);
	fputc('\n', _writer);
}


void FileOutput::OutputPairs(string sequenceName, size_t x, size_t y, size_t xSeq, size_t ySeq, StringSlice repeat, StringSlice complement)
{
	if (x >= y) {
		swap(x, y);
		swap(xSeq, ySeq);
		swap(repeat, complement);
	}

	fprintf(_writer, "%s,%zu,%zu,%zu,%zu,%zu,", _sequenceName.c_str(), x, x + repeat.size() - 1, y, y + repeat.size() - 1, complement.size());
	fwrite(&*repeat.begin(), sizeof(char), repeat.size(), _writer);
	fputc(',', _writer);
	fwrite(&*complement.begin(), sizeof(char), complement.size(), _writer);
	fputc('\n', _writer);
}

void FileOutput::BeforeAmbiguousRepeats() {
	fprintf(_writer, "Repeats with ambiguous letters:\n");
}


void FileOutput::AfterEverySequence(std::vector<size_t> repeatCount)
{
	fputc('\n', _writer);
	for (size_t i = 0; i < repeatCount.size(); i++) {
		if (repeatCount[i] != 0) {
			fprintf(_writer, "Total for length %zu is %zu.\n", i, repeatCount[i]);
		}
	}
	fputc('\n', _writer);
}

void FileOutput::Initialize(const CommandLine& c)
{
	printSequenceNames_ = !c.CombinedSequenceName().empty();
	::WriteHeader(c, _writer);
	::WriteHeader(c, cout);
}

bool FileOutput::InitializeForEverySequence(const string sequenceName, const vector<string>& containedSequences, const int version, const CommandLine& c) {
	sequences_ = containedSequences;
	_sequenceName = sequenceName;
	//fprintf(_writer, "%s ", _sequenceName.c_str());
	//if (version != 0) {
	//	fprintf(_writer, "\nComplete sequence name is %s.%d\n", sequenceName.c_str(), version);
	//	fprintf(_writer, "Version is %d.\n", version);
	//}

	return false;
}


LoadOutput::LoadOutput(string name) {
		_writer = fopen((name + ".load").c_str(), "w");
		_stat = fopen((name + ".stat").c_str(), "w");
		_id = fopen((name + ".id").c_str(), "w");

		printSequenceNames_ = false;
	}

LoadOutput::~LoadOutput() {
		fclose(_writer);
		fclose(_stat);
		fclose(_id);
	}

void LoadOutput::BeforeAmbiguousRepeats() {
	fprintf(_writer, "Repeats with ambiguous letters:\n");
}

void LoadOutput::OutputSingles(std::string sequenceName, size_t x, size_t xSeq, StringSlice repeat) {
	fprintf(_writer, "%s,%zu,%zu,%zu,", _sequenceName.c_str(), x, x + repeat.size() - 1, repeat.size());
	fwrite(&*repeat.begin(), sizeof(char), repeat.size(), _writer);
	fputc('\n', _writer);
}


void LoadOutput::OutputPairs(string sequenceName, size_t x, size_t y, size_t xSeq, size_t ySeq, StringSlice repeat, StringSlice complement)
{
	if (x >= y) {
		swap(x, y);
		swap(xSeq, ySeq);
		swap(repeat, complement);
	}

	fprintf(_writer, "%zu,%zu,%zu,%zu,%zu,", x, x + repeat.size() - 1, y, y + repeat.size() - 1, complement.size());

	fwrite(&*repeat.begin(), sizeof(char), repeat.size(), _writer);
	fputc(',', _writer);
	fwrite(&*complement.begin(), sizeof(char), complement.size(), _writer);
	fputc('\n', _writer);
}

void LoadOutput::AfterEverySequence(vector<size_t> repeatCount)
{

	for (size_t i = 0; i < repeatCount.size(); i++) {
		if (repeatCount[i] != 0) fprintf(_stat, "Total for length %zu is %zu.\n", i, repeatCount[i]);
	}
}

void LoadOutput::Initialize(const CommandLine& c)
{
	printSequenceNames_ = !c.CombinedSequenceName().empty();
	::WriteHeader(c, _stat);
	::WriteHeader(c, cout);

}

bool LoadOutput::InitializeForEverySequence(const string sequenceName, const vector<string>& containedSequences, const int version, const CommandLine& c){
	sequences_ = containedSequences;
	_sequenceName = sequenceName;
	fprintf(_id, "%s,\n", _sequenceName.c_str());
	if (version != 0) {
		fprintf(_stat, "\nComplete sequence name is %s.%d\n", sequenceName.c_str(), version);
		fprintf(_stat, "Version is %d.\n", version);
	}

	return false;
}


template<typename K>
void IncrementCount(unordered_map<K, size_t>& map, const K& key) {
	auto it = map.find(key);
	if (it == map.end()) {
		map[key] = 1;
	}
	else {
		++it->second;
	}
}

StatisticsOutput::~StatisticsOutput() {
	size_t all = 0;
	vector<tuple<size_t, size_t, size_t>> vec;
	for (auto& x : totals_) {
		vec.push_back(make_tuple(x.first.first, x.first.second, x.second));
		all += x.first.second * x.second;
	}

	sort(vec.begin(), vec.end(), [](const tuple<size_t, size_t, size_t>& x, const tuple<size_t, size_t, size_t>& y) {
		if (get<0>(x) < get<0>(y)) return true;
		else if (get<0>(x) == get<0>(y)) return get<1>(x) < get<1>(y);
		else return false;
	});

	auto first = true;
	size_t prev = 1;
	size_t prevLength = 0;

	for (auto& x : vec) {
		if (first && get<1>(x) > 1) {
			for (size_t k = prev; k < get<1>(x); k++) {
				cout << get<0>(x) << ',' << k << ',' << 0 << endl;
			}
		}
		else if (get<1>(x) - prev > 1 && prevLength == get<0>(x)) {
			for (size_t k = prev + 1; k < get<1>(x); k++) {
				cout << get<0>(x) << ',' << k << ',' << 0 << endl;
			}
		}
		else if (get<1>(x) > 1 && prevLength != get<0>(x)) {
			for (size_t k = 1; k < get<1>(x); k++) {
				cout << get<0>(x) << ',' << k << ',' << 0 << endl;
			}
		}

		cout << get<0>(x) << ',' << get<1>(x) << ',' << get<2>(x) << endl;

		prev = get<1>(x);
		prevLength = get<0>(x);
		first = false;

	}

	cout << "Total: " << all << endl;
}

void StatisticsOutput::OutputPairs(std::string sequenceName, size_t x, size_t y, size_t xSeq, size_t ySeq, StringSlice repeat, StringSlice complement)
{
	IncrementCount(counts_, repeat);
}

void StatisticsOutput::AfterEverySequence(vector<size_t> repeatCount)
{
	cout << "Sequence done" << endl;

	for (auto& x : counts_) {
		IncrementCount(totals_, pair<size_t, size_t>(x.first.size(), x.second));
	}

	counts_.clear();
}

void StatisticsOutput::Initialize(const CommandLine& c, const size_t alphabetSize)
{
	::WriteHeader(c, cout);
}

bool StatisticsOutput::InitializeForEverySequence(const std::string sequenceName, const std::vector<std::string>& containedSequences, const int version, const CommandLine& c) {
	cout << "Processing sequence " << sequenceName << "." << '\n';
	cout << "Version " << version << "." << '\n';
	cout << "Complete sequence name is " << sequenceName + "." << version << "." << '\n';
	cout << endl;
	return false;
}

SplitOutput::~SplitOutput() {
	if (_firstSequence) {
		fclose(_writer);
		fclose(_stat);
		fclose(_id);
	}
}


SplitOutput::SplitOutput()
{
	printSequenceNames_ = false;
	_firstSequence = false;
}

void SplitOutput::OutputSingles(std::string sequenceName, size_t x, size_t xSeq, StringSlice repeat) {
	if (printSequenceNames_) {
		fprintf(_writer, "%s,%s,%zu,%zu,%zu,", _sequenceName.c_str(), sequences_[xSeq].c_str(), x, x + repeat.size() - 1, repeat.size());
	}
	else {
		fprintf(_writer, "%s,%zu,%zu,%zu,", _sequenceName.c_str(), x, x + repeat.size() - 1, repeat.size());
	}

	fwrite(&*repeat.begin(), sizeof(char), repeat.size(), _writer);
	fputc('\n', _writer);

}


void SplitOutput::OutputPairs(string sequenceName, size_t x, size_t y, size_t xSeq, size_t ySeq, StringSlice repeat, StringSlice complement)
{
	if (x >= y) {
		swap(x, y);
		swap(xSeq, ySeq);
		swap(repeat, complement);
	}

	if (printSequenceNames_) {
		fprintf(_writer, "%s,%s,%s,%zu,%zu,%zu,%zu,%zu,", _sequenceName.c_str(), sequences_[xSeq].c_str(), sequences_[ySeq].c_str(), x, x + repeat.size() - 1, y, y + repeat.size() - 1, complement.size());
	}
	else {
		fprintf(_writer, "%s,%zu,%zu,%zu,%zu,%zu,", _sequenceName.c_str(), x, x + repeat.size() - 1, y, y + repeat.size() - 1, complement.size());
	}

	fwrite(&*repeat.begin(), sizeof(char), repeat.size(), _writer);
	fputc(',', _writer);
	fwrite(&*complement.begin(), sizeof(char), complement.size(), _writer);
	fputc('\n', _writer);
}

void SplitOutput::AfterEverySequence(vector<size_t> repeatCount)
{
	for (size_t i = 0; i < repeatCount.size(); i++) {
		if (repeatCount[i] != 0)  fprintf(_stat, "Total for length %zu is %zu.\n", i, repeatCount[i]);
	}
}

void SplitOutput::BeforeAmbiguousRepeats() {
	fprintf(_writer, "Repeats with ambiguous letters:\n");
}


void SplitOutput::Initialize(const CommandLine& c)
{
	printSequenceNames_ = !c.CombinedSequenceName().empty();
	::WriteHeader(c, cout);
}

bool SplitOutput::InitializeForEverySequence(const string sequenceName, const vector<string>& containedSequences, const int version, const CommandLine& c){
	_version = version;
	_sequenceName = sequenceName;
	sequences_ = containedSequences;

	_firstSequence = true;
	_writer = fopen((sequenceName + ".load").c_str(), "w");
	_stat = fopen((sequenceName + ".stat").c_str(), "w");
	_id = fopen((sequenceName + ".id").c_str(), "w");


	::WriteHeader(c, _stat);
	if (version != 0) {
		fprintf(_stat, "Complete sequence name is %s.%d\n", sequenceName.c_str(), version);
		fprintf(_stat, "Version is %d.\n", version);
	}
	fprintf(_id, "%s", sequenceName.c_str());

	return false;
}



