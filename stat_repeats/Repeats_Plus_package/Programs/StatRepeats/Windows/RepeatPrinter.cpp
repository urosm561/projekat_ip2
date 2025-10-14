#include "RepeatPrinter.h"

#ifdef ODBCOUTPUT
#include "CppOdbc.h"
#endif

#pragma warning (disable:4996)

using namespace std;

ConsoleOutput::ConsoleOutput()
{
	printSequenceNames_ = false;
}

void ConsoleOutput::OutputPairs(size_t x, size_t y, size_t xSeq, size_t ySeq, StringSlice repeat, StringSlice complement)
{
	if (printSequenceNames_) {
		cout << containedSequences_[xSeq] << ',' << x << ',' << (x + repeat.size() - 1) << ',' << containedSequences_[ySeq] << ',' << y << ',' << (y + repeat.size() - 1) << "," << repeat.size() << "," << repeat << "," << complement << '\n';
	}
	else {
		cout << x << "," << (x + repeat.size() - 1) << "," << y << "," << (y + repeat.size() - 1) << "," << repeat.size() << "," << repeat << "," << complement << '\n';
	}
}

void ConsoleOutput::AfterEverySequence(vector<size_t> repeatCount)
{
	for (size_t i = 0; i < repeatCount.size(); i++) {
		if (repeatCount[i] != 0) {
			cout << "Total for length " << i << " is " << repeatCount[i] << ". " << '\n';
		}
	}
}

void ConsoleOutput::Initialize(const CommandLine& c, const size_t alphabetSize)
{
	::WriteHeader(c, cout);

	printSequenceNames_ = !c.CombinedSequenceName().empty();
}

bool ConsoleOutput::InitializeForEverySequence(const string& sequenceName, const vector<string>& containedSequences, const int version, const CommandLine& c, const size_t alphabetSize, const size_t nonAlphabetLetters, const string& letters, const size_t excludedCount) {
	containedSequences_ = containedSequences;

	cout << endl;
		cout << "Processing sequence " << sequenceName << "." << '\n';
	if (version != 0) {
		cout << "Version " << version << "." << '\n';
		cout << "Complete sequence name is " << sequenceName + "." << version << "." << '\n';
	}
	cout << "Alphabet size is " << alphabetSize << ". " << endl;

	if (c.GetExclude() != -1) {
		cout << "Input sequence has " << excludedCount << " " << c.GetExcludeLetter() << " letters that are excluded from input sequence." << endl;
		cout << "After excluding i";
	}
	else cout << "I";

	cout << "nput sequence has " << nonAlphabetLetters  << " letters that are not in expected alphabet." << endl;
	if (nonAlphabetLetters != 0) {
		cout << "Letters are : ";
		bool first = true;
		for (char caracter: letters) {
				if (!first) cout << ", ";
				cout << caracter;
				first = false;
		}
		cout << "." << endl;
	}

	if (c.GetMotif() != "") cout<< "Finding repeats with motif " <<c.GetMotif() <<" ." << endl;
	if (!c.CombinedSequenceName().empty()) cout << "Working with msr option, finding repeats in multiple fasta files." << endl;
	cout << endl;
	cout << endl;

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

void StatisticsOutput::OutputPairs(size_t x, size_t y, size_t xSeq, size_t ySeq, StringSlice repeat, StringSlice complement)
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

bool StatisticsOutput::InitializeForEverySequence(const string& sequenceName, const vector<string>& containedSequences, const int version, const CommandLine& c, const size_t alphabetSize, const size_t nonAlphabetLetters, const string& letters, const size_t excludedCount) {
	cout << "Processing sequence " << sequenceName << "." << '\n';
	cout << "Version " << version << "." << '\n';
	cout << "Complete sequence name is " << sequenceName + "." << version << "." << '\n';
	cout << endl;
	return false;
}


FileOutput::FileOutput(string fileName) : _version(0)
	{
		_writer = fopen(fileName.c_str(), "w");
		printSequenceNames_ = false;
	}

FileOutput::~FileOutput() {
		fclose(_writer);
	}

void FileOutput::OutputPairs(size_t x, size_t y, size_t xSeq, size_t ySeq, StringSlice repeat, StringSlice complement)
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

void FileOutput::AfterEverySequence(vector<size_t> repeatCount)
{
	fputc('\n', _writer);
	for (size_t i = 0; i < repeatCount.size(); i++) {
		if (repeatCount[i] != 0) {
			fprintf(_writer, "Total for length %zu is %zu.\n", i, repeatCount[i]);
		}
	}
	fputc('\n', _writer);
}

void FileOutput::Initialize(const CommandLine& c, const size_t alphabetSize)
{
	printSequenceNames_ = !c.CombinedSequenceName().empty();
	::WriteHeader(c, _writer);
	::WriteHeader(c, cout);
}

bool FileOutput::InitializeForEverySequence(const string& sequenceName, const vector<string>& containedSequences, const int version, const CommandLine& c, const size_t alphabetSize, const size_t nonAlphabetLetters, const string& letters, const size_t excludedCount) {
	sequences_ = containedSequences;
	_sequenceName = sequenceName;
	fprintf(_writer, "%s ", _sequenceName.c_str());
	if (version != 0) {
		fprintf(_writer, "\nComplete sequence name is %s.%d\n", sequenceName.c_str(), version);
		fprintf(_writer, "Version is %d.\n", version);
	}
	fprintf(_writer, "Alphabet size is %zu.\n", alphabetSize);
	fprintf(_writer, "Input sequence has %zu letters that are not in expected alphabet.\n", nonAlphabetLetters + excludedCount);
	if (nonAlphabetLetters != 0) {
		fprintf(_writer, "Letters are : ");
		bool first = true;
		for (char caracter : letters) {
			if (!first) fprintf(_writer, ", ");
			fprintf(_writer, "%c", caracter);
			first = false;
		}
		fprintf(_writer, ".\n");
	}
	if (c.GetExclude() != -1) {
		fprintf(_writer, "Input sequence has %zu N letters that are excluded from input sequence.\n", excludedCount);
	}

	return false;
}


LoadOutput::LoadOutput(string name) : _version(0) {
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


void LoadOutput::OutputPairs(size_t x, size_t y, size_t xSeq, size_t ySeq, StringSlice repeat, StringSlice complement)
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

void LoadOutput::AfterEverySequence(vector<size_t> repeatCount)
{

	for (size_t i = 0; i < repeatCount.size(); i++) {
		if (repeatCount[i] != 0) fprintf(_stat, "Total for length %zu is %zu.\n", i, repeatCount[i]);
	}
}

void LoadOutput::Initialize(const CommandLine& c, const size_t alphabetSize)
{
	printSequenceNames_ = !c.CombinedSequenceName().empty();
	::WriteHeader(c, _stat);
	::WriteHeader(c, cout);

}

bool LoadOutput::InitializeForEverySequence(const string& sequenceName, const vector<string>& containedSequences, const int version, const CommandLine& c, const size_t alphabetSize, const size_t nonAlphabetLetters, const string& letters, const size_t excludedCount){
	sequences_ = containedSequences;
	_sequenceName = sequenceName;
	fprintf(_id, "%s,\n", _sequenceName.c_str());
	if (version != 0) {
		fprintf(_stat, "\nComplete sequence name is %s.%d\n", sequenceName.c_str(), version);
		fprintf(_stat, "Version is %d.\n", version);
	}
	fprintf(_stat, "Alphabet size is %zu.\n", alphabetSize);
	fprintf(_stat, "Input sequence has %zu letters that are not in expected alphabet.\n", nonAlphabetLetters + excludedCount);
	if (nonAlphabetLetters != 0) {
		fprintf(_stat, "Letters are : ");
		bool first = true;
		for (char caracter : letters) {
			if (!first) fprintf(_stat,", ");
			fprintf(_stat,"%c",caracter);
			first = false;
		}
		fprintf (_stat, ".\n");
	}
	if (c.GetExclude() != -1) {
		fprintf(_stat,"Input sequence has %zu N letters that are excluded from input sequence.\n", excludedCount);
	}

	return false;
}




SplitOutput::~SplitOutput() {
	if (_firstSequence) {
		fclose(_writer);
		fclose(_stat);
		fclose(_id);
	}
}


SplitOutput::SplitOutput() : _id(nullptr), _writer(nullptr), _stat(nullptr)
{
	printSequenceNames_ = false;
	_firstSequence = false;
}

void SplitOutput::OutputPairs(size_t x, size_t y, size_t xSeq, size_t ySeq, StringSlice repeat, StringSlice complement)
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


void SplitOutput::Initialize(const CommandLine& c, const size_t alphabetSize)
{
	printSequenceNames_ = !c.CombinedSequenceName().empty();
	::WriteHeader(c, cout);
}

bool SplitOutput::InitializeForEverySequence(const string& sequenceName, const vector<string>& containedSequences, const int version, const CommandLine& c, const size_t alphabetSize, const size_t nonAlphabetLetters, const string& letters, const size_t excludedCount){
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
	fprintf(_stat, "Alphabet size is %zu.\n", alphabetSize);
	fprintf(_stat, "Input sequence has %zu letters that are not in expected alphabet.\n", nonAlphabetLetters + excludedCount);
	if (nonAlphabetLetters != 0) {
		fprintf(_stat, "Letters are : ");
		bool first = true;
		for (char caracter : letters) {
			if (!first) fprintf(_stat, ", ");
			fprintf(_stat, "%c", caracter);
			first = false;
		}
		fprintf(_stat, ".\n");
	}
	if (c.GetExclude() != -1) {
		fprintf(_stat, "Input sequence has %zu N letters that are excluded from input sequence.\n", excludedCount);
	}


	return false;
}



#ifdef ODBCOUTPUT

class DataBaseOutput : public RepeatPrinter
{
	SqlEnvironment env_;
	SqlConnection con_;
	std::string type_;
	std::string sequenceName_;
	std::ofstream writer_;
	int version_;
	size_t minLength_;
	int seqId_;
	int commitNumber_;
	bool logging_;
	int commitCounter_;
	double pval_;
	bool probability_;
	shared_lifetime_unordered_map<StringSlice, int> fragments_;

	SqlStatement insertOrUpdateFragment_;
	int insertOrUpdateFragmentLength_;
	int insertOrUpdateFragmentId_;

	SqlStatement insertMatch_;
	int insertMatchFragment1Id_;
	int insertMatchFragment2Id_;
	int insertMatchLocation1_;
	int insertMatchLocation2_;

	int InsertFragmentIfNeeded(StringSlice s);
public:
	DataBaseOutput(const CommandLine& c);

	~DataBaseOutput();

	void AfterEverySequence(std::vector<size_t> repeatCount);

	virtual void Initialize(const CommandLine& c, const size_t alphabetSize);

	bool InitializeForEverySequence(const std::string& seqName, const std::vector<std::string>& containedSequences, const int version, const CommandLine& c, const size_t alphabetSize, const size_t nonAlphabetLetters, const std::string& letters, const size_t excludedCount) override;

	void OutputPairs(size_t x, size_t y, size_t xSeq, size_t ySeq, StringSlice repeat, StringSlice complement) override;
};

class DB2Output : public RepeatPrinter
{
	SqlEnvironment env_;
	SqlConnection con_;
	std::string type_;
	std::string sequenceName_;
	std::ofstream writer_;
	int version_;
	size_t minLength_;
	int seqId_;
	int commitNumber_;
	bool logging_;
	int commitCounter_;
	double pval_;
	bool probability_;
	shared_lifetime_unordered_map<StringSlice, size_t> fragments_;

	SqlStatement insertOrUpdateFragment_;
	int insertOrUpdateFragmentLength_;
	int insertOrUpdateFragmentId_;

	SqlStatement insertMatch_;
	int insertMatchFragment1Id_;
	int insertMatchFragment2Id_;
	int insertMatchLocation1_;
	int insertMatchLocation2_;

	int InsertFragmentIfNeeded(StringSlice s);

public:
	DB2Output(const CommandLine& c);
	~DB2Output();

	void AfterEverySequence(std::vector<size_t> repeatCount);

	virtual void Initialize(const CommandLine& c, const size_t alphabetSize);

	bool InitializeForEverySequence(const std::string& seqName, const std::vector<std::string>& containedSequences, const int version, const CommandLine& c, const size_t alphabetSize, const size_t nonAlphabetLetters, const std::string& letters, const size_t excludedCount) override;

	void OutputPairs(size_t x, size_t y, size_t xSeq, size_t ySeq, StringSlice repeat, StringSlice complement) override;
};

DataBaseOutput::DataBaseOutput(const CommandLine& c) : env_(), con_(env_.Connect(c.GetDatabase())), insertOrUpdateFragment_(con_.CreateStatement()),
	insertMatch_(con_.CreateStatement()), fragments_(1000, std::hash<StringSlice>(), std::equal_to<StringSlice>() , shared_lifetime_allocator<pair<StringSlice, size_t>>(true))
{
	commitNumber_ = c.GetCommitCount();
	logging_ = c.GetLogging();
	commitCounter_ = 0;
	probability_ = c.GetProbability();
	pval_ = c.GetPVal();
	minLength_ = c.GetFragmentLength();
	string s;

	c.GetIsRepeat() ? s = "d" : s = "i";
	c.GetIsMathematical() ? s = s + "n" : s = s + "c";

	type_ = s;
	writer_.open("smartrepeats.stat");
	insertOrUpdateFragment_.Prepare("{Call InsertOrUpdateFragment (?, ?, ?)}");
	insertOrUpdateFragment_.BindParameter(2, insertOrUpdateFragmentLength_);
	insertOrUpdateFragment_.BindOutputParameter(3, insertOrUpdateFragmentId_);

	insertMatch_.Prepare("INSERT INTO MATCH (SeqId,Id_fragment1,Id_fragment2,Location1,Location2) values (?,?,?,?,?)");
	insertMatch_.BindParameter(1, seqId_);
	insertMatch_.BindParameter(2, insertMatchFragment1Id_);
	insertMatch_.BindParameter(3, insertMatchFragment2Id_);
	insertMatch_.BindParameter(4, insertMatchLocation1_);
	insertMatch_.BindParameter(5, insertMatchLocation2_);
}

DataBaseOutput::~DataBaseOutput() {
	con_.Commit(false);
	fragments_.get_allocator().fast_free_on();
}

void DataBaseOutput::AfterEverySequence(vector<size_t> repeatCount)
{
	con_.Execute("UPDATE Sequence SET Finished = 1 WHERE Id = " + to_string(seqId_) + "");
	writer_ << '\n';
	for (size_t i = 0; i < repeatCount.size(); i++) {
		if (repeatCount[i] != 0)  writer_ << "Total for length " << i << " is " << repeatCount[i] << ". " << '\n';
	}
	writer_ << '\n';
}

void DataBaseOutput::Initialize(const CommandLine& c, const size_t alphabetSize)
{
	::WriteHeader(c, writer_);
}

bool DataBaseOutput::InitializeForEverySequence(const string& seqName, const vector<string>& containedSequences, const int version, const CommandLine& c, const size_t alphabetSize, const size_t nonAlphabetLetters, const string& letters, const size_t excludedCount){
	sequenceName_ = seqName;
	seqId_ = 0;

	version_ = version;
	writer_ << endl;
	writer_ << "Processing sequence " << sequenceName_ << "." << '\n';
	writer_ << "Version " << version << "." << '\n';
	writer_ << "Complete sequence name is " << sequenceName_ + "." << version << "." << '\n';
	writer_ << "Alphabet size is " << alphabetSize << ". " << endl;
	writer_ << "Input sequence has " << nonAlphabetLetters + excludedCount << " letters that are not in expected alphabet." << endl;
	if (nonAlphabetLetters != 0) {
		writer_ << "Letters are : ";
		bool first = true;
		for (char caracter : letters) {
			if (!first) writer_ << ", ";
			writer_ << caracter;
			first = false;
		}
		writer_ << "." << endl;
	}
	if (c.GetExclude() != -1) {
		writer_ << "Input sequence has " << excludedCount << " N letters that are excluded from input sequence." << endl;
	}

	int probability;
	if (probability_) {
		probability = 1;
	}
	else {
		probability = 0;
	}

	int done;
	auto sp = con_.CreateStatement();
	sp.BindOutputParameter(1, seqId_);
	sp.BindOutputParameter(2, done);
	sp.Execute("{Call Sequence ('" + sequenceName_ + "'," + to_string(version_) + "," + to_string(minLength_) + "," + to_string(pval_) + ",'" + type_ + "'," + to_string(probability) + ",?, ?)}");
	if (done) return true;
	else return false;
}

int DataBaseOutput::InsertFragmentIfNeeded(StringSlice s) {
	auto f1Id = fragments_.find(s);
	if (f1Id != fragments_.end()) {
		return f1Id->second;
	}

	string tmp(s.begin(), s.end());
	insertOrUpdateFragment_.BindParameter(1, tmp, 1000);
	insertOrUpdateFragmentLength_ = static_cast<int>(s.size());
	insertOrUpdateFragment_.ExecPrepared();
	cout << insertOrUpdateFragmentId_ << endl;
	fragments_.insert(make_pair(s, insertOrUpdateFragmentId_));
	if (commitNumber_ > 0) {
		commitCounter_++;
		if (commitCounter_ > commitNumber_) {
			con_.Commit();
			commitCounter_ = 0;
		}
	}

	return insertOrUpdateFragmentId_;
}

void DataBaseOutput::OutputPairs(size_t x, size_t y, size_t xSeq, size_t ySeq, StringSlice repeat, StringSlice complement)
{
	insertMatchFragment1Id_ = InsertFragmentIfNeeded(repeat);
	insertMatchFragment2Id_ = InsertFragmentIfNeeded(complement);

	insertMatchLocation1_ = static_cast<int> (x);
	insertMatchLocation2_ = static_cast<int> (y);

	insertMatch_.ExecPrepared();
	if (commitNumber_ > 0) {
		commitCounter_++;
		if (commitCounter_ > commitNumber_) {
			con_.Commit();
			commitCounter_ = 0;
		}
	}
}


DB2Output::DB2Output(const CommandLine& c) : env_(), con_(env_.Connect(c.GetDb2Output())), insertOrUpdateFragment_(con_.CreateStatement()),
	insertMatch_(con_.CreateStatement()), fragments_(1000, hash<StringSlice>(), equal_to<StringSlice>(), shared_lifetime_allocator<pair<StringSlice, size_t>>(true))
{
	commitNumber_ = c.GetCommitCount();
	logging_ = c.GetLogging();
	commitCounter_ = 0;
	probability_ = c.GetProbability();
	pval_ = c.GetPVal();
	minLength_ = c.GetFragmentLength();
	string s;

	c.GetIsRepeat() ? s = "d" : s = "i";
	c.GetIsMathematical() ? s = s + "n" : s = s + "c";

	type_ = s;
	writer_.open("smartrepeats.stat");
	insertOrUpdateFragment_.Prepare("{Call InsertOrUpdateFragment (?, ?, ?)}");
	insertOrUpdateFragment_.BindParameter(2, insertOrUpdateFragmentLength_);
	insertOrUpdateFragment_.BindOutputParameter(3, insertOrUpdateFragmentId_);

	insertMatch_.Prepare("INSERT INTO MATCH (SeqId,Id_fragment1,Id_fragment2,Location1,Location2) values (?,?,?,?,?)");
	insertMatch_.BindParameter(1, seqId_);
	insertMatch_.BindParameter(2, insertMatchFragment1Id_);
	insertMatch_.BindParameter(3, insertMatchFragment2Id_);
	insertMatch_.BindParameter(4, insertMatchLocation1_);
	insertMatch_.BindParameter(5, insertMatchLocation2_);
}

DB2Output::~DB2Output() {
	con_.Commit(false);
	fragments_.get_allocator().fast_free_on();
}

void DB2Output::AfterEverySequence(vector<size_t> repeatCount)
{
	con_.Execute("UPDATE Sequence SET Finished = 1 WHERE Id = " + to_string(seqId_) + "");
	writer_ << '\n';
	for (unsigned i = 0; i < repeatCount.size(); i++) {
		if (repeatCount[i] != 0)  writer_ << "Total for length " << i << " is " << repeatCount[i] << ". " << '\n';
	}
	writer_ << '\n';
}

void DB2Output::Initialize(const CommandLine& c, const size_t alphabetSize)
{
	::WriteHeader(c, writer_);
}

bool DB2Output::InitializeForEverySequence(const string& seqName, const vector<string>& containedSequences, const int version, const CommandLine& c, const size_t alphabetSize, const size_t nonAlphabetLetters, const string& letters, const size_t excludedCount){
	sequenceName_ = seqName;
	seqId_ = 0;
	version_ = version;

	writer_ << endl;
	writer_ << "Processing sequence " << sequenceName_ << "." << '\n';
	writer_ << "Version " << version << "." << '\n';
	writer_ << "Complete sequence name is " << sequenceName_ + "." << version << "." << '\n';
	writer_ << "Alphabet size is " << alphabetSize << ". " << endl;
	writer_ << "Input sequence has " << nonAlphabetLetters + excludedCount << " letters that are not in expected alphabet." << endl;
	if (nonAlphabetLetters != 0) {
		writer_ << "Letters are : ";
		bool first = true;
		for (char caracter : letters) {
			if (!first) writer_ << ", ";
			writer_ << caracter;
			first = false;
		}
		writer_ << "." << endl;
	}
	if (c.GetExclude() != -1) {
		writer_ << "Input sequence has " << excludedCount << " N letters that are excluded from input sequence." << endl;
	}

	con_.SetAutoCommit(false);
	if (!logging_) {
		con_.Execute("alter table fragment activate not logged initially");
		con_.Execute("alter table match activate not logged initially");
	}

	con_.Execute("lock table match in exclusive mode");
	con_.Execute("lock table fragment in exclusive mode");

	int probability;
	if (probability_) {
		probability = 1;
	}
	else {
		probability = 0;
	}

	int done;
	auto sp = con_.CreateStatement();
	sp.BindOutputParameter(1, seqId_);
	sp.BindOutputParameter(2, done);
	sp.Execute("{Call Sequence ('" + sequenceName_ + "'," + to_string(version_) + "," + to_string(minLength_) + "," + to_string(pval_) + ",'" + type_ + "'," + to_string(probability) + ",?, ?)}");
	if (done) return true;
	else return false;
}

int DB2Output::InsertFragmentIfNeeded(StringSlice s) {
	auto f1Id = fragments_.find(s);
	if (f1Id != fragments_.end()) {
		return static_cast<int> (f1Id->second);
	}

	string tmp(s.begin(), s.end());
	insertOrUpdateFragment_.BindParameter(1, tmp, 1000);
	insertOrUpdateFragmentLength_ = static_cast<int> (s.size());
	insertOrUpdateFragment_.ExecPrepared();
	fragments_.insert(make_pair(s, insertOrUpdateFragmentId_));
	if (commitNumber_ > 0) {
		commitCounter_++;
		if (commitCounter_ > commitNumber_) {
			con_.Commit();
			commitCounter_ = 0;
		}
	}

	return insertOrUpdateFragmentId_;
}

void DB2Output::OutputPairs(size_t x, size_t y, size_t xSeq, size_t ySeq, StringSlice repeat, StringSlice complement)
{
	insertMatchFragment1Id_ = InsertFragmentIfNeeded(repeat);
	insertMatchFragment2Id_ = InsertFragmentIfNeeded(complement);

	insertMatchLocation1_ = static_cast<int> (x);
	insertMatchLocation2_ = static_cast<int> (y);

	insertMatch_.ExecPrepared();
	if (commitNumber_ > 0) {
		commitCounter_++;
		if (commitCounter_ > commitNumber_) {
			con_.Commit();
			commitCounter_ = 0;
		}
	}
}

RepeatPrinter* CreateDataBaseOutput(const CommandLine& cl) {
	return new DataBaseOutput(cl);
}

RepeatPrinter* CreateDb2Output(const CommandLine& cl) {
	return new DB2Output(cl);
}

#endif // !NOODBC
