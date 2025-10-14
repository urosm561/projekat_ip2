static char sqla_program_id[292] = 
{
 172,0,65,69,65,85,65,73,75,66,108,51,81,80,70,104,48,49,49,49,
 49,32,50,32,32,32,32,32,32,32,32,32,8,0,65,78,65,32,32,32,
 32,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,8,0,68,66,50,69,77,66,69,68,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0
};

#include "sqladef.h"

static struct sqla_runtime_info sqla_rtinfo = 
{{'S','Q','L','A','R','T','I','N'}, sizeof(wchar_t), 0, {' ',' ',' ',' '}};


static const short sqlIsLiteral   = SQL_IS_LITERAL;
static const short sqlIsInputHvar = SQL_IS_INPUT_HVAR;


#line 1 "Db2EmbeddedSqlRepeatPrinter.sqx"
#ifdef DB2EMBEDDEDOUTPUT

#include <sqlcli.h>
#include <sqlca.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include <string>
#include <exception>
#include <cstring>
#include "CommandLine.h"
#include "SharedLifetimeAllocator.h"
#include "Db2EmbeddedSqlRepeatPrinter.h"
#include "Db2Exception.h"

using namespace std;

class Db2EmbeddedSqlRepeatPrinter : public RepeatPrinter {
	string _type;
	string _databaseName;
	string _sequenceName;
	int _minLength;
	int _commitNumber;
	int _commitCounter;
	bool _logging;
	bool _probability;
	ofstream _writer;
	shared_lifetime_unordered_map<StringSlice, int> _fragments;

	
/*
EXEC SQL BEGIN DECLARE SECTION;
*/

#line 32 "Db2EmbeddedSqlRepeatPrinter.sqx"

		sqlint32 insertOrUpdateFragmentLength_;
		sqlint32 insertOrUpdateFragmentId_;
	
		sqlint32 _finished;
		sqlint32 _seqId;
		sqlint32 insertMatchFragment1Id_;
		sqlint32 insertMatchFragment2Id_;
		sqlint32 insertMatchLocation1_;
		sqlint32 insertMatchLocation2_;

		char databaseName_[20];
		char sequenceNameChars_ [81];
		sqlint32 _version;
		sqlint32 fragmentLength_;
		char typeChars_ [7];
		double sqlConfidence_;
		sqlint32 probability_;
		sqlint32 done_;

		double _confidence;

		char fragmentText_ [2000];
	
/*
EXEC SQL END DECLARE SECTION;
*/

#line 55 "Db2EmbeddedSqlRepeatPrinter.sqx"


	
/*
EXEC SQL INCLUDE SQLCA;
*/

/* SQL Communication Area - SQLCA - structures and constants */
#include "sqlca.h"
struct sqlca sqlca;


#line 57 "Db2EmbeddedSqlRepeatPrinter.sqx"


	void CheckSqlErrorCode () {
		int code = SQLCODE;
		if (code < 0) {
			throw db2exception(code);
		}
	}

public:
	Db2EmbeddedSqlRepeatPrinter(const CommandLine& c) : _fragments(10, std::hash<StringSlice>(), std::equal_to<StringSlice>(), shared_lifetime_allocator<pair<StringSlice, int>>(true)){
		_commitNumber = c.GetCommitCount();
		_logging = c.GetLogging();
		_commitCounter = 0;
		_writer.open("smartrepeats.stat");
		_databaseName = c.GetStaticSql();

		strcpy(databaseName_,_databaseName.c_str());
		
/*
EXEC SQL CONNECT TO :databaseName_;
*/

{
#line 75 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastrt(sqla_program_id, &sqla_rtinfo, &sqlca);
#line 75 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlaaloc(2,1,1,0L);
    {
      struct sqla_setdata_list sql_setdlist[1];
#line 75 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[0].sqltype = 460; sql_setdlist[0].sqllen = 20;
#line 75 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[0].sqldata = (void*)databaseName_;
#line 75 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[0].sqlind = 0L;
#line 75 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sqlasetdata(2,0,1,sql_setdlist,0L,0L);
    }
#line 75 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlacall((unsigned short)29,4,2,0,0L);
#line 75 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastop(0L);
}

#line 75 "Db2EmbeddedSqlRepeatPrinter.sqx"

		CheckSqlErrorCode();
	}

	~Db2EmbeddedSqlRepeatPrinter() {
		
/*
EXEC SQL COMMIT;
*/

{
#line 80 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastrt(sqla_program_id, &sqla_rtinfo, &sqlca);
#line 80 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlacall((unsigned short)21,0,0,0,0L);
#line 80 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastop(0L);
}

#line 80 "Db2EmbeddedSqlRepeatPrinter.sqx"

		_fragments.get_allocator().fast_free_on();
	}

	int InsertFragmentIfNeeded (const StringSlice& s) {
		auto f1Id = _fragments.find(s);
		if (f1Id != _fragments.end()){
			return f1Id->second;
		}

		size_t textLen = s.size();
		if (textLen >= sizeof (fragmentText_)) {
			textLen = sizeof(fragmentText_) - 1;
		}
		copy (s.begin (), s.begin() + textLen, fragmentText_);
		fragmentText_ [textLen] = '\0';

		insertOrUpdateFragmentLength_ = s.size();

		
/*
EXEC SQL Call InsertOrUpdateFragment (:fragmentText_, :insertOrUpdateFragmentLength_, :insertOrUpdateFragmentId_);
*/

{
#line 99 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastrt(sqla_program_id, &sqla_rtinfo, &sqlca);
#line 99 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlaaloc(2,1,2,0L);
    {
      struct sqla_setdata_list sql_setdlist[1];
#line 99 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[0].sqltype = 460; sql_setdlist[0].sqllen = 23;
#line 99 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[0].sqldata = (void*)"INSERTORUPDATEFRAGMENT";
#line 99 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[0].sqlind = 0L;
#line 99 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sqlasetdata(2,0,1,sql_setdlist,0L,0L);
    }
#line 99 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlaaloc(3,3,3,0L);
    {
      struct sqla_setdata_list sql_setdlist[3];
#line 99 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[0].sqltype = 460; sql_setdlist[0].sqllen = 2000;
#line 99 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[0].sqldata = (void*)fragmentText_;
#line 99 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[0].sqlind = 0L;
#line 99 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[1].sqltype = 496; sql_setdlist[1].sqllen = 4;
#line 99 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[1].sqldata = (void*)&insertOrUpdateFragmentLength_;
#line 99 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[1].sqlind = 0L;
#line 99 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[2].sqltype = 496; sql_setdlist[2].sqllen = 4;
#line 99 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[2].sqldata = (void*)&insertOrUpdateFragmentId_;
#line 99 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[2].sqlind = 0L;
#line 99 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sqlasetdata(3,0,3,sql_setdlist,0L,0L);
    }
#line 99 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlacall((unsigned short)77,1,2,3,0L);
#line 99 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastop(0L);
}

#line 99 "Db2EmbeddedSqlRepeatPrinter.sqx"

		CheckSqlErrorCode();

		_fragments.insert (make_pair (s, insertOrUpdateFragmentId_));
		_commitCounter++;
		if ((_commitCounter > _commitNumber) && (_commitNumber > 0)) {
			
/*
EXEC SQL COMMIT;
*/

{
#line 105 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastrt(sqla_program_id, &sqla_rtinfo, &sqlca);
#line 105 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlacall((unsigned short)21,0,0,0,0L);
#line 105 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastop(0L);
}

#line 105 "Db2EmbeddedSqlRepeatPrinter.sqx"

			CheckSqlErrorCode();
			_commitCounter = 0;
		}

		return insertOrUpdateFragmentId_;
	}

	virtual void OutputPairs(size_t x, size_t y, size_t xSeq, size_t ySeq, StringSlice repeat, StringSlice complement) override {
		insertMatchFragment1Id_ = InsertFragmentIfNeeded (repeat);
		insertMatchFragment2Id_ = InsertFragmentIfNeeded (complement);

		insertMatchLocation1_ = x;
		insertMatchLocation2_ = y;

		
/*
EXEC SQL INSERT INTO MATCH (SeqId,Id_fragment1,Id_fragment2,Location1,Location2) values (:_seqId, :insertMatchFragment1Id_, :insertMatchFragment2Id_, :insertMatchLocation1_, :insertMatchLocation2_);
*/

{
#line 120 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastrt(sqla_program_id, &sqla_rtinfo, &sqlca);
#line 120 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlaaloc(2,5,4,0L);
    {
      struct sqla_setdata_list sql_setdlist[5];
#line 120 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[0].sqltype = 496; sql_setdlist[0].sqllen = 4;
#line 120 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[0].sqldata = (void*)&_seqId;
#line 120 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[0].sqlind = 0L;
#line 120 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[1].sqltype = 496; sql_setdlist[1].sqllen = 4;
#line 120 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[1].sqldata = (void*)&insertMatchFragment1Id_;
#line 120 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[1].sqlind = 0L;
#line 120 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[2].sqltype = 496; sql_setdlist[2].sqllen = 4;
#line 120 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[2].sqldata = (void*)&insertMatchFragment2Id_;
#line 120 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[2].sqlind = 0L;
#line 120 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[3].sqltype = 496; sql_setdlist[3].sqllen = 4;
#line 120 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[3].sqldata = (void*)&insertMatchLocation1_;
#line 120 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[3].sqlind = 0L;
#line 120 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[4].sqltype = 496; sql_setdlist[4].sqllen = 4;
#line 120 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[4].sqldata = (void*)&insertMatchLocation2_;
#line 120 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[4].sqlind = 0L;
#line 120 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sqlasetdata(2,0,5,sql_setdlist,0L,0L);
    }
#line 120 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlacall((unsigned short)24,2,2,0,0L);
#line 120 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastop(0L);
}

#line 120 "Db2EmbeddedSqlRepeatPrinter.sqx"

		
		CheckSqlErrorCode();
		_commitCounter++;
		if ((_commitCounter > _commitNumber) && (_commitNumber > 0)) {
			
/*
EXEC SQL COMMIT;
*/

{
#line 125 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastrt(sqla_program_id, &sqla_rtinfo, &sqlca);
#line 125 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlacall((unsigned short)21,0,0,0,0L);
#line 125 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastop(0L);
}

#line 125 "Db2EmbeddedSqlRepeatPrinter.sqx"

			CheckSqlErrorCode();
			_commitCounter = 0;
		}
	}

	virtual void AfterEverySequence(std::vector<size_t> repeatCount) override {
		
/*
EXEC SQL UPDATE Sequence SET Finished = 1 WHERE Id = :_seqId;
*/

{
#line 132 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastrt(sqla_program_id, &sqla_rtinfo, &sqlca);
#line 132 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlaaloc(2,1,5,0L);
    {
      struct sqla_setdata_list sql_setdlist[1];
#line 132 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[0].sqltype = 496; sql_setdlist[0].sqllen = 4;
#line 132 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[0].sqldata = (void*)&_seqId;
#line 132 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[0].sqlind = 0L;
#line 132 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sqlasetdata(2,0,1,sql_setdlist,0L,0L);
    }
#line 132 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlacall((unsigned short)24,3,2,0,0L);
#line 132 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastop(0L);
}

#line 132 "Db2EmbeddedSqlRepeatPrinter.sqx"

		CheckSqlErrorCode();

		_writer << '\n';
		for(unsigned i = 0; i < repeatCount.size(); i++){
			if (repeatCount[i] != 0)  _writer << "Total for length " << i << " is "<< repeatCount[i]<< ". " << '\n';
		}
		_writer << '\n';
	}

	virtual void Initialize(const CommandLine& c, const size_t alphabetSize) override {
		_probability = c.GetProbability();
		_confidence = c.GetPVal();
		_minLength = c.GetFragmentLength();
		string s;
		c.GetIsMathematical() ? s = "n": s = "c";
		c.GetIsRepeat() ? s = "d" + s : s = "i" + s;
		_type = s;
	}

	virtual bool InitializeForEverySequence(const std::string& sequenceName, const std::vector<std::string>& containedSequences, const int version, const CommandLine& c, const size_t alphabetSize, const size_t nonAlphabetLetters, const std::string& letters, const size_t excludedCount) override {
		_sequenceName = sequenceName;
		_seqId = 0;
		string rType = "";
		_version = version;

		_writer << "Processing sequence " << _sequenceName << "." << '\n';
		_writer << "Version "<< _version << "." << '\n';
		_writer << "Complete sequence name is " << _sequenceName + "." << _version << "." << '\n';
	    _writer << '\n';

		if (!_logging) {
			
/*
EXEC SQL alter table fragment activate not logged initially;
*/

{
#line 164 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastrt(sqla_program_id, &sqla_rtinfo, &sqlca);
#line 164 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlacall((unsigned short)24,4,0,0,0L);
#line 164 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastop(0L);
}

#line 164 "Db2EmbeddedSqlRepeatPrinter.sqx"

			CheckSqlErrorCode();
			
/*
EXEC SQL alter table match activate not logged initially;
*/

{
#line 166 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastrt(sqla_program_id, &sqla_rtinfo, &sqlca);
#line 166 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlacall((unsigned short)24,5,0,0,0L);
#line 166 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastop(0L);
}

#line 166 "Db2EmbeddedSqlRepeatPrinter.sqx"

			CheckSqlErrorCode();
		}

		
/*
EXEC SQL lock table match in exclusive mode;
*/

{
#line 170 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastrt(sqla_program_id, &sqla_rtinfo, &sqlca);
#line 170 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlacall((unsigned short)24,6,0,0,0L);
#line 170 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastop(0L);
}

#line 170 "Db2EmbeddedSqlRepeatPrinter.sqx"

		CheckSqlErrorCode();
		
/*
EXEC SQL lock table fragment in exclusive mode;
*/

{
#line 172 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastrt(sqla_program_id, &sqla_rtinfo, &sqlca);
#line 172 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlacall((unsigned short)24,7,0,0,0L);
#line 172 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastop(0L);
}

#line 172 "Db2EmbeddedSqlRepeatPrinter.sqx"

		CheckSqlErrorCode();


		strcpy(sequenceNameChars_, _sequenceName.c_str());
		fragmentLength_ = c.GetFragmentLength();

		strcpy(typeChars_, _type.c_str());

		if (_probability) {
			probability_=1;
		}
		else {
			probability_=0;
		} 

		
/*
EXEC SQL Call Sequence (:sequenceNameChars_, :_version, :fragmentLength_, :_confidence, :typeChars_, :probability_, :_seqId, :done_);
*/

{
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastrt(sqla_program_id, &sqla_rtinfo, &sqlca);
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlaaloc(2,1,6,0L);
    {
      struct sqla_setdata_list sql_setdlist[1];
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[0].sqltype = 460; sql_setdlist[0].sqllen = 9;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[0].sqldata = (void*)"SEQUENCE";
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[0].sqlind = 0L;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sqlasetdata(2,0,1,sql_setdlist,0L,0L);
    }
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlaaloc(3,8,7,0L);
    {
      struct sqla_setdata_list sql_setdlist[8];
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[0].sqltype = 460; sql_setdlist[0].sqllen = 81;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[0].sqldata = (void*)sequenceNameChars_;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[0].sqlind = 0L;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[1].sqltype = 496; sql_setdlist[1].sqllen = 4;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[1].sqldata = (void*)&_version;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[1].sqlind = 0L;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[2].sqltype = 496; sql_setdlist[2].sqllen = 4;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[2].sqldata = (void*)&fragmentLength_;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[2].sqlind = 0L;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[3].sqltype = 480; sql_setdlist[3].sqllen = 8;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[3].sqldata = (void*)&_confidence;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[3].sqlind = 0L;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[4].sqltype = 460; sql_setdlist[4].sqllen = 7;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[4].sqldata = (void*)typeChars_;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[4].sqlind = 0L;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[5].sqltype = 496; sql_setdlist[5].sqllen = 4;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[5].sqldata = (void*)&probability_;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[5].sqlind = 0L;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[6].sqltype = 496; sql_setdlist[6].sqllen = 4;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[6].sqldata = (void*)&_seqId;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[6].sqlind = 0L;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[7].sqltype = 496; sql_setdlist[7].sqllen = 4;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[7].sqldata = (void*)&done_;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sql_setdlist[7].sqlind = 0L;
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
      sqlasetdata(3,0,8,sql_setdlist,0L,0L);
    }
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlacall((unsigned short)77,8,2,3,0L);
#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"
  sqlastop(0L);
}

#line 188 "Db2EmbeddedSqlRepeatPrinter.sqx"

		CheckSqlErrorCode();
		if (done_) return true;
		else return false; 
	}
};

RepeatPrinter* CreateDb2EmbeddedSqlRepeatPrinter(const CommandLine& c) {
	return new Db2EmbeddedSqlRepeatPrinter(c);
}

#ifdef _MSC_VER
#pragma comment(lib, "db2api.lib")
#endif

#endif
