#pragma once

#if defined (WIN32)
#include <Windows.h>
#include <tchar.h>
#endif

#include <odbcinst.h>
#include <sql.h>
#include <sqltypes.h>
#include <iostream>
#include <sql.h>
#include <sqlext.h>
#include <stdexcept>
#include <tuple>
#include <string>


inline char* convert (SQLCHAR* x) {
	return reinterpret_cast<char*> (x);
}

inline SQLCHAR* convert (char* x) {
	return reinterpret_cast<SQLCHAR*> (x);
}

template <SQLSMALLINT TypeValue>
class SqlHandle {
public:
	SqlHandle (SqlHandle<TypeValue>&& other) {
		handle_ = other.handle_;
		other.handle_ = nullptr;
	}

	~SqlHandle() {
		if (handle_ != nullptr) {
			SQLFreeHandle (TypeValue, handle_);
		}
	}

	friend void operator >> (SQLRETURN retVal, const SqlHandle<TypeValue>& handle) {
		if (retVal < 0) {
			SQLCHAR sqlstate[1024];
			SQLCHAR message[1024];
			auto diagResult = SQLGetDiagRec(TypeValue, handle.handle_, 1, sqlstate, NULL, message, 1024, NULL);
			if(SQL_SUCCESS != diagResult) {
				throw std::logic_error ("Could not obtain SQL server diagnostics.");
			}

			throw std::logic_error (convert (message));
		}
	}

protected:
	SqlHandle () {
	}

	SqlHandle (SQLHANDLE handle) : handle_ (handle) {}

	SQLHANDLE handle_;

private:
	SqlHandle (const SqlHandle<TypeValue>& other) {}
	void operator = (const SqlHandle<TypeValue>& other) {}
};

class SqlStatement : public SqlHandle<SQL_HANDLE_STMT> {
public:
	SqlStatement (SQLHANDLE handle) : SqlHandle (handle) {}

	SqlStatement (SqlStatement&& other) : SqlHandle (std::move (other)) {}

	~SqlStatement () {
		if (handle_ != nullptr) {
			SQLFreeStmt (handle_, SQL_CLOSE | SQL_UNBIND | SQL_RESET_PARAMS);
		}
	};

	void BindColumn (SQLSMALLINT columnNumber, int& var, SQLLEN& len) const {
		SQLBindCol (handle_, columnNumber, SQL_INTEGER, &var, sizeof (int), &len) >> *this;
	}

	template<size_t N>
	void BindColumn(SQLSMALLINT columnNumber, char(&var)[N], SQLLEN& len) const {
		SQLBindCol (handle_, columnNumber, SQL_C_CHAR, var, sizeof (var), &len) >> *this;
	}

	void BindParameter (SQLSMALLINT columnNumber, int& var) const {
		SQLBindParameter(handle_, columnNumber, SQL_PARAM_INPUT, SQL_C_LONG, SQL_INTEGER, 0, 0, &var, 0, NULL) >> *this;
	}

	void BindParameter (SQLSMALLINT columnNumber, const std::string& var, int columnSize) const {
		SQLBindParameter(handle_, columnNumber, SQL_PARAM_INPUT, SQL_C_CHAR, SQL_LONGVARCHAR, columnSize, 0, const_cast<char*> (var.c_str()), var.length(), NULL) >> *this;
	}

	void BindParameter (SQLSMALLINT columnNumber, const char* begin, const char* end, int columnSize) const {
		SQLBindParameter(handle_, columnNumber, SQL_PARAM_INPUT, SQL_C_CHAR, SQL_LONGVARCHAR, columnSize, 0, const_cast<char*> (begin), end - begin, NULL) >> *this;
	}

	void BindOutputParameter(SQLSMALLINT columnNumber, int& var) const {
		SQLBindParameter(handle_, columnNumber, SQL_PARAM_OUTPUT, SQL_C_LONG, SQL_INTEGER, 0, 0, &var, 0, NULL) >> *this;
	}

	std::pair<bool, std::string> GetNullableString (int column) {
        char x [8192];
        SQLLEN xLen;
        SQLGetData (handle_, column, SQL_C_CHAR, x, sizeof (x), &xLen) >> *this;
        if (xLen == SQL_NULL_DATA) return std::pair<bool, std::string> (false, std::string());

        if (xLen <= sizeof (x)) {
			return std::pair<bool, std::string> (true, std::string (x, x + xLen));
        }
        else {
			char*c = new char [xLen];
			SQLGetData (handle_, column, SQL_C_CHAR, c, xLen, &xLen) >> *this;
			auto ret = std::pair<bool, std::string> (true, std::string (x, x + xLen));
			delete [] c;
			return ret;
        }
	}

    std::string GetString (int column) {
        char x [8192];
        SQLLEN xLen;
        SQLGetData (handle_, column, SQL_C_CHAR, x, sizeof (x), &xLen) >> *this;
        if (xLen == SQL_NULL_DATA) throw std::logic_error ("Attempting to read a nullable field with a non-null method.");

        if (xLen <= sizeof (x)) {
            return std::string (x, x + xLen);
        }
        else {
            char*c = new char [xLen];
            SQLGetData (handle_, column, SQL_C_CHAR, c, xLen, &xLen) >> *this;
            if (xLen < 0) throw std::logic_error ("Internal error reading a string field.");
            auto ret = std::string (x, x + xLen);
            delete [] c;
            return ret;
        }
    }

    std::pair<bool, int> GetNullableInt (int column) {
        int x;
        SQLLEN xLen;
        SQLGetData (handle_, column, SQL_INTEGER, &x, sizeof (x), &xLen) >> *this;
        if (xLen == SQL_NULL_DATA) {
            return std::pair<bool, int> (false, 0);
        }
        else {
            return std::pair<bool, int> (true, x);
        }
    }

    int GetInt (int column) {
        int x;
        SQLLEN xLen;
        SQLGetData (handle_, column, SQL_INTEGER, &x, sizeof (x), &xLen) >> *this;
        if (xLen == SQL_NULL_DATA) {
			throw std::logic_error ("Attempting to read a nullable field with a non-null method.");
        }
        else {
			return x;
        }
    }

	double GetDouble(int column) {
		double x;
		SQLLEN xLen;
		SQLGetData(handle_, column, SQL_C_DOUBLE, &x, sizeof (x), &xLen) >> *this;
		if (xLen == SQL_NULL_DATA) {
			throw std::logic_error("Attempting to read a nullable field with a non-null method.");
		}
		else {
			return x;
		}
	}

	std::pair<bool, double> GetNullableDouble(int column) {
		double x;
		SQLLEN xLen;
		SQLGetData(handle_, column, SQL_C_DOUBLE, &x, sizeof (x), &xLen) >> *this;
		if (xLen == SQL_NULL_DATA) {
			return std::pair<bool, double>(false, 0);
		}
		else {
			return std::pair<bool, double>(true, x);
		}
	}

	bool GetNext() {
		auto ret = SQLFetch (handle_);
		if (ret == SQL_NO_DATA) return false;
		ret >> *this;
		return true;
	}

	void ExecPrepared () {
		SQLExecute(handle_) >> *this;
	}

	void Prepare (const std::string& sql) {
		SQLPrepare (handle_, convert (const_cast<char*> (sql.c_str())), SQL_NTS) >> *this;
	}

	void Execute (const std::string& sql) {
		SQLExecDirect (handle_, convert (const_cast<char*> (sql.c_str())), SQL_NTS) >> *this;
	}
};

class SqlConnection : public SqlHandle<SQL_HANDLE_DBC> {
public:
	SqlConnection (SQLHANDLE handle) : SqlHandle (handle) {}

	SqlConnection (SqlConnection&& other) : SqlHandle (std::move (other)) {}

	SqlStatement CreateStatement () const {
		SQLHSTMT hStatement;
		SQLAllocStmt (handle_, &hStatement) >> *this;

		SqlStatement statement (hStatement);

		return statement;
	}

	SqlStatement Execute (const std::string& sql) const {
		auto statement = CreateStatement ();
		statement.Execute (sql);
		return statement;
	}

	~SqlConnection () {
		if (handle_ != nullptr) {
			SQLFreeConnect (handle_);
		}
	}

	void SetAutoCommit(bool autoCommit)
	{
		SQLSetConnectAttr(handle_, SQL_ATTR_AUTOCOMMIT, autoCommit ? (SQLPOINTER)SQL_AUTOCOMMIT_ON : (SQLPOINTER)SQL_AUTOCOMMIT_OFF, SQL_NTS) >> *this;
	}

	void Commit(bool check = true) {
		auto result = SQLEndTran(SQL_HANDLE_DBC, handle_, SQL_COMMIT);
		if (check) {
			result >> *this;
		}
	}

	void RollBack() {
		SQLEndTran(SQL_HANDLE_DBC, handle_, SQL_ROLLBACK) >> *this;
	}
};

class SqlEnvironment : SqlHandle<SQL_HANDLE_ENV> {
public :
	SqlEnvironment () {
		SQLAllocEnv (&handle_);
	}

	SqlEnvironment (SqlEnvironment&& other) : SqlHandle (std::move (other)) {}

	~SqlEnvironment () {
		if (handle_ != nullptr) {
			SQLFreeEnv (handle_);
		}
	}

	SqlConnection Connect (const std::string& connectionString) const {
		SQLHDBC hCon;
		SQLAllocConnect (handle_, &hCon) >> *this;

		SqlConnection connection (hCon);

		SQLSetConnectAttr(hCon, SQL_LOGIN_TIMEOUT, (SQLPOINTER)10, 0) >> connection;

		SQLCHAR outConnectionString [1024];

		SQLDriverConnect (hCon, NULL, convert (const_cast<char*> (connectionString.c_str())), SQL_NTS, outConnectionString, 1024, NULL, SQL_DRIVER_NOPROMPT) >> connection;

		return connection;
	}
};

