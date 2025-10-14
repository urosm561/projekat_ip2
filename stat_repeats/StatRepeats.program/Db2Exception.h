#include <string>
#include <exception>

class db2exception : public std::exception {
	int code_;
	std::string msg_;

public:
	db2exception(int code) {
		code_ = code;
		msg_ = "DB2 error code " + std::to_string(code_);
	}

	int getCode() const {
		return code_;
	}

	virtual const char* what() const throw() {
		return msg_.c_str();
	}
};