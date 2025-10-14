#ifndef STRINGSLICE_H
#define STRINGSLICE_H

#include <algorithm>
#include <string>
#include <iostream>

class StringSlice
{
	std::string::const_iterator first_;
	std::string::const_iterator last_;

public:
	StringSlice()
	{
	}

	StringSlice(const std::string &s, std::string::size_type offset, std::string::size_type len)
	{
		first_ = s.begin() + offset;
		last_ = first_ + len;
	}

	StringSlice(std::string::const_iterator x, std::string::const_iterator y)
	{
		first_ = x;
		last_ = y;
	}

	StringSlice(const StringSlice& other)
	{
		first_ = other.first_;
		last_ = other.last_;
	}

	std::string::const_iterator begin() const
	{
		return first_;
	}

	std::string::const_iterator end() const
	{
		return last_;
	}

	std::string::size_type size() const {
		return last_ - first_;
	}

	StringSlice& operator = (const StringSlice& other) {
		first_ = other.first_;
		last_ = other.last_;
		return *this;
	}

};

bool operator == (const StringSlice& x, const StringSlice& y);
bool operator < (const StringSlice& x, const StringSlice& y);
std::ostream& operator << (std::ostream& out, const StringSlice& slice);

namespace std {

	template<>
	struct hash < StringSlice > {
		size_t operator() (const StringSlice& slice) const;
	};

}


#endif