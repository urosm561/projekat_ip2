#include "StringSlice.h"

bool operator == (const StringSlice& x, const StringSlice& y)
{
	if (x.end() - x.begin() != y.end() - y.begin()) return false;

	return std::equal(x.begin(), x.end(), y.begin());
}

bool operator < (const StringSlice& x, const StringSlice& y)
{
	return lexicographical_compare(x.begin(), x.end(), y.begin(), y.end());
}

std::ostream& operator << (std::ostream& out, const StringSlice& slice)
{
	out.write(&*slice.begin(), slice.size());
	return out;
}

size_t std::hash<StringSlice>::operator() (const StringSlice& slice) const
{
	size_t h1 = 0;
	for (auto x = slice.begin(); x != slice.end(); ++x) {
		auto h2 = static_cast<size_t> (*x);
		h1 = h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h2 >> 2));
	}
	return h1;
}
