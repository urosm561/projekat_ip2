#ifndef REPEAT_COEFFICIENTS_H
#define REPEAT_COEFFICIENTS_H

#include <vector>
#include <array>

extern std::array<std::array<const double*, 17>, 17> DnaNonRepeatCoefficients;

extern std::array<std::array<const double*, 31>, 23> GenericRepeatCoefficients;
extern std::array<std::array<const double*, 33>, 23> GenericMatPalCoefficients;

#endif