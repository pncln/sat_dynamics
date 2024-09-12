#include "magnitude.h"
#include <cmath>

double magnitude(const std::vector<double>& v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}
