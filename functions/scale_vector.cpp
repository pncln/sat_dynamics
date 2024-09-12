#include "scale_vector.h"

std::vector<double> scale_vector(const std::vector<double>& v, double scalar) {
    return {v[0] * scalar, v[1] * scalar, v[2] * scalar};
}

