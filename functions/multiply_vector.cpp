#include "multiply_vector.h"

std::vector<double> multiply_vector(const std::vector<double>& vec, double scalar) {
    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = vec[i] * scalar;
    }
    return result;
}
