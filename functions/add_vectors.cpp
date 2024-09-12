#include "add_vectors.h"

std::vector<double> add_vectors(const std::vector<double>& vec1, const std::vector<double>& vec2, const std::vector<double>& vec3) {
    std::vector<double> result(vec1.size());
    for (size_t i = 0; i < vec1.size(); ++i) {
        result[i] = vec1[i] + vec2[i];
        if (!vec3.empty()) {
            result[i] += vec3[i];
        }
    }
    return result;
}
