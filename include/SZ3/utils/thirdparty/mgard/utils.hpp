#ifndef SZ3_THIRDPARTY_MGARD_UTILS_HPP
#define SZ3_THIRDPARTY_MGARD_UTILS_HPP

#include <cstddef>
#include <vector>

namespace SZ3 {
namespace MGARD {

inline std::vector<std::vector<size_t>> init_levels(const std::vector<size_t>& dims, size_t target_level) {
    std::vector<std::vector<size_t>> level_dims;
    for (size_t i = 0; i <= target_level; i++) {
        level_dims.push_back(std::vector<size_t>(dims.size()));
    }
    for (size_t i = 0; i < dims.size(); i++) {
        size_t n = dims[i];
        for (size_t j = 0; j <= target_level; j++) {
            level_dims[target_level - j][i] = n;
            n = (n >> 1) + 1;
        }
    }
    return level_dims;
}

}  // namespace MGARD
}  // namespace SZ3

#endif
