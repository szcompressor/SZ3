#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>

#include "SZ3/api/sz.hpp"
#include "gtest/gtest.h"

// Every decompressed integer lies within floor(bound) of the input, including values at the type's limits and
// data whose range the type cannot hold.
template <class T>
void runIntegers() {
    const T lo = std::numeric_limits<T>::lowest(), hi = std::numeric_limits<T>::max();
    std::mt19937_64 g(7);
    std::vector<std::vector<T>> inputs(3, std::vector<T>(20000));
    for (size_t i = 0; i < 20000; i++) {
        inputs[0][i] = static_cast<T>(g());                                                            // full range
        inputs[1][i] = (g() & 1) ? static_cast<T>(hi - T(g() % 4)) : static_cast<T>(lo + T(g() % 4));  // both ends
        inputs[2][i] = static_cast<T>(hi - T((i * 7) % 50));  // smooth, at the top
    }
    if (sizeof(T) == 8) {  // within +-2^53
        for (auto &in : inputs)
            for (auto &v : in) v = static_cast<T>(v >> 12);
    }
    for (const auto &in : inputs) {
        for (int algo : {SZ3::ALGO_LORENZO_REG, SZ3::ALGO_INTERP_LORENZO, SZ3::ALGO_INTERP, SZ3::ALGO_NOPRED}) {
            for (auto e : {std::make_pair(SZ3::EB_ABS, 0.0), std::make_pair(SZ3::EB_ABS, 1.0),
                           std::make_pair(SZ3::EB_ABS, 3.7), std::make_pair(SZ3::EB_REL, 1e-2)}) {
                SZ3::Config conf(in.size());
                conf.cmprAlgo = algo;
                conf.errorBoundMode = e.first;
                conf.absErrorBound = e.second;
                conf.relErrorBound = e.second;
                size_t cmpSize = 0;
                std::unique_ptr<char[]> cmp(SZ_compress(conf, in.data(), cmpSize));
                SZ3::Config dconf;
                std::unique_ptr<T[]> out(SZ_decompress<T>(dconf, cmp.get(), cmpSize));
                const double range =
                    double(*std::max_element(in.begin(), in.end())) - double(*std::min_element(in.begin(), in.end()));
                const double bound = std::floor(e.first == SZ3::EB_ABS ? e.second : e.second * range);
                for (size_t i = 0; i < in.size(); i++) {
                    ASSERT_LE(std::fabs(double(out[i]) - double(in[i])), bound)
                        << "algo " << algo << " mode " << e.first << " bound " << e.second << " at " << i;
                }
            }
        }
    }
}

TEST(SZ3_IntegerData, WithinTheBound) {
    runIntegers<int8_t>();
    runIntegers<uint8_t>();
    runIntegers<int16_t>();
    runIntegers<uint16_t>();
    runIntegers<int32_t>();
    runIntegers<uint32_t>();
    runIntegers<int64_t>();
    runIntegers<uint64_t>();
}

TEST(SZ3_IntegerData, EightByteValuesPast2To53AreRefused) {
    std::vector<int64_t> in(100, int64_t(1) << 53);
    SZ3::Config conf(in.size());
    size_t cmpSize = 0;
    EXPECT_THROW(SZ_compress(conf, in.data(), cmpSize), std::invalid_argument);
}
