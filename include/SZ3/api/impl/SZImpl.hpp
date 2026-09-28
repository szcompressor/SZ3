#ifndef SZ3_IMPL_SZ_HPP
#define SZ3_IMPL_SZ_HPP

/**
 * @file SZImpl.hpp
 * @ingroup API
 * @brief Internal entry point for compression and decompression.
 *
 * Selects between OpenMP-parallel and single-threaded execution paths
 * based on runtime configuration and OpenMP availability.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "SZ3/api/impl/SZDispatcher.hpp"
#include "SZ3/api/impl/SZImplOMP.hpp"
#include "SZ3/def.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"

namespace SZ3 {
template <class T>
using SZ_float_of = std::conditional_t<sizeof(T) <= 2, float, double>;

/**
 * @brief Internal compress entry point (dimension-templated).
 *
 * Converts integer data to `SZ_float_of<T>`, then routes to `SZ_compress_OMP` or `SZ_compress_dispatcher`
 * depending on whether OpenMP is enabled in the configuration (and available at compile time).
 *
 * @tparam T Data type
 * @tparam N Dimension
 * @param conf Configuration (may be modified to disable openmp if unavailable)
 * @param data Input data pointer
 * @param cmpData Output compressed data buffer
 * @param cmpCap Buffer capacity
 * @return Compressed size in bytes
 */
template <class T, uint N>
size_t SZ_compress_impl(Config &conf, const T *data, uchar *cmpData, size_t cmpCap) {
    // Integers are compressed as SZ_float_of<T>, which holds every value of T (8-byte ones within +-2^53), within
    // floor(bound) + 0.49, so rounding each decompressed value lands within the bound.
    if constexpr (std::is_integral<T>::value) {
        std::vector<SZ_float_of<T>> values(data, data + conf.num);
        if (sizeof(T) == 8 &&
            std::any_of(values.begin(), values.end(), [](double v) { return std::fabs(v) >= 0x1p53; }))
            throw std::invalid_argument("SZ3: 8-byte integers must lie within +-2^53, where doubles hold them exactly");
        calAbsErrorBound(conf, values.data());
        conf.absErrorBound = std::floor(conf.absErrorBound) + 0.49;
        return SZ_compress_impl<SZ_float_of<T>, N>(conf, values.data(), cmpData, cmpCap);
    } else {
        // Floating-point data is split into OpenMP chunks, or compressed as a whole by the dispatcher.
#ifndef _OPENMP
        conf.openmp = false;
#endif
        if (conf.openmp) {
            return SZ_compress_OMP<T, N>(conf, data, cmpData, cmpCap);
        } else {
            return SZ_compress_dispatcher<T, N>(conf, data, cmpData, cmpCap);
        }
    }
}

/**
 * @brief Internal decompress entry point (dimension-templated).
 *
 * Routes to `SZ_decompress_OMP` or `SZ_decompress_dispatcher` depending on
 * whether OpenMP is enabled in the configuration (and available at compile time); integer data is decompressed
 * as `SZ_float_of<T>` and rounded.
 *
 * @tparam T Data type
 * @tparam N Dimension
 * @param conf Configuration (read from compressed data header)
 * @param cmpData Compressed data buffer
 * @param cmpSize Compressed data size
 * @param decData Output decompressed data buffer
 */
template <class T, uint N>
void SZ_decompress_impl(Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
    // Integers are decompressed as SZ_float_of<T>, then rounded; floating-point data as SZ_compress_impl split it.
    if constexpr (std::is_integral<T>::value) {
        using F = SZ_float_of<T>;
        std::vector<F> values(conf.num);
        SZ_decompress_impl<F, N>(conf, cmpData, cmpSize, values.data());
        const T lo = std::numeric_limits<T>::lowest(), hi = std::numeric_limits<T>::max();
        // Clamped first: converting a value outside T is undefined. A NaN from damaged data becomes lo.
        std::transform(values.begin(), values.end(), decData, [&](F v) {
            v = std::round(v);
            return !(v > F(lo)) ? lo : v >= F(hi) ? hi : static_cast<T>(v);
        });
    } else if (conf.openmp) {
        SZ_decompress_OMP<T, N>(conf, cmpData, cmpSize, decData);
    } else {
        SZ_decompress_dispatcher<T, N>(conf, cmpData, cmpSize, decData);
    }
}

/**
 * @brief Upper bound on the compressed data size.
 *
 * Returns the minimum pre-allocated buffer size sufficient for `SZ_compress`.
 * Accounts for both the compressed payload (via ZSTD's bound) and metadata overhead.
 *
 * @tparam T Data type
 * @param conf Configuration
 * @return Maximum compressed size in bytes
 */
template <class T>
size_t SZ_compress_size_bound(const Config &conf) {
    if constexpr (std::is_integral<T>::value) {
        return SZ_compress_size_bound<SZ_float_of<T>>(conf);
    } else {
        bool omp = conf.openmp;
#ifndef _OPENMP
        omp = false;
#endif
        if (omp) {
            return 4096 + SZ_compress_size_bound_omp<T>(conf);
        } else {
            return 4096 + conf.size_est() + Lossless_zstd::compress_bound(conf.num * sizeof(T));
        }
    }
}

}  // namespace SZ3
#endif
