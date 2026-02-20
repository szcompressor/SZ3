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

#include "SZ3/api/impl/SZDispatcher.hpp"
#include "SZ3/api/impl/SZImplOMP.hpp"
#include "SZ3/def.hpp"

namespace SZ3 {

/**
 * @brief Internal compress entry point (dimension-templated).
 *
 * Routes to `SZ_compress_OMP` or `SZ_compress_dispatcher` depending on
 * whether OpenMP is enabled in the configuration (and available at compile time).
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
#ifndef _OPENMP
    conf.openmp = false;
#endif
    if (conf.openmp) {
        return SZ_compress_OMP<T, N>(conf, data, cmpData, cmpCap);
    } else {
        return SZ_compress_dispatcher<T, N>(conf, data, cmpData, cmpCap);
    }
}

/**
 * @brief Internal decompress entry point (dimension-templated).
 *
 * Routes to `SZ_decompress_OMP` or `SZ_decompress_dispatcher` depending on
 * whether OpenMP is enabled in the configuration (and available at compile time).
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
#ifndef _OPENMP
    conf.openmp = false;
#endif
    if (conf.openmp) {
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
    bool omp = conf.openmp;
#ifndef _OPENMP
    omp = false;
#endif
    if (omp) {
        return 4096 + SZ_compress_size_bound_omp<T>(conf);
    } else {
        return 4096 + conf.size_est() + ZSTD_compressBound(conf.num * sizeof(T));
    }
}

}  // namespace SZ3
#endif
