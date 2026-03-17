#ifndef SZ3_SZALGO_SPERR_HPP
#define SZ3_SZALGO_SPERR_HPP

#include <cstdint>

#include "SZ3/compressor/SZGenericCompressor.hpp"
#include "SZ3/decomposition/SPERRDecomposition.hpp"
#include "SZ3/encoder/SPERREncoder.hpp"
#include "SZ3/lossless/Lossless_bypass.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Statistic.hpp"

namespace SZ3 {

inline SZ3::SPERR::dims_type make_sperr_dims_from_conf(const Config &conf) {
    if (conf.dims.size() == 3) {
        // SZ3 dims: [z, y, x], SPERR dims: [x, y, z].
        return SZ3::SPERR::dims_type{conf.dims[2], conf.dims[1], conf.dims[0]};
    }
    if (conf.dims.size() == 2) {
        return SZ3::SPERR::dims_type{conf.dims[1], conf.dims[0], 1};
    }
    if (conf.dims.size() == 1) {
        return SZ3::SPERR::dims_type{conf.dims[0], 1, 1};
    }
    return SZ3::SPERR::dims_type{0, 0, 0};
}

template <class T, uint N>
size_t SZ_compress_SPERR(Config &conf, T *data, uchar *cmpData, size_t cmpCap) {
    assert(conf.cmprAlgo == ALGO_SPERR);
    assert(conf.N == N);
    // Preserve PSNR mode. Normalize other modes to ABS.
    if (conf.errorBoundMode != EB_PSNR) {
        calAbsErrorBound(conf, data);
    }

    auto sperr = make_compressor_sz_generic<T, N>(SPERRDecomposition<T, N>(),
                                                  SPERREncoder<int64_t, N>(make_sperr_dims_from_conf(conf)),
                                                  Lossless_bypass());
    return sperr->compress(conf, data, cmpData, cmpCap);
}

template <class T, uint N>
void SZ_decompress_SPERR(const Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
    assert(conf.cmprAlgo == ALGO_SPERR);
    auto sperr = make_compressor_sz_generic<T, N>(SPERRDecomposition<T, N>(),
                                                  SPERREncoder<int64_t, N>(make_sperr_dims_from_conf(conf)),
                                                  Lossless_bypass());
    sperr->decompress(conf, cmpData, cmpSize, decData);
}

}  // namespace SZ3

#endif
