#ifndef SZ3_SZALGO_SPERR_HPP
#define SZ3_SZALGO_SPERR_HPP

#include "SZ3/compressor/SPERRCompressor.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Statistic.hpp"

namespace SZ3 {

template <class T, uint N>
size_t SZ_compress_SPERR(Config &conf, T *data, uchar *cmpData, size_t cmpCap) {
    assert(conf.cmprAlgo == ALGO_SPERR);
    assert(conf.N == N);
    // Keep PSNR mode intact for SPERR's native q estimation path.
    // Other modes are normalized to ABS as usual.
    if (conf.errorBoundMode != EB_PSNR) {
        calAbsErrorBound(conf, data);
    }

    auto sperr = make_compressor_sperr<T, N>();
    return sperr->compress(conf, data, cmpData, cmpCap);
}

template <class T, uint N>
void SZ_decompress_SPERR(const Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
    assert(conf.cmprAlgo == ALGO_SPERR);
    auto sperr = make_compressor_sperr<T, N>();
    sperr->decompress(conf, cmpData, cmpSize, decData);
}

}  // namespace SZ3

#endif
