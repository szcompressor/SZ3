#ifndef SZ3_SZALGO_ZFP_HPP
#define SZ3_SZALGO_ZFP_HPP

/**
 * @file SZAlgoZFP.hpp
 * @ingroup API
 * @brief Compression algorithm using the ZFP compressor.
 *
 * Provides a standalone ZFP compression path, selectable via `ALGO_ZFP`.
 * 
 * Note on execution workflow:
 * This algorithm wrapper invokes the native ZFP compressor API directly. 
 * We also have ZFP modules, if you want decomposition then encoding, you can use ZFPDecomposition and ZFPEncoder yourself.
 * 
 * ZFP uses a block-based transform approach and supports 1D, 2D, and 3D data.
 * Only floating-point types (float, double) are supported.
 */

#include "SZ3/compressor/ZFPCompressor.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Statistic.hpp"

namespace SZ3 {

template <class T, uint N>
size_t SZ_compress_ZFP(Config &conf, T *data, uchar *cmpData, size_t cmpCap) {
    assert(N == conf.N);
    assert(conf.cmprAlgo == ALGO_ZFP);
    calAbsErrorBound(conf, data);

    if constexpr (std::is_floating_point<T>::value) {
        auto zfp = make_compressor_zfp<T, N>();
        return zfp->compress(conf, data, cmpData, cmpCap);
    } else {
        throw std::invalid_argument("ZFP algorithm only supports floating-point data types.");
    }
}

template <class T, uint N>
void SZ_decompress_ZFP(const Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
    assert(conf.cmprAlgo == ALGO_ZFP);

    if constexpr (std::is_floating_point<T>::value) {
        auto zfp = make_compressor_zfp<T, N>();
        zfp->decompress(conf, cmpData, cmpSize, decData);
    } else {
        throw std::invalid_argument("ZFP algorithm only supports floating-point data types.");
    }
}

}  // namespace SZ3
#endif
