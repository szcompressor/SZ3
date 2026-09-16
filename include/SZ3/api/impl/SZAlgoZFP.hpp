#ifndef SZ3_SZALGO_ZFP_HPP
#define SZ3_SZALGO_ZFP_HPP

/**
 * @file SZAlgoZFP.hpp
 * @ingroup API
 * @brief Compression algorithm using the ZFP compressor.
 *
 * `ALGO_ZFP` composes ZFPDecomposition with ZFPEncoder and Lossless_bypass, so wiring those
 * three together by hand gives the same pipeline.
 *
 * zfp is a block transform: 1D, 2D and 3D, float and double.
 */

#include <cassert>

#include "SZ3/compressor/SZGenericCompressor.hpp"
#include "SZ3/decomposition/ZFPDecomposition.hpp"
#include "SZ3/encoder/ZFPEncoder.hpp"
#include "SZ3/lossless/Lossless_bypass.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Statistic.hpp"

namespace SZ3 {

/// ALGO_ZFP is ZFPDecomposition -> ZFPEncoder -> bypass; the coder's output is already packed.
template <class T, uint N>
auto make_compressor_zfp(const Config &conf) {
    using Int = std::conditional_t<sizeof(T) == 4, int32_t, int64_t>;
    return make_compressor_sz_generic<T, N>(ZFPDecomposition<T, Int, N>(), ZFPEncoder<Int, N>(conf),
                                            Lossless_bypass());
}


template <class T, uint N>
size_t SZ_compress_ZFP(Config &conf, T *data, uchar *cmpData, size_t cmpCap) {
    assert(N == conf.N);
    assert(conf.cmprAlgo == ALGO_ZFP);
    calAbsErrorBound(conf, data);

    if constexpr (std::is_floating_point<T>::value && N <= 3) {
        auto zfp = make_compressor_zfp<T, N>(conf);
        return zfp->compress(conf, data, cmpData, cmpCap);
    } else if constexpr (!std::is_floating_point<T>::value) {
        throw std::invalid_argument("ZFP algorithm only supports floating-point data types.");
    } else {
        throw std::invalid_argument("ZFP algorithm supports 1D, 2D and 3D data.");
    }
}

template <class T, uint N>
void SZ_decompress_ZFP(const Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
    assert(conf.cmprAlgo == ALGO_ZFP);

    if constexpr (std::is_floating_point<T>::value && N <= 3) {
        auto zfp = make_compressor_zfp<T, N>(conf);
        zfp->decompress(conf, cmpData, cmpSize, decData);
    } else if constexpr (!std::is_floating_point<T>::value) {
        throw std::invalid_argument("ZFP algorithm only supports floating-point data types.");
    } else {
        throw std::invalid_argument("ZFP algorithm supports 1D, 2D and 3D data.");
    }
}

}  // namespace SZ3
#endif
