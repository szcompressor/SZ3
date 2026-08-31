#ifndef SZ3_SZALGO_MGARD_HPP
#define SZ3_SZALGO_MGARD_HPP

/**
 * @file SZAlgoMGARD.hpp
 * @ingroup API
 * @brief MGARD pipeline composed from fz's modular building blocks.
 *
 * Pipeline:
 *   MGARDFusedDecomposition (multigrid + per-level LinearQuantizer
 *                               with a geometric eb schedule)
 *   -> HuffmanEncoder<int>     (one Huffman tree across all per-level indices)
 *   -> Lossless_zstd
 *
 * Each stage is a standalone module — developers can substitute any module
 * to build a different pipeline. Bundled MGARD supports floating-point
 * 1D/2D/3D data.
 */

#include "SZ3/compressor/SZGenericCompressor.hpp"
#include "SZ3/decomposition/MGARDFusedDecomposition.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Statistic.hpp"

namespace SZ3 {

template <class T, uint N>
size_t SZ_compress_MGARD(Config& conf, T* data, uchar* cmpData, size_t cmpCap) {
    assert(N == conf.N);
    assert(conf.cmprAlgo == ALGO_MGARD);
    calAbsErrorBound(conf, data);

    const int radius = conf.quantbinCnt / 2;
    auto sz = make_compressor_sz_generic<T, N>(
        MGARDFusedDecomposition<T, N>(conf.absErrorBound, radius),
        HuffmanEncoder<int>(), Lossless_zstd());
    return sz->compress(conf, data, cmpData, cmpCap);
}

template <class T, uint N>
void SZ_decompress_MGARD(const Config& conf, const uchar* cmpData, size_t cmpSize, T* decData) {
    assert(conf.cmprAlgo == ALGO_MGARD);
    auto cmpDataPos = cmpData;
    const int radius = conf.quantbinCnt / 2;
    auto sz = make_compressor_sz_generic<T, N>(
        MGARDFusedDecomposition<T, N>(conf.absErrorBound, radius),
        HuffmanEncoder<int>(), Lossless_zstd());
    sz->decompress(conf, cmpDataPos, cmpSize, decData);
}

}  // namespace SZ3

#endif
