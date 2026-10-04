#ifndef SZ3_SZ_BIOMD_HPP
#define SZ3_SZ_BIOMD_HPP

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "SZ3/compressor/SZGenericCompressor.hpp"
#include "SZ3/decomposition/SZBioMDDecomposition.hpp"
#include "SZ3/decomposition/SZBioMDXtcDecomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/encoder/SegmentedEncoder.hpp"
#include "SZ3/encoder/XtcBasedEncoder.hpp"
#include "SZ3/lossless/Lossless_bypass.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Statistic.hpp"

namespace SZ3 {

// Data BIOMD does not code:
//  * Other shapes, and chunks of more values than an int counts (a stream holds at most one symbol per value): an
//    exception.
//  * NaN or Inf outside trailing fill frames, coordinates beyond the lattice, and a bound that is not positive and
//    finite (an ABS one, or a relative one over data with Inf): the chunk is stored losslessly. Frames appended to a
//    chunk BIOMD coded can bring the first two; lossless storage keeps the values BIOMD decoded, so a rewritten chunk
//    stays within the bound.
template <class T, uint N>
size_t SZ_compress_bioMD(Config &conf, const T *data, uchar *cmpData, size_t cmpCap) {
    assert(N == conf.N);
    assert(conf.cmprAlgo == ALGO_BIOMD);
    calAbsErrorBound(conf, data);

    auto lossless = [&] {
        conf.cmprAlgo = ALGO_LOSSLESS;
        return Lossless_zstd().compress(reinterpret_cast<const uchar *>(data), conf.num * sizeof(T), cmpData, cmpCap);
    };
    if (!(conf.absErrorBound > 0) || !std::isfinite(conf.absErrorBound)) return lossless();
    if (N > 3 || conf.dims[N - 1] != 3)
        throw std::invalid_argument("SZ3 ALGO_BIOMD: data must be of shape (atoms, 3) or (frames, atoms, 3)");
    if (conf.num > size_t(std::numeric_limits<int>::max()))
        throw std::invalid_argument("SZ3 ALGO_BIOMD: more values than an int counts; compress fewer frames at a time");
    auto sz =
        make_compressor_sz_generic<T, N>(make_decomposition_biomd<T, N>(conf),
                                         SegmentedEncoder<HuffmanEncoder<int>>(biomd::NUM_STREAMS), Lossless_bypass());
    try {
        return sz->compress(conf, const_cast<T *>(data), cmpData, cmpCap);
    } catch (const biomd::Fallback &) {
        return lossless();
    }
}

template <class T, uint N>
void SZ_decompress_bioMD(const Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
    assert(conf.cmprAlgo == ALGO_BIOMD);

    auto sz =
        make_compressor_sz_generic<T, N>(make_decomposition_biomd<T, N>(conf),
                                         SegmentedEncoder<HuffmanEncoder<int>>(biomd::NUM_STREAMS), Lossless_bypass());
    sz->decompress(conf, cmpData, cmpSize, decData);
}

template <class T, uint N>
size_t SZ_compress_bioMDXtcBased(Config &conf, T *data, uchar *cmpData, size_t cmpCap) {
    assert(N == conf.N);
    assert(conf.cmprAlgo == ALGO_BIOMDXTC);
    calAbsErrorBound(conf, data);

    // Not strict, the same behavior as GROMACS's xtc: rounding the stored integer back to float can put a coordinate
    // slightly past the bound.
    auto quantizer = LinearQuantizer<T>(conf.absErrorBound, XTC_radius, false);
    auto sz = make_compressor_sz_generic<T, N>(make_decomposition_biomdxtc<T, N>(conf, quantizer),
                                               XtcBasedEncoder<int>(), Lossless_bypass());
    return sz->compress(conf, data, cmpData, cmpCap);
}

template <class T, uint N>
void SZ_decompress_bioMDXtcBased(const Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
    assert(conf.cmprAlgo == ALGO_BIOMDXTC);

    LinearQuantizer<T> quantizer;
    auto sz = make_compressor_sz_generic<T, N>(make_decomposition_biomdxtc<T, N>(conf, quantizer),
                                               XtcBasedEncoder<int>(), Lossless_bypass());
    sz->decompress(conf, cmpData, cmpSize, decData);
}

}  // namespace SZ3
#endif
