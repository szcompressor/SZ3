#ifndef SZ3_SZ_BIOMD_HPP
#define SZ3_SZ_BIOMD_HPP

#include <cmath>
#include <limits>
#include <vector>

#include "SZ3/api/impl/SZAlgoLorenzoReg.hpp"
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
//  * Other shapes, chunks of more values than an int counts (a stream holds at most one symbol per value), and
//    coordinates beyond the lattice go to LORENZO_REG with first-order Lorenzo alone, which codes coordinates best of
//    its predictors and only reads data. A dataset keeps its shape; a chunk BIOMD coded is compressed again by
//    LORENZO_REG only if appended frames take its coordinates past the lattice, which then puts the frames BIOMD
//    decoded up to twice the bound from the original data, once.
//  * NaN or Inf outside trailing fill frames, which can turn up in any appended frame: the chunk is stored losslessly,
//    which keeps the values BIOMD decoded within the bound. So is the data under a bound that is not positive and
//    finite (an ABS one, or a relative one over data with Inf), which bounds nothing.
template <class T, uint N>
size_t SZ_compress_bioMD(Config &conf, const T *data, uchar *cmpData, size_t cmpCap) {
    assert(N == conf.N);
    assert(conf.cmprAlgo == ALGO_BIOMD);
    calAbsErrorBound(conf, data);

    ALGO fallback = ALGO_LOSSLESS;
    if (conf.absErrorBound > 0 && std::isfinite(conf.absErrorBound)) {
        fallback = ALGO_LORENZO_REG;
        if (N <= 3 && conf.dims[N - 1] == 3 && conf.num <= size_t(std::numeric_limits<int>::max())) {
            auto sz = make_compressor_sz_generic<T, N>(make_decomposition_biomd<T, N>(conf),
                                                       SegmentedEncoder<HuffmanEncoder<int>>(biomd::NUM_STREAMS),
                                                       Lossless_bypass());
            try {
                return sz->compress(conf, const_cast<T *>(data), cmpData, cmpCap);
            } catch (const biomd::Fallback &e) {
                fallback = e.algo;
            }
        }
    }
    if (fallback == ALGO_LORENZO_REG) {
        conf.cmprAlgo = ALGO_LORENZO_REG;
        conf.lorenzo = true;
        conf.lorenzo2 = false;
        conf.regression = false;
        return SZ_compress_LorenzoReg<T, N>(conf, const_cast<T *>(data), cmpData, cmpCap);
    }
    conf.cmprAlgo = ALGO_LOSSLESS;
    return Lossless_zstd().compress(reinterpret_cast<const uchar *>(data), conf.num * sizeof(T), cmpData, cmpCap);
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
