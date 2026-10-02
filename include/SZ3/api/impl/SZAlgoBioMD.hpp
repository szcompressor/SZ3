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

template <class T, uint N>
std::shared_ptr<concepts::CompressorInterface<T>> make_compressor_biomd(const Config &conf) {
    return make_compressor_sz_generic<T, N>(make_decomposition_biomd<T, N>(conf),
                                            SegmentedEncoder<HuffmanEncoder<int>>(biomd::NUM_STREAMS),
                                            Lossless_bypass());
}

// Input BIOMD does not code goes to LORENZO_REG with first-order Lorenzo alone, which codes coordinates best of its
// predictors and keeps NaN and Inf exactly; neither writes into data. A bound that is not positive and finite (an ABS
// one, or a relative one over data with Inf) bounds nothing, and the data is stored losslessly.
template <class T, uint N>
size_t SZ_compress_bioMD(Config &conf, const T *data, uchar *cmpData, size_t cmpCap) {
    assert(N == conf.N);
    assert(conf.cmprAlgo == ALGO_BIOMD);
    calAbsErrorBound(conf, data);

    if (!(conf.absErrorBound > 0) || !std::isfinite(conf.absErrorBound)) {
        conf.cmprAlgo = ALGO_LOSSLESS;
        return Lossless_zstd().compress(reinterpret_cast<const uchar *>(data), conf.num * sizeof(T), cmpData, cmpCap);
    }
    // other shapes, and chunks of more values than an int counts (a stream holds at most one symbol per value), are
    // known up front; NaN, Inf or coordinates beyond the lattice only once compress() has scanned the values
    if (N <= 3 && conf.dims[N - 1] == 3 && conf.num <= size_t(std::numeric_limits<int>::max())) {
        try {
            return make_compressor_biomd<T, N>(conf)->compress(conf, const_cast<T *>(data), cmpData, cmpCap);
        } catch (const biomd::Fallback &) {
        }
    }
    conf.cmprAlgo = ALGO_LORENZO_REG;
    conf.lorenzo = true, conf.lorenzo2 = false, conf.regression = false;
    return SZ_compress_LorenzoReg<T, N>(conf, const_cast<T *>(data), cmpData, cmpCap);
}

template <class T, uint N>
void SZ_decompress_bioMD(const Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
    assert(conf.cmprAlgo == ALGO_BIOMD);

    auto sz = make_compressor_biomd<T, N>(conf);
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
