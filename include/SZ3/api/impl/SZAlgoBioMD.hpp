#ifndef SZ3_SZ_BIOMD_HPP
#define SZ3_SZ_BIOMD_HPP

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

// From a bound of 1e-2 (nm) on, most residuals are 0 and zstd at level 1 shortens their runs (by 2% at 1e-2, 23% at
// 1e-1 with one frame per chunk); below it, it gains under 2% and costs about 10% of the compression time.
template <class T, uint N>
std::shared_ptr<concepts::CompressorInterface<T>> make_compressor_biomd(const Config &conf) {
    auto encoder = SegmentedEncoder<HuffmanEncoder<int>>(biomd::NUM_STREAMS);
    if (conf.absErrorBound >= 1e-2)
        return make_compressor_sz_generic<T, N>(make_decomposition_biomd<T, N>(conf), encoder, Lossless_zstd(1));
    return make_compressor_sz_generic<T, N>(make_decomposition_biomd<T, N>(conf), encoder, Lossless_bypass());
}

template <class T, uint N>
size_t SZ_compress_bioMD(Config &conf, const T *data, uchar *cmpData, size_t cmpCap) {
    assert(N == conf.N);
    assert(conf.cmprAlgo == ALGO_BIOMD);
    if (N > 3 || conf.dims[N - 1] != 3)
        throw biomd::Fallback(ALGO_INTERP_LORENZO, "SZ3 BioMD: data must be {frames, atoms, 3}");
    calAbsErrorBound(conf, data);

    auto sz = make_compressor_biomd<T, N>(conf);
    // BIOMD only reads the data the compressor interface takes as T *
    return sz->compress(conf, const_cast<T *>(data), cmpData, cmpCap);
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
