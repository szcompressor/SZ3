#ifndef SZ3_SZ_BIOMD_HPP
#define SZ3_SZ_BIOMD_HPP

#include <cstring>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "SZ3/compressor/SZGenericCompressor.hpp"
#include "SZ3/compressor/specialized/biomd/BioMDCodec.hpp"
#include "SZ3/decomposition/SZBioMDXtcDecomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/encoder/XtcBasedEncoder.hpp"
#include "SZ3/lossless/Lossless_bypass.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Statistic.hpp"

namespace SZ3 {

// ALGO_BIOMD takes molecular-dynamics coordinates {frames, atoms, 3}, in nm, with an absolute bound. Config drops
// dimensions of 1, so one frame, or one atom, arrives as {atoms, 3} or {frames, 3}, and one atom of one frame as {3}.
template <class T, uint N>
size_t SZ_compress_bioMD(Config &conf, const T *data, uchar *cmpData, size_t cmpCap) {
    assert(N == conf.N);
    assert(conf.cmprAlgo == ALGO_BIOMD);
    if (N > 3 || conf.dims[N - 1] != 3) throw std::invalid_argument("SZ3 BioMD: data must be {frames, atoms, 3}");
    if constexpr (!std::is_floating_point<T>::value) {
        throw std::invalid_argument("SZ3 BioMD: data must be float or double");
    } else {
        calAbsErrorBound(conf, data);
        const size_t frames = N == 3 ? conf.dims[0] : 1, atoms = N >= 2 ? conf.dims[N - 2] : 1;
        const size_t bound = biomd::compress_bound(frames, atoms);
        if (cmpCap >= bound) return biomd::compress(data, frames, atoms, conf.absErrorBound, cmpData);
        thread_local std::vector<uint8_t> scratch;
        if (scratch.size() < bound) scratch.resize(bound);
        const size_t size = biomd::compress(data, frames, atoms, conf.absErrorBound, scratch.data());
        if (size > cmpCap) throw std::length_error(SZ3_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH);
        memcpy(cmpData, scratch.data(), size);
        return size;
    }
}

template <class T, uint N>
void SZ_decompress_bioMD(const Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
    assert(conf.cmprAlgo == ALGO_BIOMD);
    if constexpr (std::is_floating_point<T>::value) {
        if (biomd::stored_values(cmpData, cmpSize) != conf.num)
            throw std::runtime_error("SZ3 BioMD: stream does not match the configured dimensions");
        biomd::decompress(cmpData, cmpSize, decData);
    } else {
        throw std::invalid_argument("SZ3 BioMD: only float and double data");
    }
}

template <class T, uint N>
size_t SZ_compress_bioMDXtcBased(Config &conf, const T *data, uchar *cmpData, size_t cmpCap) {
    assert(N == conf.N);
    assert(conf.cmprAlgo == ALGO_BIOMDXTC);
    calAbsErrorBound(conf, data);

    // Not strict, the same behavior as GROMACS's xtc: rounding the stored integer back to float can put a coordinate
    // slightly past the bound.
    auto quantizer = LinearQuantizer<T>(conf.absErrorBound, XTC_radius, false);
    auto sz = make_compressor_sz_generic<T, N>(make_decomposition_biomdxtc<T, N>(conf, quantizer),
                                               XtcBasedEncoder<int>(), Lossless_bypass());
    // SZBioMDXtcDecomposition reads the data without writing it
    return sz->compress(conf, const_cast<T *>(data), cmpData, cmpCap);
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
