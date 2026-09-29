#ifndef SZ3_SZ_BIOMD_HPP
#define SZ3_SZ_BIOMD_HPP

#include <cstring>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "SZ3/api/impl/SZAlgoInterp.hpp"
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

// ALGO_BIOMD takes molecular-dynamics coordinates as {frames, atoms, 3} or, for a one-frame chunk, {atoms, 3}, in nm
// with an absolute bound, and does not modify them. Other 1D to 3D data -- another shape, integers, non-finite values,
// a bound too small for the coordinate range -- goes to ALGO_INTERP_LORENZO, and the stored configuration says so.
template <class T, uint N>
size_t SZ_compress_bioMD(Config &conf, const T *data, uchar *cmpData, size_t cmpCap) {
    assert(N == conf.N);
    assert(conf.cmprAlgo == ALGO_BIOMD);
    if (N > 3) throw std::invalid_argument("SZ3 BioMD: only 1D, 2D or 3D data");
    calAbsErrorBound(conf, data);
    const bool shape = (N == 2 && conf.dims[1] == 3) || (N == 3 && conf.dims[2] == 3);
    if constexpr (std::is_floating_point<T>::value) {
        if (shape && conf.absErrorBound > 0) {
            const size_t frames = N == 3 ? conf.dims[0] : 1, atoms = conf.num / (3 * frames);
            const size_t bound = biomd::compress_bound(frames, atoms);
            try {
                if (cmpCap >= bound) return biomd::compress(data, frames, atoms, conf.absErrorBound, cmpData);
                thread_local std::vector<uint8_t> scratch;
                if (scratch.size() < bound) scratch.resize(bound);
                const size_t size = biomd::compress(data, frames, atoms, conf.absErrorBound, scratch.data());
                if (size > cmpCap) throw std::length_error(SZ3_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH);
                memcpy(cmpData, scratch.data(), size);
                return size;
            } catch (std::runtime_error &) {
                // not representable by ALGO_BIOMD
            }
        }
    }
    conf.cmprAlgo = ALGO_INTERP_LORENZO;
    std::vector<T> dataCopy(data, data + conf.num);  // ALGO_INTERP_LORENZO overwrites its input
    return SZ_compress_Interp_lorenzo<T, N>(conf, dataCopy.data(), cmpData, cmpCap);
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
