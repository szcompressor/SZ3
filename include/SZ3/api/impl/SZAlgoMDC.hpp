#ifndef SZ3_SZ_MDC_HPP
#define SZ3_SZ_MDC_HPP

#include <cstring>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "SZ3/api/impl/SZAlgoBioMD.hpp"
#include "SZ3/compressor/specialized/mdc/MDCodec.hpp"
#include "SZ3/def.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Statistic.hpp"

namespace SZ3 {

// ALGO_MDC takes molecular-dynamics coordinates as {frames, atoms, 3} or, for a one-frame chunk, {atoms, 3}, in nm
// with an absolute bound. Anything else -- another shape, integer data, non-finite values, a bound too fine for the
// coordinate range -- goes to ALGO_BIOMD, and the stored configuration says so.
template <class T, uint N>
size_t SZ_compress_MDC(Config &conf, const T *data, uchar *cmpData, size_t cmpCap) {
    assert(N == conf.N);
    calAbsErrorBound(conf, data);
    const bool shape = (N == 2 && conf.dims[1] == 3) || (N == 3 && conf.dims[2] == 3);
    if constexpr (std::is_floating_point<T>::value) {
        if (shape && conf.absErrorBound > 0) {
            const size_t frames = N == 3 ? conf.dims[0] : 1, atoms = conf.num / (3 * frames);
            const size_t bound = mdc::compress_bound(frames, atoms);
            try {
                if (cmpCap >= bound)
                    return mdc::compress(data, frames, atoms, conf.absErrorBound, cmpData, mdc::Options());
                thread_local std::unique_ptr<uint8_t[]> scratch;
                thread_local size_t scratchCap = 0;
                if (scratchCap < bound) {
                    scratch.reset(new uint8_t[bound]);
                    scratchCap = bound;
                }
                size_t size = mdc::compress(data, frames, atoms, conf.absErrorBound, scratch.get(), mdc::Options());
                if (size > cmpCap) throw std::length_error(SZ3_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH);
                memcpy(cmpData, scratch.get(), size);
                return size;
            } catch (std::runtime_error &) {
                // not representable by ALGO_MDC (non-finite input, bound too fine for the range): ALGO_BIOMD
            }
        }
    }
    conf.cmprAlgo = ALGO_BIOMD;
    std::vector<T> dataCopy(data, data + conf.num);  // ALGO_BIOMD overwrites its input
    return SZ_compress_bioMD<T, N>(conf, dataCopy.data(), cmpData, cmpCap);
}

template <class T, uint N>
void SZ_decompress_MDC(const Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
    assert(conf.cmprAlgo == ALGO_MDC);
    if constexpr (std::is_floating_point<T>::value) {
        const size_t frames = N == 3 ? conf.dims[0] : 1;
        if (mdc::stored_values(cmpData, cmpSize) != conf.num || conf.num != frames * (conf.num / frames))
            throw std::runtime_error("SZ3 MDC: stream does not match the configured dimensions");
        mdc::decompress(cmpData, cmpSize, decData);
    } else {
        throw std::invalid_argument("SZ3 MDC: only float and double data");
    }
}

}  // namespace SZ3
#endif
