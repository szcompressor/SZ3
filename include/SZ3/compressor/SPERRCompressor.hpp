#ifndef SZ3_SPERR_COMPRESSOR_HPP
#define SZ3_SPERR_COMPRESSOR_HPP

#include <cstdlib>
#include <memory>
#include <type_traits>
#include <vector>

#include "SZ3/compressor/Compressor.hpp"
#include "SZ3/decomposition/SPERRDecomposition.hpp"
#include "SZ3/encoder/SPERREncoder.hpp"
#include "SZ3/lossless/Lossless_bypass.hpp"

namespace SZ3 {

/**
 * @file SPERRCompressor.hpp
 * @brief SPERR compressor wired through SZ3 module interfaces.
 *
 * The execution flow is strictly interface-driven:
 * - Compression: `decomposition.compress` -> `encoder.encode` -> `lossless.compress`
 * - Decompression: `lossless.decompress` -> `encoder.decode` -> `decomposition.decompress`
 */

template <class T, uint N, class Decomposition, class Encoder, class Lossless>
class SPERRCompressor : public concepts::CompressorInterface<T> {
   public:
    SPERRCompressor(Decomposition decomposition, Encoder encoder, Lossless lossless)
        : decomposition(decomposition), encoder(encoder), lossless(lossless) {
        static_assert(std::is_base_of<concepts::DecompositionInterface<T, uchar, N>, Decomposition>::value,
                      "SPERR decomposition must implement DecompositionInterface<T, uchar, N>.");
        static_assert(std::is_base_of<concepts::EncoderInterface<uchar>, Encoder>::value,
                      "SPERR encoder must implement EncoderInterface<uchar>.");
        static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                      "SPERR lossless must implement LosslessInterface.");
    }

    size_t compress(const Config &conf, T *data, uchar *cmpData, size_t cmpCap) override {
        auto bins = decomposition.compress(conf, data);
        encoder.preprocess_encode(bins, 0);
        auto encoded = std::vector<uchar>(bins.size() + encoder.size_est());
        auto encoded_pos = encoded.data();
        const auto encoded_len = encoder.encode(bins, encoded_pos);
        encoder.postprocess_encode();
        encoded.resize(encoded_len);

        if (encoded.size() > cmpCap) {
            throw std::length_error(SZ3_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH);
        }
        return lossless.compress(encoded.data(), encoded.size(), cmpData, cmpCap);
    }

    T *decompress(const Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) override {
        uchar *buffer = nullptr;
        size_t buffer_size = 0;
        lossless.decompress(cmpData, cmpSize, buffer, buffer_size);
        auto buffer_pos = static_cast<const uchar *>(buffer);
        encoder.preprocess_decode();
        auto bins = encoder.decode(buffer_pos, buffer_size);
        encoder.postprocess_decode();
        decomposition.decompress(conf, bins, decData);
        std::free(buffer);
        return decData;
    }

   private:
    Decomposition decomposition;
    Encoder encoder;
    Lossless lossless;
};

template <class T, uint N, class Decomposition, class Encoder, class Lossless>
std::shared_ptr<SPERRCompressor<T, N, Decomposition, Encoder, Lossless>> make_compressor_sperr(
    Decomposition decomposition, Encoder encoder, Lossless lossless) {
    return std::make_shared<SPERRCompressor<T, N, Decomposition, Encoder, Lossless>>(decomposition, encoder,
                                                                                      lossless);
}

template <class T, uint N>
std::shared_ptr<SPERRCompressor<T, N, SPERRDecomposition<T, N>, SPERREncoder<T>, Lossless_bypass>>
make_compressor_sperr() {
    return make_compressor_sperr<T, N>(SPERRDecomposition<T, N>(), SPERREncoder<T>(), Lossless_bypass());
}

}  // namespace SZ3

#endif
