#ifndef SZ3_COMPRESSOR_TYPE_ONE_HPP
#define SZ3_COMPRESSOR_TYPE_ONE_HPP

#include <cstdlib>
#include <cstring>
#include <memory>
#include <type_traits>

#include "SZ3/compressor/Compressor.hpp"
#include "SZ3/decomposition/Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Timer.hpp"

namespace SZ3 {

/// Detects the optional (non-virtual) `set_decode_bound()` an encoder may expose so the compressor can
/// hand it the exact number of bytes its encoded stream may read. Encoders without it keep whatever
/// bound their own `load()` recorded.
template <class E, class = void>
struct encoder_has_decode_bound : std::false_type {};
template <class E>
struct encoder_has_decode_bound<E, std::void_t<decltype(std::declval<E &>().set_decode_bound(size_t{}))>>
    : std::true_type {};

/**
 * SZGenericCompressor glues together decomposition, encoder, and lossless modules to form the compressor.
 * It only takes Decomposition, not Predictor.
 * @tparam T original data type
 * @tparam N original data dimension
 * @tparam Decomposition decomposition module
 * @tparam Encoder encoder module
 * @tparam Lossless lossless module
 */
template <class T, uint N, class Decomposition, class Encoder, class Lossless>
class SZGenericCompressor : public concepts::CompressorInterface<T> {
   public:
    SZGenericCompressor(Decomposition decomposition, Encoder encoder, Lossless lossless)
        : decomposition(decomposition), encoder(encoder), lossless(lossless) {
        static_assert(std::is_base_of<concepts::DecompositionInterface<T, int, N>, Decomposition>::value,
                      "must implement the frontend interface");
        static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                      "must implement the encoder interface");
        static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                      "must implement the lossless interface");
    }

    size_t compress(const Config &conf, T *data, uchar *cmpData, size_t cmpCap) override {
        std::vector<int> quant_inds = decomposition.compress(conf, data);

        if (decomposition.get_out_range().first != 0) {
            throw std::runtime_error("The output range of the decomposition must start from 0 for this compressor");
        }
        encoder.preprocess_encode(quant_inds, decomposition.get_out_range().second);
        size_t bufferSize = std::max<size_t>(
            1000, 2 * (decomposition.size_est() + encoder.size_est() + sizeof(T) * quant_inds.size()));

        // Own the scratch buffer so it is released on every path: the encoder and the lossless layer below
        // can throw (e.g. Lossless_zstd::compress throws std::length_error when the destination capacity is
        // too small for poorly-compressible data), and the caller catches and continues, so a bare free() at
        // the end leaks the buffer on each failed compression. `new uchar[]` leaves the bytes uninitialized,
        // which is what we want here -- they are all written below.
        std::unique_ptr<uchar[]> buffer_owner(new uchar[bufferSize]);
        uchar *const buffer = buffer_owner.get();
        uchar *buffer_pos = buffer;

        decomposition.save(buffer_pos);
        encoder.save(buffer_pos);

        // store the size of quant_inds is necessary as it is not always equal to conf.num
        write<size_t>(quant_inds.size(), buffer_pos);
        encoder.encode(quant_inds, buffer_pos);
        encoder.postprocess_encode();
        
        auto cmpSize = lossless.compress(buffer, buffer_pos - buffer, cmpData, cmpCap);

        return cmpSize;
    }

    T *decompress(const Config &conf, uchar const *cmpData, size_t cmpSize, T *decData) override {
        uchar *buffer = nullptr;
        // No bound is passed to the lossless layer here. The internal buffer compress() produced is
        // sized max(1000, 2 * (decomposition.size_est() + encoder.size_est() + sizeof(Q) * bins)), which
        // for a wide bin type exceeds any bound derivable from conf alone -- bounding it by
        // SZ_compress_size_bound rejects valid streams. Every module in this tree emits int bins so the
        // bound happens to hold here, but an out-of-tree module with 64-bit bins would be rejected.
        // A corrupted declared size is still caught by the zstd frame check and the size comparison in
        // Lossless_zstd::decompress, after the allocation.
        size_t bufferSize = 0;
        lossless.decompress(cmpData, cmpSize, buffer, bufferSize);

        // The lossless layer allocated `buffer` with malloc, so this one is owned with a free() deleter rather
        // than as unique_ptr<uchar[]> like the scratch buffers we allocate ourselves. Own it either way so it
        // is released on every path below - including the parsing steps that operate on untrusted data and can
        // throw before we are done with it. decompress() is reached repeatedly for corrupted blocks (fuzzing).
        std::unique_ptr<uchar, void (*)(void *)> buffer_owner(buffer, &free);

        uchar const *bufferPos = buffer;

        decomposition.load(bufferPos, bufferSize);
        encoder.load(bufferPos, bufferSize);

        size_t quant_inds_size = 0;
        // Read the count with the bounded overload so a truncated buffer can not be read past its end.
        read(quant_inds_size, bufferPos, bufferSize);
        // The count field sits between the encoder's tree and its encoded stream, so the bound load()
        // recorded for decode() is that many bytes too large. Hand the encoder the exact remaining length
        // now that the field has been consumed. Optional: encoders without the hook keep load()'s bound.
        if constexpr (encoder_has_decode_bound<Encoder>::value) {
            encoder.set_decode_bound(bufferSize);
        }
        auto quant_inds = encoder.decode(bufferPos, quant_inds_size);
        encoder.postprocess_decode();

        // The remaining work uses `quant_inds` and `decData` only, so release the internal buffer now.
        buffer_owner.reset();

        decomposition.decompress(conf, quant_inds, decData);
        return decData;
    }

   private:
    Decomposition decomposition;
    Encoder encoder;
    Lossless lossless;
};

template <class T, uint N, class Decomposition, class Encoder, class Lossless>
std::shared_ptr<SZGenericCompressor<T, N, Decomposition, Encoder, Lossless>> make_compressor_sz_generic(
    Decomposition decomposition, Encoder encoder, Lossless lossless) {
    return std::make_shared<SZGenericCompressor<T, N, Decomposition, Encoder, Lossless>>(decomposition, encoder,
                                                                                         lossless);
}

}  // namespace SZ3
#endif
