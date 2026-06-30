#ifndef SZ3_COMPRESSOR_TYPE_ONE_HPP
#define SZ3_COMPRESSOR_TYPE_ONE_HPP

#include <cstdlib>
#include <cstring>
#include <memory>

#include "SZ3/compressor/Compressor.hpp"
#include "SZ3/decomposition/Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#include "zstd.h"

namespace SZ3 {
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

        auto buffer = static_cast<uchar *>(malloc(bufferSize));
        uchar *buffer_pos = buffer;

        decomposition.save(buffer_pos);
        encoder.save(buffer_pos);

        //store the size of quant_inds is necessary as it is not always equal to conf.num
        write<size_t>(quant_inds.size(), buffer_pos);
        encoder.encode(quant_inds, buffer_pos);
        encoder.postprocess_encode();
        
        auto cmpSize = lossless.compress(buffer, buffer_pos - buffer, cmpData, cmpCap);
        free(buffer);

        return cmpSize;
    }

    T *decompress(const Config &conf, uchar const *cmpData, size_t cmpSize, T *decData) override {
        uchar *buffer = nullptr;
        // The lossless layer reads the size of this internal buffer from the untrusted payload and would
        // otherwise allocate it unbounded. Bound it by the largest internal buffer this configuration could
        // have produced: during compression the buffer is losslessly (zstd) compressed, and
        // ZSTD_compressBound(B) >= B, so a block that was actually stored with this compressor satisfies
        // B < SZ_compress_size_bound = 4096 + conf.size_est() + ZSTD_compressBound(conf.num * sizeof(T)).
        // Passing this as the capacity lets the lossless decoder reject a corrupted payload that declares a
        // larger internal size before allocating it. conf.num is validated against the trusted output size by
        // the caller, so this bound can not be inflated by corrupted input.
        size_t bufferSize = 4096 + conf.size_est() + ZSTD_compressBound(conf.num * sizeof(T));
        lossless.decompress(cmpData, cmpSize, buffer, bufferSize);

        // The lossless layer allocated `buffer` with malloc. Own it with RAII so it is released on every path
        // below - including the parsing steps that operate on untrusted data and can throw before we are done
        // with it - instead of being leaked. decompress() is reached repeatedly for corrupted blocks (fuzzing).
        std::unique_ptr<uchar, void (*)(void *)> buffer_owner(buffer, &free);

        uchar const *bufferPos = buffer;

        decomposition.load(bufferPos, bufferSize);
        encoder.load(bufferPos, bufferSize);

        size_t quant_inds_size = 0;
        // Read the count with the bounded overload so a truncated buffer can not be read past its end.
        read(quant_inds_size, bufferPos, bufferSize);
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
