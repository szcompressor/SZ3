/**
 * @file SZGenericCompressor.hpp
 * @ingroup Compressor
 */

#ifndef SZ3_COMPRESSOR_TYPE_ONE_HPP
#define SZ3_COMPRESSOR_TYPE_ONE_HPP

#include <cstring>

#include "SZ3/compressor/Compressor.hpp"
#include "SZ3/decomposition/Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/Config.hpp"

namespace SZ3 {

/**
 * @brief The default SZ3 compression pipeline.
 * 
 * This class coordinates the standard three-stage compression pipeline:
 * 1. **Decomposition** — Transforms original floating-point data into quantized integers.
 *    This is the primary extension point. Two options are available:
 *    - `BlockwiseDecomposition`: Uses a `Predictor` (e.g., Lorenzo, Regression) per block.
 *    -  `DecompositionInterface` implementations (e.g., `InterpolationDecomposition`, `ZFPDecomposition`) which operate on the whole data array.
 * 2. **Encoder** — Encodes the quantized integers (e.g., Huffman, Arithmetic).
 * 3. **Lossless** — Applies lossless compression to the encoded byte stream (e.g., Zstd).
 * 
 * @tparam T Original data type (e.g., float, double)
 * @tparam N Data dimension (e.g., 1, 2, 3)
 * @tparam Decomposition Class implementing `DecompositionInterface`
 * @tparam Encoder Class implementing `EncoderInterface`
 * @tparam Lossless Class implementing `LosslessInterface`
 */
template <class T, uint N, class Decomposition, class Encoder, class Lossless>
class SZGenericCompressor : public concepts::CompressorInterface<T> {
   public:
    /**
     * @brief Construct a new SZGenericCompressor object
     * 
     * @param decomposition Decomposition module instance
     * @param encoder Encoder module instance
     * @param lossless Lossless module instance
     */
    SZGenericCompressor(Decomposition decomposition, Encoder encoder, Lossless lossless)
        : decomposition(decomposition), encoder(encoder), lossless(lossless) {
        static_assert(std::is_base_of<concepts::DecompositionInterface<T, int, N>, Decomposition>::value,
                      "must implement the frontend interface");
        static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                      "must implement the encoder interface");
        static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                      "must implement the lossless interface");
    }

    /**
     * @brief Compress data using the configured pipeline
     * 
     * @param conf Compression configuration
     * @param data Pointer to input data
     * @param cmpData Output buffer for compressed data
     * @param cmpCap Capacity of the output buffer
     * @return size_t Size of the compressed data in bytes
     * @throw std::runtime_error
     */
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

    /**
     * @brief Decompress data using the configured pipeline
     * 
     * @param conf Compression configuration
     * @param cmpData Pointer to compressed data
     * @param cmpSize Size of compressed data
     * @param decData Buffer for decompressed data (must be pre-allocated)
     * @return T* Pointer to the decompressed data (same as decData)
     */
    T *decompress(const Config &conf, uchar const *cmpData, size_t cmpSize, T *decData) override {
        uchar *buffer = nullptr;
        size_t bufferSize = 0;
        lossless.decompress(cmpData, cmpSize, buffer, bufferSize);

        uchar const *bufferPos = buffer;

        decomposition.load(bufferPos, bufferSize);
        encoder.load(bufferPos, bufferSize);

        size_t quant_inds_size = 0;
        read(quant_inds_size, bufferPos);
        auto quant_inds = encoder.decode(bufferPos, quant_inds_size);
        encoder.postprocess_decode();

        free(buffer);

        decomposition.decompress(conf, quant_inds, decData);
        return decData;
    }

   private:
    Decomposition decomposition;
    Encoder encoder;
    Lossless lossless;
};

/**
 * @brief Factory function to create a shared_ptr to SZGenericCompressor
 * 
 * @tparam T Data type
 * @tparam N Dimension
 * @tparam Decomposition Type of decomposition module
 * @tparam Encoder Type of encoder module
 * @tparam Lossless Type of lossless module
 * @param decomposition Decomposition module instance
 * @param encoder Encoder module instance
 * @param lossless Lossless module instance
 * @return std::shared_ptr<SZGenericCompressor<T, N, Decomposition, Encoder, Lossless>> 
 */
template <class T, uint N, class Decomposition, class Encoder, class Lossless>
std::shared_ptr<SZGenericCompressor<T, N, Decomposition, Encoder, Lossless>> make_compressor_sz_generic(
    Decomposition decomposition, Encoder encoder, Lossless lossless) {
    return std::make_shared<SZGenericCompressor<T, N, Decomposition, Encoder, Lossless>>(decomposition, encoder,
                                                                                         lossless);
}

}  // namespace SZ3
#endif
