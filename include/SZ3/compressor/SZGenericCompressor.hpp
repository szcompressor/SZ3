#ifndef SZ_COMPRESSOR_TYPE_ONE_HPP
#define SZ_COMPRESSOR_TYPE_ONE_HPP

#include "SZ3/compressor/Compressor.hpp"
#include "SZ3/decomposition/Decomposition.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Timer.hpp"
#include "SZ3/def.hpp"
#include <cstring>

/**
 * SZGenericCompressor glues together decomposition, encoder, and lossless modules to form the compression pipeline
 * it doesn't contains the logic to iterate through the input data. The logic is handled inside decomposition
 */

namespace SZ3 {
    template<class T, uint N, class Decomposition, class Encoder, class Lossless>
    class SZGenericCompressor : public concepts::CompressorInterface<T> {
    public:


        SZGenericCompressor(Decomposition decomposition, Encoder encoder, Lossless lossless) :
                decomposition(decomposition), encoder(encoder), lossless(lossless) {
            static_assert(std::is_base_of<concepts::DecompositionInterface<T, N>, Decomposition>::value,
                          "must implement the frontend interface");
            static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                          "must implement the encoder interface");
            static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                          "must implement the lossless interface");
        }

        size_t compress(const Config &conf, T *data, uchar *cmpData, size_t cmpCap) {

            std::vector<int> quant_inds = decomposition.compress(conf, data);

            encoder.preprocess_encode(quant_inds, decomposition.get_radius() * 2);
            size_t bufferSize = std::max<size_t>(1000, 1.2 * (decomposition.size_est() + encoder.size_est() + sizeof(T) * quant_inds.size()));

            auto buffer = (uchar *) malloc(bufferSize);
            uchar *buffer_pos = buffer;

            write(conf.num, buffer_pos);

            decomposition.save(buffer_pos);

            encoder.save(buffer_pos);
            encoder.encode(quant_inds, buffer_pos);
            encoder.postprocess_encode();

//            assert(buffer_pos - buffer < bufferSize);
            
            auto cmpSize = lossless.compress(buffer, buffer_pos - buffer, cmpData, cmpCap);
            free(buffer);
//            lossless.postcompress_data(buffer);

            return cmpSize;
        }

        T *decompress(const Config &conf, uchar const *cmpData, size_t cmpSize, T *decData) {
//            Timer timer(true);
            size_t bufferCap = conf.num * sizeof(T);
            auto buffer = (uchar *) malloc(bufferCap);
            lossless.decompress(cmpData, cmpSize, buffer, bufferCap);
            size_t remaining_length = bufferCap;
            uchar const *buffer_pos = buffer;
//            timer.stop("Lossless");
            size_t num = 0;
            read(num, buffer_pos, remaining_length);

            decomposition.load(buffer_pos, remaining_length);

            encoder.load(buffer_pos, remaining_length);

//            timer.start();
            auto quant_inds = encoder.decode(buffer_pos, num);
            encoder.postprocess_decode();
//            timer.stop("Decoder");

            free(buffer);
//            lossless.postdecompress_data(buffer);

//            timer.start();
            decomposition.decompress(conf, quant_inds, decData);
//            timer.stop("Prediction & Recover");
            return decData;
        }


    private:
        Decomposition decomposition;
        Encoder encoder;
        Lossless lossless;
    };

    template<class T, uint N, class Decomposition, class Encoder, class Lossless>
    std::shared_ptr<SZGenericCompressor<T, N, Decomposition, Encoder, Lossless>>
    make_compressor_sz_generic(Decomposition decomposition, Encoder encoder, Lossless lossless) {
        return std::make_shared<SZGenericCompressor<T, N, Decomposition, Encoder, Lossless>>(decomposition, encoder, lossless);
    }


}
#endif
