#ifndef SZ_BLOCK_COMPRESSOR_HPP
#define SZ_BLOCK_COMPRESSOR_HPP

#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/compressor/Compressor.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Timer.hpp"
#include "SZ3/def.hpp"
#include <cstring>

namespace SZ {
    template<class T, uint N, class Predictor, class Encoder, class Lossless>
    class SZBlockCompressor : public concepts::CompressorInterface<T> {
    public:


        SZBlockCompressor(const Config &conf, Predictor predictor, Encoder encoder, Lossless lossless) :
                fallback_predictor(LinearQuantizer<T>(conf.absErrorBound)),
                predictor(predictor), encoder(encoder), lossless(lossless),
                block_size(conf.blockSize), num(conf.num), dims(conf.dims) {
            static_assert(std::is_base_of<concepts::PredictorInterface<T, N>, Predictor>::value,
                          "must implement the Predictor interface");
            static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                          "must implement the Encoder interface");
            static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                          "must implement the Lossless interface");

//            std::copy_n(conf.dims.begin(), N, dims.begin());
        }

        uchar *compress(const Config &conf, const T *data, size_t &compressed_size) {


            auto mddata = std::make_shared<SZ::multi_dimensional_data<T, N>>(data, dims, predictor.get_padding());
            auto block = mddata->block_iter(block_size);
            std::vector<int> quant_inds;
            quant_inds.reserve(num);

            do {
                auto isvalid = predictor.precompress(block);
                if (isvalid) {
                    predictor.compress(block, quant_inds);
                } else {
                    fallback_predictor.compress(block, quant_inds);
                }

            } while (block.next());


            encoder.preprocess_encode(quant_inds, 0);
            size_t bufferSize = 1.2 * (encoder.size_est() + sizeof(T) * quant_inds.size() + predictor.size_est());

            uchar *buffer = new uchar[bufferSize];
            uchar *buffer_pos = buffer;

            fallback_predictor.save(buffer_pos);
            predictor.save(buffer_pos);

            encoder.save(buffer_pos);
            encoder.encode(quant_inds, buffer_pos);
            encoder.postprocess_encode();
            assert(buffer_pos - buffer < bufferSize);

            uchar *lossless_data = lossless.compress(buffer, buffer_pos - buffer, compressed_size);
            lossless.postcompress_data(buffer);

            return lossless_data;
        }

        T *decompress(uchar const *cmpData, const size_t &cmpSize, size_t num) {
            T *dec_data = new T[num];
            return decompress(cmpData, cmpSize, dec_data);
        }

        T *decompress(uchar const *cmpData, const size_t &cmpSize, T *decData) {
            size_t remaining_length = cmpSize;

            auto compressed_data = lossless.decompress(cmpData, remaining_length);
            uchar const *compressed_data_pos = compressed_data;

            fallback_predictor.load(compressed_data_pos, remaining_length);

            predictor.load(compressed_data_pos, remaining_length);

            encoder.load(compressed_data_pos, remaining_length);
            auto quant_inds = encoder.decode(compressed_data_pos, num);
            encoder.postprocess_decode();

            lossless.postdecompress_data(compressed_data);

            int *quant_inds_pos = &quant_inds[0];

            auto mddata = std::make_shared<SZ::multi_dimensional_data<T, N>>(decData, dims, predictor.get_padding());
            auto block = mddata->block_iter(block_size);
            do {
                auto isvalid = predictor.predecompress(block);
                if (isvalid) {
                    predictor.decompress(block, quant_inds_pos);
                } else {
                    fallback_predictor.decompress(block, quant_inds_pos);
                }
            } while (block.next());

            mddata->copy_data_out(decData);

            return decData;
        }


    private:
        Predictor predictor;
        LorenzoPredictor<T, N, 1, LinearQuantizer<T>> fallback_predictor;
        Encoder encoder;
        Lossless lossless;
        uint block_size;
        size_t num;
        std::vector<size_t> dims;
    };

    template<class T, uint N, class Predictor, class Encoder, class Lossless>
    std::shared_ptr<SZBlockCompressor<T, N, Predictor, Encoder, Lossless>>
    make_sz_block_compressor(const Config &conf, Predictor predictor, Encoder encoder, Lossless lossless) {
        return std::make_shared<SZBlockCompressor<T, N, Predictor, Encoder, Lossless>>(conf, predictor, encoder, lossless);
    }


}
#endif
