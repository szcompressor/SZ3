#ifndef SZ_Truncate_COMPRESSOR_HPP
#define SZ_Truncate_COMPRESSOR_HPP

#include "compressor/Compressor.hpp"
#include "frontend/Frontend.hpp"
#include "encoder/Encoder.hpp"
#include "lossless/Lossless.hpp"
#include "utils/FileUtil.h"
#include "utils/Config.hpp"
#include "utils/Timer.hpp"
#include "def.hpp"
#include <cstring>

namespace SZ {
    template<class T, uint N, class Lossless>
    class SZTruncateCompressor : public concepts::CompressorInterface<T> {
    public:


        SZTruncateCompressor(const Config <T, N> &conf, Lossless lossless) :
                lossless(lossless), conf(conf) {
            static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                          "must implement the lossless interface");
        }

        uchar *compress(T *data, size_t &compressed_size) {

            uchar *compressed_data = new uchar[conf.num * sizeof(T)];
            uchar *compressed_data_pos = compressed_data;
            Timer timer(true);
            uchar *data_bytes = (uchar *) data;
            uchar *data_bytes_end = data_bytes + conf.num * sizeof(T);
            while (data_bytes < data_bytes_end) {
                *compressed_data_pos++ = *data_bytes++;
                *compressed_data_pos++ = *data_bytes++;
                data_bytes += 2;
            }
            timer.stop("Prediction & Quantization");


            uchar *lossless_data = lossless.compress(compressed_data,
                                                     compressed_data_pos - compressed_data,
                                                     compressed_size);
            lossless.postcompress_data(compressed_data);
            return lossless_data;
        }

        T *decompress(uchar const *lossless_compressed_data, const size_t length) {
            size_t remaining_length = length;

            auto compressed_data = lossless.decompress(lossless_compressed_data, remaining_length);
            uchar const *compressed_data_pos = compressed_data;

            Timer timer(true);
            auto dec_data = new T[conf.num];
            uchar *data_bytes = (uchar *) dec_data;
            uchar *data_bytes_end = data_bytes + conf.num * sizeof(T);
            while (data_bytes < data_bytes_end) {
                *data_bytes++ = *compressed_data_pos++;
                *data_bytes++ = *compressed_data_pos++;
                data_bytes += 2;
            }
            lossless.postdecompress_data(compressed_data);
            timer.stop("Prediction & Recover");
            return dec_data;
        }


    private:
        Lossless lossless;
        Config <T, N> conf;
    };

    template<class T, uint N, class Lossless>
    SZTruncateCompressor<T, N, Lossless>
    make_sz_truncate_compressor(const Config <T, N> &conf, Lossless lossless) {
        return SZTruncateCompressor<T, N, Lossless>(conf, lossless);
    }
}
#endif
