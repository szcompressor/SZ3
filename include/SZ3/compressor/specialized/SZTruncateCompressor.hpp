#ifndef SZ_Truncate_COMPRESSOR_HPP
#define SZ_Truncate_COMPRESSOR_HPP

#include <cstring>

#include "SZ3/compressor/Compressor.hpp"
#include "SZ3/decomposition/Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Timer.hpp"

/**
 * TODO
 * add truncate as Decomposition
 */
namespace SZ3 {
template <class T, uint N, class Lossless>
class SZTruncateCompressor : public concepts::CompressorInterface<T> {
   public:
    SZTruncateCompressor(const Config &conf, Lossless lossless, int byteLens)
        : lossless(lossless), conf(conf), byteLen(byteLens) {
        static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                      "must implement the lossless interface");
    }

    size_t compress(const Config &conf, T *data, uchar *cmpData, size_t cmpCap) override {
        auto buffer = static_cast<uchar *>(malloc(conf.num * sizeof(T)));
        auto buffer_pos = buffer;

        //            Timer timer(true);
        truncateArray(data, conf.num, byteLen, buffer_pos);
        //            timer.stop("Prediction & Quantization");

        auto cmpSize = lossless.compress(buffer, buffer_pos - buffer, cmpData, cmpCap);
        free(buffer);
        return cmpSize;
        //            lossless.postcompress_data(buffer);
        //            return lossless_data;
    }

    T *decompress(const Config &conf, uchar const *cmpData, size_t cmpSize, T *decData) override {
        uchar *buffer = nullptr;
        size_t bufferSize = 0;
        lossless.decompress(cmpData, cmpSize, buffer, bufferSize);
        // size_t remaining_length = bufferCap;
        uchar const *buffer_pos = buffer;

        //            Timer timer(true);
        //            auto dec_data = new T[conf.num];
        truncateArrayRecover(buffer_pos, conf.num, byteLen, decData);

        lossless.postdecompress_data(buffer);
        //            timer.stop("Prediction & Recover");
        return decData;
    }

   private:
    Lossless lossless;
    Config conf;
    int byteLen = 2;
};

template <class T, uint N, class Lossless>
SZTruncateCompressor<T, N, Lossless> make_sz_truncate_compressor(const Config &conf, Lossless lossless, int byteLens) {
    return SZTruncateCompressor<T, N, Lossless>(conf, lossless, byteLens);
}
}  // namespace SZ3
#endif
