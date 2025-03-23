#ifndef SZ_COMPRESSOR_LOSSLESS_HPP
#define SZ_COMPRESSOR_LOSSLESS_HPP

#include <cstring>

#include "SZ3/compressor/Compressor.hpp"
#include "SZ3/decomposition/Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Timer.hpp"

namespace SZ3 {
/**
 * SZEncodingLosslessCompressor glues together encoder, and lossless modules to form the compressor.
 * @tparam Encoder encoder module
 * @tparam Lossless lossless module
 */
template <  class Encoder, class Lossless>
class SZEncodingLosslessCompressor  {
   public:
    SZEncodingLosslessCompressor(Encoder encoder, Lossless lossless, int radius = 32768)
        : encoder(encoder), lossless(lossless), radius (radius) {
        static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                      "must implement the encoder interface");
        static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                      "must implement the lossless interface");
    }

    size_t compress( std::vector<int> &quant_inds, uchar *cmpData, size_t cmpCap)  {
        //std::vector<int> quant_inds = decomposition.compress(conf, data);

        encoder.preprocess_encode(quant_inds, radius * 2);
        size_t bufferSize = std::max<size_t>(
            1000, 2.0 * ( encoder.size_est() + sizeof(int) * quant_inds.size()));

        auto buffer = static_cast<uchar *>(malloc(bufferSize));
        uchar *buffer_pos = buffer;

        encoder.save(buffer_pos);

        //store the size of quant_inds is necessary as it is not always equal to conf.num
        write<size_t>(quant_inds.size(), buffer_pos);
        encoder.encode(quant_inds, buffer_pos);
        encoder.postprocess_encode();
        
        auto cmpSize = lossless.compress(buffer, buffer_pos - buffer, cmpData, cmpCap);
        free(buffer);

        return cmpSize;
    }

     std::vector<int> decompress( uchar const *cmpData, size_t cmpSize)  {
        uchar *buffer = nullptr;
        size_t bufferSize = 0;
        lossless.decompress(cmpData, cmpSize, buffer, bufferSize);

        uchar const *bufferPos = buffer;
        encoder.load(bufferPos, bufferSize);

        size_t quant_inds_size = 0;
        read(quant_inds_size, bufferPos);
        auto quant_inds = encoder.decode(bufferPos, quant_inds_size);
        encoder.postprocess_decode();

        free(buffer);
        return quant_inds;
    }

   private:
    Encoder encoder;
    Lossless lossless;
    int radius = 32768;
};

template < class Encoder, class Lossless>
std::shared_ptr<SZEncodingLosslessCompressor<Encoder, Lossless>> make_compressor_sz_encodinglossless(
    Encoder encoder, Lossless lossless) {
    return std::make_shared<SZEncodingLosslessCompressor<Encoder, Lossless>>(encoder,lossless);
}

}  // namespace SZ3
#endif
