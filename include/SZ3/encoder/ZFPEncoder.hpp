#ifndef SZ3_ZFP_ENCODER_HPP
#define SZ3_ZFP_ENCODER_HPP

#include <vector>

#include "Encoder.hpp"
#include "SZ3/def.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/zfp/bitstream.h"
#include "SZ3/utils/zfp/intcodec04.h"
#include "SZ3/utils/zfp/intcodec16.h"
#include "SZ3/utils/zfp/intcodec64.h"

namespace SZ3 {

template <class Int, uint N>
class ZFPEncoder : public concepts::EncoderInterface<Int> {
   public:
    using UInt = std::make_unsigned_t<Int>;
    static const uint ebits = sizeof(Int) == 4 ? 8 : 11;  // number of exponent bits
    static const int ebias = (1 << (ebits - 1)) - 1;      // floating-point exponent bias

    ZFPEncoder(const Config &conf) {
        int expmin = INT_MIN;
        if (conf.absErrorBound > 0) {
            frexp(conf.absErrorBound, &expmin);
            expmin--;
        }
        minexp = std::max(expmin, std::numeric_limits<Int>::min_exponent - std::numeric_limits<Int>::digits);
    }

    size_t encode(const std::vector<Int> &data, uchar *&bytes) override {
        ZFP::MemoryBitStream stream;
        stream.open(bytes, data.size() * sizeof(Int));

        size_t n_blocks = data[0];
        uint maxprec = CHAR_BIT * sizeof(Int);
        auto emax_pos = &data[1];
        auto int_pos = &data[1 + n_blocks];
        int block_size = (N == 3 ? 64 : (N == 2 ? 16 : 4));

        stream.write(n_blocks, sizeof(size_t) * 8);
        for (auto i = 0; i < n_blocks; i++) {
            stream.write(*emax_pos + ebias, ebits);
            if constexpr (N == 1) {
                uint precision = std::min(maxprec, static_cast<uint>(std::max(0, *emax_pos - minexp + 4)));
                ZFP::IntCodec04<ZFP::MemoryBitStream, Int, UInt>::encode(stream, int_pos, 0, UINT_MAX, precision);
                int_pos += block_size;
                ++emax_pos;
            } else if constexpr (N == 2) {
                auto precision = std::min(maxprec, static_cast<uint>(std::max(0, *emax_pos - minexp + 6)));
                ZFP::IntCodec16<ZFP::MemoryBitStream, Int, UInt>::encode(stream, int_pos, 0, UINT_MAX, precision);
                int_pos += block_size;
                ++emax_pos;
            } else if constexpr (N == 3) {
                auto precision = std::min(maxprec, static_cast<uint>(std::max(0, *emax_pos - minexp + 8)));
                ZFP::IntCodec64<ZFP::MemoryBitStream, Int, UInt>::encode(stream, int_pos, 0, UINT_MAX, precision);
                int_pos += block_size;
                ++emax_pos;
            }
        }
        stream.flush();
        bytes += stream.size();
        // printf("outSize=%zu\n", stream.size());
        return stream.size();
    }

    std::vector<Int> decode(const uchar *&bytes, size_t targetLength) override {
        ZFP::MemoryBitStream stream;
        stream.open(const_cast<uchar *>(bytes), targetLength);

        size_t n_blocks = stream.read(sizeof(size_t) * 8);
        uint maxprec = CHAR_BIT * sizeof(Int);
        int block_size = (N == 3 ? 64 : (N == 2 ? 16 : 4));

        std::vector<int> output(1 + n_blocks + n_blocks * block_size);
        output[0] = n_blocks;
        auto emax_pos = &output[1];
        auto int_pos = &output[1 + n_blocks];

        for (auto i = 0; i < n_blocks; i++) {
            *emax_pos = stream.read(ebits) - ebias;
            if constexpr (N == 1) {
                uint precision = std::min(maxprec, static_cast<uint>(std::max(0, *emax_pos - minexp + 4)));
                ZFP::IntCodec04<ZFP::MemoryBitStream, Int, UInt>::decode(stream, int_pos, 0, UINT_MAX, precision);
                int_pos += block_size;
                ++emax_pos;
            } else if constexpr (N == 2) {
                auto precision = std::min(maxprec, static_cast<uint>(std::max(0, *emax_pos - minexp + 6)));
                ZFP::IntCodec16<ZFP::MemoryBitStream, Int, UInt>::decode(stream, int_pos, 0, UINT_MAX, precision);
                int_pos += block_size;
                ++emax_pos;
            } else if constexpr (N == 3) {
                auto precision = std::min(maxprec, static_cast<uint>(std::max(0, *emax_pos - minexp + 8)));
                ZFP::IntCodec64<ZFP::MemoryBitStream, Int, UInt>::decode(stream, int_pos, 0, UINT_MAX, precision);
                int_pos += block_size;
                ++emax_pos;
            }
        }
        return output;
    }
    void preprocess_encode(const std::vector<Int> &bins, int stateNum) override {}
    void postprocess_encode() override {}

    void preprocess_decode() override {}

    void postprocess_decode() override {}

    void save(uchar *&c) override {}

    void load(const uchar *&c, size_t &remaining_length) override {}

   private:
    int minexp = INT_MIN;
};

}  // namespace SZ3
#endif
