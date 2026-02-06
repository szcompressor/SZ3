
/**
 * @file BitshuffleEncoder.hpp
 * @ingroup Encoder
 */

#ifndef SZ3_BITSHUFFLEENCODER_HPP
#define SZ3_BITSHUFFLEENCODER_HPP

#include "Encoder.hpp"
#include "SZ3/def.hpp"
#include "SZ3/utils/ByteUtil.hpp"

namespace SZ3 {

    template<class T>
    class BitshuffleEncoder : public concepts::EncoderInterface<T> {
    public:
        BitshuffleEncoder(uint bits) : bits(bits) {}

        void preprocess_encode(const std::vector<T> &bins, int stateNum) override {}

        size_t encode(const std::vector<T> &bins, uchar *&bytes) override {
            size_t num = bins.size();
            size_t total_bits = num * sizeof(T) * 8;
            size_t byte_size = (total_bits + 7) / 8;
            unsigned char* compressed_data = bytes;
            bytes += byte_size;
            size_t out_offset = 0;

            for (uint i = 0; i < sizeof(T) * 8; i += bits) {
                for (size_t j = 0; j < num; j++) {
                    for (uint k = 0; k < bits; k++) {
                        if (i + k < sizeof(T) * 8) {
                            if constexpr (std::is_same_v<T, float>) {
                                if (((*(reinterpret_cast<const uint32_t*>(&bins[j]))) >> (i + k)) & 1) {
                                    compressed_data[out_offset / 8] |= (1 << (out_offset % 8));
                                }
                            } else {
                                if ((bins[j] >> (i + k)) & 1) {
                                    compressed_data[out_offset / 8] |= (1 << (out_offset % 8));
                                }
                            }
                            out_offset++;
                        }
                    }
                }
            }
            return byte_size;
        }

        std::vector<T> decode(const uchar *&bytes, size_t targetLength) override {
            std::vector<T> data(targetLength);
            const unsigned char* compressed_data = bytes;
            size_t in_offset = 0;

            for (uint i = 0; i < sizeof(T) * 8; i += bits) {
                for (size_t j = 0; j < targetLength; j++) {
                    for (uint k = 0; k < bits; k++) {
                        if (i + k < sizeof(T) * 8) {
                            if constexpr (std::is_same_v<T, float>) {
                                if ((compressed_data[in_offset / 8] >> (in_offset % 8)) & 1) {
                                    *(reinterpret_cast<uint32_t*>(&data[j])) |= (1 << (i + k));
                                }
                            } else {
                                if ((compressed_data[in_offset / 8] >> (in_offset % 8)) & 1) {
                                    data[j] |= (1 << (i + k));
                                }
                            }
                            in_offset++;
                        }
                    }
                }
            }
            bytes += (targetLength * sizeof(T) * 8 + 7) / 8;
            return data;
        }

        void save(uchar *&c) override {}

        void load(const uchar *&c, size_t &remaining_length) override {}

        void postprocess_decode() override {}

        void postprocess_encode() override {}

        void preprocess_decode() override {}

    private:
        uint bits;
    };
}
#endif //SZ3_BITSHUFFLEENCODER_HPP
