#ifndef SZ3_BIOMD_COMPRESSOR_HPP
#define SZ3_BIOMD_COMPRESSOR_HPP

#include <cstring>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "SZ3/compressor/Compressor.hpp"
#include "SZ3/decomposition/SZBioMDDecomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {

/**
 * Specialized compressor for molecular-dynamics coordinates (ALGO_BIOMD): SZBioMDDecomposition turns a chunk into
 * biomd::NUM_STREAMS streams of integer symbols, and each stream is coded with its own Encoder, without a lossless
 * stage. The Encoder's save() and encode() of n symbols take at most size_est() + 4 n + 40 bytes (HuffmanEncoder: codes
 * of at most 32 bits, four 8-byte part headers, each part padded to a byte).
 *
 * The output: the decomposition's save(), then for each stream its length as a LEB128 varint followed, if the stream
 * is not empty, by its encoder's save() and encode().
 */
template <class T, uint N, class Encoder>
class SZBioMDCompressor : public concepts::CompressorInterface<T> {
    static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                  "must implement the encoder interface");
    static constexpr int NUM_STREAMS = biomd::NUM_STREAMS;

   public:
    explicit SZBioMDCompressor(const Config &conf) : decomposition_(conf) {}

    size_t compress(const Config & /*conf*/, T *data, uchar *cmpData, size_t cmpCap) override {
        const biomd::Streams streams = decomposition_.compress(data);
        Encoder encoders[NUM_STREAMS];
        // the decomposition's header, and per stream a varint of at most 10 bytes and what its encoder can write
        size_t bound = decomposition_.size_est();
        for (int s = 0; s < NUM_STREAMS; s++) {
            bound += 10;
            if (streams[s].empty()) continue;
            encoders[s].preprocess_encode(streams[s], 0);
            bound += encoders[s].size_est() + 4 * streams[s].size() + 40;
        }
        // straight into cmpData if it holds the bound; else into a buffer of the bound, not cleared, copied if it fits
        std::unique_ptr<uchar[]> spill(bound > cmpCap ? new uchar[bound] : nullptr);
        uchar *const out = spill ? spill.get() : cmpData;
        uchar *p = out;
        decomposition_.save(p);
        for (int s = 0; s < NUM_STREAMS; s++) {
            write_varint(streams[s].size(), p);
            if (streams[s].empty()) continue;
            encoders[s].save(p);
            encoders[s].encode(streams[s], p);
            encoders[s].postprocess_encode();
        }
        const size_t size = size_t(p - out);
        if (out != cmpData) {
            if (size > cmpCap) throw std::length_error(SZ3_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH);
            memcpy(cmpData, out, size);
        }
        return size;
    }

    T *decompress(const Config & /*conf*/, uchar const *cmpData, size_t cmpSize, T *decData) override {
        const uchar *p = cmpData;
        size_t remaining = cmpSize;
        decomposition_.load(p, remaining);
        biomd::Streams streams;
        for (int s = 0; s < NUM_STREAMS; s++) {
            const uint64_t n = read_varint(p, remaining);
            if (n > decomposition_.max_stream_size(s))
                throw std::out_of_range("SZ3 BioMD: corrupt stream, a stream longer than the data can hold");
            if (n == 0) continue;
            Encoder encoder;
            encoder.load(p, remaining);
            streams[s] = encoder.decode(p, size_t(n), remaining);
            encoder.postprocess_decode();
        }
        return decomposition_.decompress(streams, decData);
    }

   private:
    SZBioMDDecomposition<T, N> decomposition_;
};

}  // namespace SZ3
#endif
