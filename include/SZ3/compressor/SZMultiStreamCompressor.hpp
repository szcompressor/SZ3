#ifndef SZ3_MULTI_STREAM_COMPRESSOR_HPP
#define SZ3_MULTI_STREAM_COMPRESSOR_HPP

#include <cstring>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

#include "SZ3/compressor/Compressor.hpp"
#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {

/**
 * The workflow of a decomposition that turns the data into several streams of integer symbols, each coded with its own
 * Encoder, without a lossless stage. The decomposition provides:
 *  - Streams, an std::array of std::vector<int>;
 *  - Streams compress(const T *data), and T *decompress(const Streams &, T *out) after load();
 *  - save(uchar *&), load(const uchar *&, size_t &) and size_est() for the rest of what it needs;
 *  - max_stream_size(s) after load(): the most symbols stream s can hold, which bounds what the decoder allocates.
 * The Encoder's save() and encode() of n symbols take at most size_est() + 4 n + 40 bytes (HuffmanEncoder: codes of
 * at most 32 bits, four 8-byte part headers, each part padded to a byte).
 *
 * The output: the decomposition's save(), then for each stream its length as a LEB128 varint followed, if the stream
 * is not empty, by its encoder's save() and encode().
 */
template <class T, class Decomposition, class Encoder>
class SZMultiStreamCompressor : public concepts::CompressorInterface<T> {
    static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                  "must implement the encoder interface");
    using Streams = typename Decomposition::Streams;
    static constexpr int NUM_STREAMS = int(std::tuple_size<Streams>::value);

   public:
    explicit SZMultiStreamCompressor(Decomposition decomposition) : decomposition_(std::move(decomposition)) {}

    size_t compress(const Config & /*conf*/, T *data, uchar *cmpData, size_t cmpCap) override {
        const Streams streams = decomposition_.compress(data);
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
            uint64_t n = streams[s].size();  // LEB128: 7 bits a byte, low first, high bit set if more follow
            for (; n >= 128; n >>= 7) *p++ = uchar(n | 128);
            *p++ = uchar(n);
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
        Streams streams;
        for (int s = 0; s < NUM_STREAMS; s++) {
            uint64_t n = 0;
            uint8_t byte = 128;
            for (unsigned shift = 0; byte >= 128; shift += 7) {
                if (shift > 63)
                    throw std::out_of_range("SZ3 multi-stream: corrupt stream, a stream length beyond 64 bits");
                read(byte, p, remaining);
                n |= uint64_t(byte & 127) << shift;
            }
            if (n > decomposition_.max_stream_size(s))
                throw std::out_of_range("SZ3 multi-stream: corrupt stream, a stream longer than the data can hold");
            if (n == 0) continue;
            Encoder encoder;
            encoder.load(p, remaining);
            streams[s] = encoder.decode(p, size_t(n), remaining);
            encoder.postprocess_decode();
        }
        return decomposition_.decompress(streams, decData);
    }

   private:
    Decomposition decomposition_;
};

}  // namespace SZ3
#endif
