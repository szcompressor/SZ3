#ifndef SZ3_MULTI_STREAM_DECOMPOSITION_HPP
#define SZ3_MULTI_STREAM_DECOMPOSITION_HPP

#include <array>
#include <cstddef>
#include <vector>

#include "SZ3/def.hpp"

namespace SZ3::concepts {

/**
 * A decomposition that turns the data into NumStreams streams of integer symbols, for SZMultiStreamCompressor, which
 * gives each stream its own encoder. What the streams do not hold (a header, raw bits) goes through save() and load().
 * @tparam T input data type
 * @tparam NumStreams number of symbol streams
 */
template <class T, size_t NumStreams>
class MultiStreamDecompositionInterface {
   public:
    using Streams = std::array<std::vector<int>, NumStreams>;

    virtual ~MultiStreamDecompositionInterface() = default;

    /// The symbol streams of the data.
    virtual Streams compress(const T *data) = 0;

    /// The data from its streams, after load().
    virtual T *decompress(const Streams &streams, T *dec_data) = 0;

    /// What compress() kept besides the streams.
    virtual void save(uchar *&c) = 0;

    /// What save() wrote.
    virtual void load(const uchar *&c, size_t &remaining_length) = 0;

    /// A bound on what save() writes.
    virtual size_t size_est() const = 0;

    /// After load(): the most symbols stream s can hold, which bounds what a decoder allocates for a corrupt stream.
    virtual size_t max_stream_size(size_t s) const = 0;
};

}  // namespace SZ3::concepts

#endif
