#include <stdexcept>
#include <vector>

#include "SZ3/compressor/SZMultiStreamCompressor.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "gtest/gtest.h"

namespace {

/// Two streams: the values at even and at odd positions; the header is the number of values.
class EvenOddDecomposition : public SZ3::concepts::MultiStreamDecompositionInterface<int, 2> {
   public:
    explicit EvenOddDecomposition(size_t n) : n_(n) {}
    Streams compress(const int *data) override {
        Streams streams;
        for (size_t i = 0; i < n_; i++) streams[i % 2].push_back(data[i]);
        return streams;
    }
    int *decompress(const Streams &streams, int *dec_data) override {
        if (streams[0].size() != max_stream_size(0) || streams[1].size() != max_stream_size(1))
            throw std::invalid_argument("stream sizes do not match");
        for (size_t i = 0; i < n_; i++) dec_data[i] = streams[i % 2][i / 2];
        return dec_data;
    }
    void save(SZ3::uchar *&c) override { SZ3::write_varint(n_, c); }
    void load(const SZ3::uchar *&c, size_t &remaining_length) override { n_ = SZ3::read_varint(c, remaining_length); }
    size_t size_est() const override { return 10; }
    size_t max_stream_size(size_t s) const override { return (n_ + 1 - s) / 2; }

   private:
    size_t n_;
};

using EvenOddCompressor = SZ3::SZMultiStreamCompressor<int, EvenOddDecomposition, SZ3::HuffmanEncoder<int>>;

std::vector<SZ3::uchar> compress(std::vector<int> data, size_t capacity) {
    SZ3::Config conf(1);  // the compressor takes its sizes from the decomposition
    EvenOddCompressor compressor{EvenOddDecomposition(data.size())};
    std::vector<SZ3::uchar> out(capacity);
    out.resize(compressor.compress(conf, data.data(), out.data(), out.size()));
    return out;
}

std::vector<int> decompress(const std::vector<SZ3::uchar> &cmp, size_t n) {
    SZ3::Config conf(1);
    EvenOddCompressor compressor{EvenOddDecomposition(0)};
    std::vector<int> dec(n);
    compressor.decompress(conf, cmp.data(), cmp.size(), dec.data());
    return dec;
}

}  // namespace

TEST(MultiStreamCompressorTest, RoundTrip) {
    for (size_t n : {size_t(0), size_t(1), size_t(2), size_t(1001)}) {
        std::vector<int> data(n);
        for (size_t i = 0; i < n; i++) data[i] = i % 2 ? int(i % 7) - 3 : int(i % 100);
        EXPECT_EQ(decompress(compress(data, 1 << 20), n), data);
        // a capacity below the bound goes through the spill buffer and gives the same bytes
        const std::vector<SZ3::uchar> direct = compress(data, 1 << 20);
        EXPECT_EQ(compress(data, direct.size()), direct);
    }
}

TEST(MultiStreamCompressorTest, RejectsSmallBufferAndCorruptStreams) {
    std::vector<int> data(1000, 5);
    const std::vector<SZ3::uchar> cmp = compress(data, 1 << 20);
    EXPECT_THROW(compress(data, cmp.size() - 1), std::length_error);
    // truncated
    EXPECT_ANY_THROW(decompress(std::vector<SZ3::uchar>(cmp.begin(), cmp.begin() + cmp.size() / 2), data.size()));
    // a stream length beyond what the header allows
    std::vector<SZ3::uchar> corrupt = cmp;
    corrupt[0] = 4;  // header: 4 values, so streams of at most 2
    EXPECT_THROW(decompress(corrupt, data.size()), std::out_of_range);
}
