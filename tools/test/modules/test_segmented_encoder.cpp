// SegmentedEncoder: bins that are a sequence of segments [n, n bins], each with a Huffman code of its own.

#include <cstdint>
#include <cstring>
#include <random>
#include <stdexcept>
#include <vector>

#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/encoder/SegmentedEncoder.hpp"
#include "gtest/gtest.h"

namespace {

std::vector<int> segmented_round_trip(const std::vector<int> &bins) {
    SZ3::SegmentedEncoder<SZ3::HuffmanEncoder<int>> enc;
    enc.preprocess_encode(bins, 0);
    std::vector<SZ3::uchar> buf(enc.size_est() + 8 * bins.size() + 1024);
    SZ3::uchar *p = buf.data();
    enc.save(p);
    enc.encode(bins, p);
    const size_t n = p - buf.data();
    SZ3::SegmentedEncoder<SZ3::HuffmanEncoder<int>> dec;
    const SZ3::uchar *q = buf.data();
    size_t rem = n;
    dec.load(q, rem);
    return dec.decode(q, bins.size(), rem);
}

TEST(SZ3_SegmentedEncoder, RoundTrip) {
    std::mt19937 g(5);
    std::vector<int> bins;
    for (int s = 0; s < 30; s++) {
        const int n = s % 3 == 0 ? 0 : int(g() % 6000);  // empty segments too
        bins.push_back(n);
        for (int i = 0; i < n; i++) bins.push_back(int(g() % (s + 2)) - s / 2);
    }
    EXPECT_EQ(segmented_round_trip(bins), bins);
    EXPECT_EQ(segmented_round_trip({}), std::vector<int>());
    EXPECT_EQ(segmented_round_trip({0, 0, 0}), std::vector<int>({0, 0, 0}));
}

// Corrupt segment tables: sizes beyond an int, sizes that do not add up (even with wrap-around), and many empty
// segments, which must cost their bytes rather than an encoder each.
TEST(SZ3_SegmentedEncoder, CorruptSegmentTables) {
    auto stream = [](uint32_t count, std::vector<SZ3::uchar> rest) {
        std::vector<SZ3::uchar> b(4);
        memcpy(b.data(), &count, 4);
        b.insert(b.end(), rest.begin(), rest.end());
        return b;
    };
    auto load = [](const std::vector<SZ3::uchar> &b, SZ3::SegmentedEncoder<SZ3::HuffmanEncoder<int>> &dec) {
        const SZ3::uchar *p = b.data();
        size_t rem = b.size();
        dec.load(p, rem);
        return rem;
    };
    {
        SZ3::SegmentedEncoder<SZ3::HuffmanEncoder<int>> dec;
        EXPECT_THROW(load(stream(1, {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x01}), dec),
                     std::out_of_range);                                                        // 2^64 - 1
        EXPECT_THROW(load(stream(1, {0x80, 0x80, 0x80, 0x80, 0x08}), dec), std::out_of_range);  // 2^31
    }
    {
        SZ3::SegmentedEncoder<SZ3::HuffmanEncoder<int>> dec;
        load(stream(3, {0, 0, 0}), dec);
        const SZ3::uchar *p = nullptr;
        size_t rem = 0;
        EXPECT_THROW(dec.decode(p, 2, rem), std::out_of_range);  // three counts are three bins
        EXPECT_EQ(dec.decode(p, 3, rem), std::vector<int>({0, 0, 0}));
    }
    {
        const uint32_t count = 1000000;
        SZ3::SegmentedEncoder<SZ3::HuffmanEncoder<int>> dec;
        load(stream(count, std::vector<SZ3::uchar>(count, 0)), dec);
        const SZ3::uchar *p = nullptr;
        size_t rem = 0;
        EXPECT_EQ(dec.decode(p, count, rem).size(), count);
    }
    {
        // more segments than the caller allows: refused before an encoder is built for each
        SZ3::SegmentedEncoder<SZ3::HuffmanEncoder<int>> dec(12);
        EXPECT_THROW(load(stream(13, std::vector<SZ3::uchar>(13 * 5, 0)), dec), std::out_of_range);
        EXPECT_NO_THROW(load(stream(12, std::vector<SZ3::uchar>(12, 0)), dec));
    }
}

}  // namespace
