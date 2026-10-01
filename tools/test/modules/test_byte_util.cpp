#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

#include "SZ3/utils/ByteUtil.hpp"
#include "gtest/gtest.h"

TEST(ByteUtilTest, Zigzag) {
    const std::vector<int64_t> values = {0, -1, 1, -2, 2, INT64_MAX, INT64_MIN};
    const std::vector<uint64_t> expected = {0, 1, 2, 3, 4, UINT64_MAX - 1, UINT64_MAX};
    for (size_t i = 0; i < values.size(); i++) {
        EXPECT_EQ(SZ3::zigzag(values[i]), expected[i]);
        EXPECT_EQ(SZ3::unzigzag(expected[i]), values[i]);
    }
}

TEST(ByteUtilTest, VarintRoundTrip) {
    const std::vector<uint64_t> values = {0, 1, 127, 128, 16383, 16384, uint64_t(1) << 35, UINT64_MAX};
    std::vector<SZ3::uchar> buf(10 * values.size());
    SZ3::uchar *w = buf.data();
    for (uint64_t v : values) SZ3::write_varint(v, w);
    EXPECT_EQ(w - buf.data(), 1 + 1 + 1 + 2 + 2 + 3 + 6 + 10);
    const SZ3::uchar *r = buf.data();
    size_t remaining = size_t(w - buf.data());
    for (uint64_t v : values) EXPECT_EQ(SZ3::read_varint(r, remaining), v);
    EXPECT_EQ(remaining, 0u);
}

TEST(ByteUtilTest, VarintRejectsTruncatedAndOverlong) {
    const SZ3::uchar truncated[] = {0x80, 0x80};
    const SZ3::uchar *r = truncated;
    size_t remaining = sizeof(truncated);
    EXPECT_THROW(SZ3::read_varint(r, remaining), std::out_of_range);

    std::vector<SZ3::uchar> overlong(11, 0x80);
    overlong.back() = 0;
    r = overlong.data();
    remaining = overlong.size();
    EXPECT_THROW(SZ3::read_varint(r, remaining), std::out_of_range);

    std::vector<SZ3::uchar> past64(10, 0x80);  // a 10th byte above 1 holds bits past 64
    past64.back() = 0x02;
    r = past64.data();
    remaining = past64.size();
    EXPECT_THROW(SZ3::read_varint(r, remaining), std::out_of_range);
}

TEST(ByteUtilTest, BitsRoundTrip) {
    std::mt19937_64 gen(7);
    std::vector<std::pair<uint64_t, int>> fields;
    std::vector<SZ3::uchar> bytes;
    SZ3::BitAppender appender(bytes);
    size_t total_bits = 0;
    for (int i = 0; i < 1000; i++) {
        const int bits = int(gen() % 57);
        const uint64_t v = bits == 0 ? 0 : gen() >> (64 - bits);
        fields.emplace_back(v, bits);
        appender.put(v, bits);
        total_bits += bits;
    }
    appender.flush();
    EXPECT_EQ(bytes.size(), (total_bits + 7) / 8);

    SZ3::BitConsumer consumer(bytes.data(), bytes.data() + bytes.size());
    for (const auto &f : fields) EXPECT_EQ(consumer.get(f.second), f.first);
    EXPECT_TRUE(consumer.at_end());
    EXPECT_THROW(consumer.get(8), std::out_of_range);
}
