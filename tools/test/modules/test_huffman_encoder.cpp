// HuffmanEncoder: exact round trips, save() == size_est(), size_bound() covers save() + encode(),
// optimal payloads, stable bytes, and decoders that throw rather than read out of bounds on corrupted streams.
//
// Every buffer handed to save(), encode(), load() and decode() is an exactly sized heap allocation, so a sanitizer
// build reports an access one byte past what the stream holds.

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <map>
#include <memory>
#include <queue>
#include <random>
#include <stdexcept>
#include <vector>

#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/encoder/SegmentedEncoder.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "gtest/gtest.h"

namespace {

using SZ3::uchar;

template <class T>
struct Encoded {
    std::vector<uchar> tree, data;
    size_t est = 0;
};

template <class T>
Encoded<T> encode(const std::vector<T> &bins) {
    Encoded<T> e;
    SZ3::HuffmanEncoder<T> enc;
    enc.preprocess_encode(bins, 0);
    e.est = enc.size_est();
    {
        std::unique_ptr<uchar[]> t(new uchar[e.est]);
        uchar *p = t.get();
        enc.save(p);
        e.tree.assign(t.get(), p);
    }
    // First into a roomy buffer to learn the size, then again into an exact one.
    size_t size;
    {
        std::vector<uchar> big(32 + bins.size() * 5);
        uchar *p = big.data();
        size = enc.encode(bins, p);
        EXPECT_EQ(size, static_cast<size_t>(p - big.data()));
    }
    std::unique_ptr<uchar[]> d(new uchar[size]);
    uchar *p = d.get();
    EXPECT_EQ(enc.encode(bins, p), size);
    e.data.assign(d.get(), p);
    enc.postprocess_encode();
    return e;
}

template <class T>
std::vector<T> decode(const Encoded<T> &e, size_t n) {
    std::unique_ptr<uchar[]> t(new uchar[e.tree.size()]);
    std::unique_ptr<uchar[]> d(new uchar[e.data.size()]);
    if (!e.tree.empty()) memcpy(t.get(), e.tree.data(), e.tree.size());
    if (!e.data.empty()) memcpy(d.get(), e.data.data(), e.data.size());
    SZ3::HuffmanEncoder<T> dec;
    const uchar *tp = t.get();
    size_t trem = e.tree.size();
    dec.load(tp, trem);
    EXPECT_EQ(trem, 0u);
    dec.preprocess_decode();
    const uchar *dp = d.get();
    size_t drem = e.data.size();
    auto out = dec.decode(dp, n, drem);
    dec.postprocess_decode();
    EXPECT_EQ(drem, 0u);
    return out;
}

// The parts encode() writes, each with its bit count in front: one below 4096 bins, else four, unless the top bit of
// the first count marks a single part.
size_t payload_parts(size_t n, size_t distinct) { return n >= (size_t(1) << 12) && distinct >= 2 ? 4 : 1; }
uint64_t payload_bits(const std::vector<uchar> &v, size_t parts = 1, uint64_t *bytes = nullptr) {
    uint64_t bits = 0, b = 0;
    const uchar *p = v.data();
    for (size_t k = 0; k < parts; k++) {
        SZ3::read(b, p);
        if (k == 0 && parts == 4 && (b >> 63)) {
            b &= ~(uint64_t(1) << 63);
            parts = 1;
        }
        bits += b;
        if (bytes) *bytes += (b + 7) / 8;
    }
    if (bytes) *bytes += parts * sizeof(uint64_t);
    return bits;
}

// Cost in bits of an optimal prefix code for these bins, and the depth of one optimal tree.
template <class T>
uint64_t optimal_bits(const std::vector<T> &bins, int *depth = nullptr) {
    std::map<T, uint64_t> freq;
    for (T b : bins) freq[b]++;
    if (depth) *depth = 0;
    if (freq.size() < 2) return 0;
    std::priority_queue<std::pair<uint64_t, int>, std::vector<std::pair<uint64_t, int>>, std::greater<>> q;
    for (auto &kv : freq) q.push({kv.second, 0});
    uint64_t cost = 0;
    int d = 0;
    while (q.size() > 1) {
        auto a = q.top();
        q.pop();
        auto b = q.top();
        q.pop();
        cost += a.first + b.first;
        d = std::max(a.second, b.second) + 1;
        q.push({a.first + b.first, d});
    }
    if (depth) *depth = d;
    return cost;
}

template <class T>
size_t distinct(std::vector<T> v) {
    std::sort(v.begin(), v.end());
    return std::unique(v.begin(), v.end()) - v.begin();
}

// Everything a valid stream must satisfy.
template <class T>
Encoded<T> check(const std::vector<T> &bins, bool expect_optimal = true) {
    Encoded<T> e = encode(bins);
    EXPECT_EQ(e.tree.size(), e.est) << "save() must write exactly size_est()";
    EXPECT_LE(e.tree.size() + e.data.size(), SZ3::HuffmanEncoder<T>::size_bound(bins.size(), distinct(bins)));
    uint64_t bytes = 0;
    const uint64_t bits = payload_bits(e.data, payload_parts(bins.size(), distinct(bins)), &bytes);
    EXPECT_EQ(e.data.size(), bytes);
    if (expect_optimal) EXPECT_EQ(bits, optimal_bits(bins));
    auto out = decode(e, bins.size());
    EXPECT_TRUE(out == bins) << "round trip differs, n = " << bins.size();
    return e;
}

std::vector<int> laplace(size_t n, double scale, int center, uint64_t seed) {
    std::mt19937_64 g(seed);
    std::exponential_distribution<double> ex(1.0 / scale);
    std::vector<int> v(n);
    for (size_t i = 0; i < n; i++) {
        double x = ex(g);
        v[i] = center + static_cast<int>(std::lround((g() & 1) ? x : -x));
    }
    return v;
}

TEST(SZ3_HuffmanEncoder, Empty) {
    std::vector<int> bins;
    auto e = check(bins);
    EXPECT_EQ(e.tree.size(), sizeof(uint32_t));
    EXPECT_EQ(e.data.size(), sizeof(uint64_t));
    // An empty table cannot produce values.
    SZ3::HuffmanEncoder<int> dec;
    const uchar *tp = e.tree.data();
    size_t trem = e.tree.size();
    dec.load(tp, trem);
    const uchar *dp = e.data.data();
    size_t drem = e.data.size();
    EXPECT_THROW(dec.decode(dp, 1, drem), std::out_of_range);
}

TEST(SZ3_HuffmanEncoder, SingleSymbol) {
    for (int v : {0, 1, -1, 32768, INT_MAX, INT_MIN}) {
        for (size_t n : {1, 2, 1000}) {
            auto e = check(std::vector<int>(n, v));
            EXPECT_EQ(e.data.size(), sizeof(uint64_t));
        }
    }
}

TEST(SZ3_HuffmanEncoder, SmallAlphabetsEveryBitBoundary) {
    std::mt19937_64 g(1);
    for (int k = 2; k <= 9; k++) {
        for (size_t n = 1; n <= 130; n++) {
            std::vector<int> bins(n);
            for (auto &b : bins) b = 100 + static_cast<int>(g() % k);
            check(bins);
        }
    }
}

TEST(SZ3_HuffmanEncoder, QuantizationLike) {
    for (double scale : {0.3, 1.0, 5.0, 60.0, 3000.0}) {
        for (size_t n : {10, 1000, 7956, 100000}) {
            auto bins = laplace(n, scale, 32768, static_cast<uint64_t>(scale * 1000) + n);
            for (size_t i = 0; i < n; i += 97) bins[i] = 0;  // unpredictable
            check(bins);
        }
    }
}

// A dominant symbol gets the one-bit all-zero code, and the decoder takes its runs in one step.
TEST(SZ3_HuffmanEncoder, LongRunsOfTheCommonestSymbol) {
    std::mt19937_64 g(16);
    for (int rarity : {3, 20, 200, 5000}) {
        for (size_t n : {2, 9, 63, 64, 65, 200, 4099, 60000}) {
            std::vector<int> bins(n, 32768);
            for (auto &b : bins)
                if (g() % rarity == 0) b = 32768 + static_cast<int>(g() % 9) - 4;
            check(bins);
            bins.assign(n, 7);
            bins.back() = 8;  // a run up to the last code
            check(bins);
            bins.assign(n, 7);
            bins.front() = -8;
            check(bins);
        }
    }
}

// Table-length codes and longer ones, around every length the decoder switches on.
TEST(SZ3_HuffmanEncoder, CodeLengthsAroundTheTable) {
    for (int top = 8; top <= 22; top++) {
        std::vector<int> bins;
        // Geometric weights 2^-k give code lengths 1..top.
        for (int k = 0; k <= top; k++)
            for (int j = 0; j < (1 << (top - k)); j++) bins.push_back(k * 3 - 20);
        std::shuffle(bins.begin(), bins.end(), std::mt19937_64(top));
        check(bins);
    }
}

TEST(SZ3_HuffmanEncoder, NegativeAndExtremeBins) {
    check(std::vector<int>{INT_MIN, INT_MAX, INT_MIN, 0, -1, 1, INT_MAX});
    check(laplace(5000, 20, -100000, 3));
    std::vector<int64_t> w{INT64_MIN, INT64_MAX, 0, INT64_MIN, -5, 5};
    check(w);
    std::vector<uint64_t> u{0, UINT64_MAX, UINT64_MAX, 12345};
    check(u);
}

TEST(SZ3_HuffmanEncoder, OtherTypes) {
    std::mt19937_64 g(4);
    std::vector<uint8_t> a(3000);
    for (auto &x : a) x = static_cast<uint8_t>(g() % 7 + 250 * (g() % 50 == 0));
    check(a);
    std::vector<int16_t> b(3000);
    for (auto &x : b) x = static_cast<int16_t>(static_cast<int>(g() % 11) - 5 + (g() % 100 == 0 ? -32768 : 0));
    check(b);
    std::vector<uint16_t> c(3000);
    for (auto &x : c) x = static_cast<uint16_t>(g() % 3 ? 65535 : g() % 40);
    check(c);
}

// A range far wider than the input takes the hash path.
TEST(SZ3_HuffmanEncoder, SparseWideRange) {
    std::mt19937_64 g(5);
    std::vector<int> vals(3000);
    for (auto &x : vals) x = static_cast<int>(g());
    std::geometric_distribution<int> geo(0.003);
    for (size_t n : {2, 50, 4000, 200000}) {
        std::vector<int> bins(n);
        for (auto &b : bins) b = vals[std::min(geo(g), 2999)];
        check(bins);
    }
    // All distinct.
    std::vector<int> all(20000);
    for (size_t i = 0; i < all.size(); i++) all[i] = static_cast<int>(i * 104729u) - 1000000000;
    check(all);
}

// Fibonacci weights: an optimal code is 34 bits deep, so the lengths are limited to kMaxLen.
TEST(SZ3_HuffmanEncoder, LengthLimited) {
    std::vector<int> bins;
    uint64_t a = 1, b = 1;
    for (int k = 0; k < 35; k++) {
        bins.insert(bins.end(), a, k);
        uint64_t c = a + b;
        a = b;
        b = c;
    }
    std::shuffle(bins.begin(), bins.end(), std::mt19937_64(6));
    int depth;
    const uint64_t opt = optimal_bits(bins, &depth);
    ASSERT_GT(depth, 32);
    auto e = check(bins, false);
    const uint64_t bits = payload_bits(e.data, payload_parts(bins.size(), distinct(bins)));
    EXPECT_GE(bits, opt);
    EXPECT_LE(bits, opt + opt / 1000000);
}

// The same bins give the same bytes on every host; these were produced by this implementation.
TEST(SZ3_HuffmanEncoder, StableBytes) {
    std::vector<int> bins{5, 5, 5, 7, 7, 9, 5, 12, 5, 7};
    auto e = check(bins);
    // Lengths 5:1 7:2 9:3 12:3, so codes 0, 10, 110, 111 and 17 payload bits. After the count 4 and the offset 5
    // (little-endian), the table holds gaps 2 2 3 and length changes +1 +1 +1 0.
    const std::vector<uchar> tree{4, 0, 0, 0, 5, 0, 0, 0, 0x69, 0xa6, 0xe0};
    const std::vector<uchar> data{17, 0, 0, 0, 0, 0, 0, 0, 0x15, 0x9d, 0x00};
    EXPECT_EQ(e.tree, tree);
    EXPECT_EQ(e.data, data);
    // Independent of what the encoder was used for before.
    SZ3::HuffmanEncoder<int> enc;
    enc.preprocess_encode(laplace(1000, 3, 0, 7), 0);
    enc.preprocess_encode(bins, 0);
    std::vector<uchar> t(enc.size_est());
    uchar *p = t.data();
    enc.save(p);
    EXPECT_EQ(t, tree);
}

TEST(SZ3_HuffmanEncoder, CopiesAreIndependent) {
    auto bins = laplace(3000, 4, 10, 8);
    SZ3::HuffmanEncoder<int> enc;
    enc.preprocess_encode(bins, 0);
    auto copy = enc;
    enc.preprocess_encode(std::vector<int>{1, 2, 3}, 0);
    std::vector<uchar> buf(copy.size_est() + 16 + bins.size() * 4);
    uchar *p = buf.data();
    copy.save(p);
    copy.encode(bins, p);
    SZ3::HuffmanEncoder<int> dec;
    const uchar *c = buf.data();
    size_t rem = p - buf.data();
    dec.load(c, rem);
    EXPECT_EQ(dec.decode(c, bins.size(), rem), bins);
}

TEST(SZ3_HuffmanEncoder, RefusesBinsItWasNotBuiltFor) {
    SZ3::HuffmanEncoder<int> enc;
    enc.preprocess_encode(std::vector<int>{1, 3, 3}, 0);
    std::vector<uchar> buf(64);
    uchar *p = buf.data();
    EXPECT_THROW(enc.encode(std::vector<int>{1, 3, 9}, p), std::invalid_argument);  // outside the range
    p = buf.data();
    EXPECT_THROW(enc.encode(std::vector<int>{1, 2, 3}, p), std::invalid_argument);  // inside, never counted
    enc.preprocess_encode(std::vector<int>{4, 4}, 0);
    p = buf.data();
    EXPECT_THROW(enc.encode(std::vector<int>{4, 5}, p), std::invalid_argument);  // one symbol
    enc.preprocess_encode(std::vector<int>{}, 0);
    p = buf.data();
    EXPECT_THROW(enc.encode(std::vector<int>{0}, p), std::invalid_argument);  // no symbol
    enc.preprocess_encode(std::vector<int>{0, 1000000, 1000000}, 0);
    p = buf.data();
    EXPECT_THROW(enc.encode(std::vector<int>{0, 1000000, 7}, p), std::invalid_argument);  // hashed
}

TEST(SZ3_HuffmanEncoder, WrongValueCountThrows) {
    auto bins = laplace(5000, 6, 0, 9);
    auto e = encode(bins);
    for (size_t n : {size_t(0), size_t(1), bins.size() - 1, bins.size() + 1, bins.size() * 40}) {
        SZ3::HuffmanEncoder<int> dec;
        const uchar *tp = e.tree.data();
        size_t trem = e.tree.size();
        dec.load(tp, trem);
        const uchar *dp = e.data.data();
        size_t drem = e.data.size();
        EXPECT_THROW(dec.decode(dp, n, drem), std::out_of_range) << n;
    }
}

// Decode a possibly corrupted stream from exactly sized buffers; true if it decoded.
bool try_decode(const std::vector<uchar> &tree, const std::vector<uchar> &data, size_t n) {
    std::unique_ptr<uchar[]> te(new uchar[tree.size()]);
    std::unique_ptr<uchar[]> de(new uchar[data.size()]);
    if (!tree.empty()) memcpy(te.get(), tree.data(), tree.size());
    if (!data.empty()) memcpy(de.get(), data.data(), data.size());
    try {
        SZ3::HuffmanEncoder<int> dec;
        const uchar *tp = te.get();
        size_t trem = tree.size();
        dec.load(tp, trem);
        if (tp > te.get() + tree.size() || trem > tree.size()) ADD_FAILURE() << "load overran";
        const uchar *dp = de.get();
        size_t drem = data.size();
        auto out = dec.decode(dp, n, drem);
        if (out.size() != n) ADD_FAILURE() << "decode returned " << out.size() << " values";
        if (dp > de.get() + data.size() || drem > data.size()) ADD_FAILURE() << "decode overran";
        return true;
    } catch (const std::out_of_range &) {
        return false;
    } catch (const std::length_error &) {
        return false;
    } catch (const std::bad_alloc &) {
        ADD_FAILURE() << "allocation driven by the stream";
        return false;
    }
}

TEST(SZ3_HuffmanEncoder, TruncationsThrow) {
    for (auto bins : {laplace(3000, 2, 32768, 10), laplace(3000, 300, 32768, 11)}) {
        auto e = encode(bins);
        for (size_t k = 0; k < e.tree.size(); k++)
            EXPECT_FALSE(try_decode(std::vector<uchar>(e.tree.begin(), e.tree.begin() + k), e.data, bins.size()));
        for (size_t k = 0; k < e.data.size(); k++)
            EXPECT_FALSE(try_decode(e.tree, std::vector<uchar>(e.data.begin(), e.data.begin() + k), bins.size()));
    }
}

TEST(SZ3_HuffmanEncoder, MutationsNeverOverrun) {
    std::mt19937_64 g(12);
    std::vector<std::vector<int>> inputs{laplace(2000, 1, 32768, 13), laplace(2000, 50, 32768, 14),
                                         laplace(300, 3000, 0, 15), std::vector<int>(100, 7)};
    {
        std::vector<int> fib;  // long codes
        uint64_t a = 1, b = 1;
        for (int k = 0; k < 20; k++) {
            fib.insert(fib.end(), a, k * 1000);
            uint64_t c = a + b;
            a = b;
            b = c;
        }
        inputs.push_back(fib);
    }
    size_t decoded = 0, total = 0;
    for (auto &bins : inputs) {
        auto e = encode(bins);
        for (int it = 0; it < 3000; it++) {
            auto t = e.tree;
            auto d = e.data;
            int flips = 1 + static_cast<int>(g() % 4);
            for (int f = 0; f < flips; f++) {
                bool in_tree = (g() & 1) && !t.empty();
                auto &v = in_tree ? t : d;
                if (v.empty()) continue;
                switch (g() % 3) {
                    case 0:
                        v[g() % v.size()] ^= static_cast<uchar>(1u << (g() % 8));
                        break;
                    case 1:
                        v[g() % v.size()] = static_cast<uchar>(g());
                        break;
                    default:
                        if (v.size() > 1) v.resize(g() % v.size());
                        break;
                }
            }
            total++;
            decoded += try_decode(t, d, bins.size());
        }
    }
    // Some mutations still decode (a flipped payload bit can map to another valid code); none may overrun.
    RecordProperty("decoded", static_cast<int>(decoded));
    RecordProperty("total", static_cast<int>(total));
}

// A value count the payload cannot hold is refused before the output is allocated.
TEST(SZ3_HuffmanEncoder, ValueCountBeyondThePayloadIsRefusedFirst) {
    std::vector<int> bins(5000, 0);
    for (size_t i = 0; i < bins.size(); i += 7) bins[i] = int(i % 5);
    const auto e = encode(bins);
    SZ3::HuffmanEncoder<int> dec;
    const SZ3::uchar *p = e.tree.data();
    size_t rem = e.tree.size();
    dec.load(p, rem);
    p = e.data.data();
    rem = e.data.size();
    EXPECT_THROW(dec.decode(p, size_t(1) << 40, rem), std::out_of_range);  // would be 4 TB of ints
}

// A bin preprocess_encode() did not see has no code.
TEST(SZ3_HuffmanEncoder, UnseenBinIsRefused) {
    SZ3::HuffmanEncoder<int> enc;
    enc.preprocess_encode({0, 0, 1, 3}, 0);
    std::vector<SZ3::uchar> buf(1024);
    SZ3::uchar *p = buf.data();
    EXPECT_THROW(enc.encode({1, 1, 1, 2}, p), std::invalid_argument);
}

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
