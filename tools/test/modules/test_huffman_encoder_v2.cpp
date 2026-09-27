// HuffmanEncoderV2: exact round trips, save() within size_est(), payloads of optimal length, stable bytes,
// and the inputs it refuses.
//
// Buffers handed to save() and to the decoder are exactly sized heap allocations, so a sanitizer build reports a write
// past size_est() or a read past what save()/encode() produced at the offending access.

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

#include "SZ3/decomposition/SZBioMDDecomposition.hpp"
#include "SZ3/encoder/HuffmanEncoderV2.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include "gtest/gtest.h"

namespace {

using V2 = SZ3::HuffmanEncoderV2<int>;

struct Encoded {
    std::vector<uint8_t> tree;  // save()
    std::vector<uint8_t> data;  // encode()
    size_t est = 0;             // size_est() after preprocess_encode()
};

Encoded encode_v2(const std::vector<int>& bins, int stateNum, uint8_t flag = 0x00) {
    Encoded e;
    V2 enc;
    enc.preprocess_encode(bins.data(), bins.size(), stateNum, flag);
    e.est = enc.size_est();
    std::unique_ptr<uint8_t[]> t(new uint8_t[e.est]);
    uint8_t* tp = t.get();
    enc.save(tp);
    e.tree.assign(t.get(), tp);

    std::unique_ptr<uint8_t[]> d(new uint8_t[16 + 8 * bins.size()]);
    uint8_t* dp = d.get();
    const size_t ret = enc.encode(bins, dp);
    e.data.assign(d.get(), dp);
    EXPECT_EQ(ret, e.data.size());
    enc.postprocess_encode();
    return e;
}

// Tree and payload in separate, exactly sized buffers.
std::vector<int> decode_separate(const Encoded& e, size_t n) {
    std::unique_ptr<uint8_t[]> t(new uint8_t[e.tree.size()]);
    std::unique_ptr<uint8_t[]> d(new uint8_t[e.data.size()]);
    memcpy(t.get(), e.tree.data(), e.tree.size());
    memcpy(d.get(), e.data.data(), e.data.size());
    V2 dec;
    const uint8_t* tp = t.get();
    size_t trem = e.tree.size();
    dec.load(tp, trem);
    EXPECT_EQ(trem, 0u);
    EXPECT_TRUE(tp == t.get() + e.tree.size());
    dec.preprocess_decode();
    const uint8_t* dp = d.get();
    size_t drem = e.data.size();
    auto out = dec.decode(dp, n, drem);
    dec.postprocess_decode();
    EXPECT_EQ(drem, 0u);
    EXPECT_TRUE(dp == d.get() + e.data.size());
    return out;
}

// Code bits the payload declares in its 8-byte header.
uint64_t payload_bits(const Encoded& e) { return SZ3::bytesToInt64_bigEndian(e.data.data()) ^ 0x1234abcd; }

// Cost of an optimal prefix code for these bins (sum of internal-node weights), and its longest code.
uint64_t optimal_bits(const std::vector<int>& bins, int* max_len = nullptr) {
    std::map<int, uint64_t> freq;
    for (int b : bins) freq[b]++;
    if (freq.size() < 2) {
        if (max_len) *max_len = 0;
        return 0;
    }
    // (weight, depth of the deepest leaf below)
    std::priority_queue<std::pair<uint64_t, int>, std::vector<std::pair<uint64_t, int>>, std::greater<>> q;
    for (auto& kv : freq) q.push({kv.second, 0});
    uint64_t cost = 0;
    int depth = 0;
    while (q.size() > 1) {
        auto a = q.top();
        q.pop();
        auto b = q.top();
        q.pop();
        cost += a.first + b.first;
        depth = std::max(a.second, b.second) + 1;
        q.push({a.first + b.first, depth});
    }
    if (max_len) *max_len = depth;  // one optimal tree's depth; ties may give another optimal depth
    return cost;
}

size_t distinct(const std::vector<int>& bins) {
    std::vector<int> s(bins);
    std::sort(s.begin(), s.end());
    return std::unique(s.begin(), s.end()) - s.begin();
}

// Everything a valid stream must satisfy. Returns the encoding for further checks.
Encoded check_round_trip(const std::vector<int>& bins, int stateNum, uint8_t flag = 0x00) {
    Encoded e = encode_v2(bins, stateNum, flag);
    EXPECT_LE(e.tree.size(), e.est) << "save() wrote more than size_est()";

    const size_t k = distinct(bins);
    const bool raw = (flag & 0x01) != 0;
    if (k >= 2) {
        const uint64_t bits = payload_bits(e);
        if (!raw) {
            EXPECT_EQ(bits, optimal_bits(bins)) << "payload is not an optimal prefix code";
        }
        EXPECT_EQ(e.data.size(), 8 + (bits + 7) / 8);
    } else {
        EXPECT_EQ(e.data.size(), 8u);
    }

    EXPECT_TRUE(decode_separate(e, bins.size()) == bins) << "n=" << bins.size();
    return e;
}

uint64_t fnv1a(const Encoded& e) {
    uint64_t h = 1469598103934665603ull;
    for (auto* v : {&e.tree, &e.data})
        for (uint8_t c : *v) h = (h ^ c) * 1099511628211ull;
    return h;
}

// Symbol s appears F(s+1) times (1, 1, 2, 3, 5, ...), which forces a code of k-1 bits.
std::vector<int> fibonacci_stream(int k, const std::vector<int>& values) {
    std::vector<int> v;
    uint64_t a = 1, b = 1;
    for (int s = 0; s < k; s++) {
        v.insert(v.end(), a, values[s]);
        const uint64_t c = a + b;
        a = b;
        b = c;
    }
    return v;
}

std::vector<int> iota_values(int k, int base) {
    std::vector<int> v(k);
    for (int i = 0; i < k; i++) v[i] = base + i;
    return v;
}

// Random symbol streams, and the stateNum to encode them with.
std::vector<int> random_stream(std::mt19937_64& rng, int dist, size_t n, int& stateNum) {
    std::vector<int> v(n);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    const int kind = static_cast<int>(rng() % 3);  // 0: offset 0 + stateNum, 1: offset + stateNum 0, 2: far offset
    const int base = kind == 0 ? 0 : kind == 1 ? static_cast<int>(rng() % 70000) : 1000000000;
    int k = 2 + static_cast<int>(rng() % 5000);
    switch (dist) {
        case 0:  // uniform over k symbols
            for (auto& x : v) x = static_cast<int>(rng() % k);
            break;
        case 1: {  // two-sided geometric around a centre, like quantization bins
            const double p = 0.02 + 0.9 * u01(rng);
            std::geometric_distribution<int> g(p);
            k = 32768;
            for (auto& x : v) {
                const int m = std::min(g(rng), k / 2 - 1);
                x = k / 2 + ((rng() & 1) ? m : -m);
            }
            break;
        }
        case 2: {  // Zipf over k symbols
            const double s = 0.5 + 2.5 * u01(rng);
            std::vector<double> w(k);
            for (int i = 0; i < k; i++) w[i] = 1.0 / std::pow(i + 1.0, s);
            std::discrete_distribution<int> z(w.begin(), w.end());
            std::vector<int> perm(k);
            for (int i = 0; i < k; i++) perm[i] = i;
            std::shuffle(perm.begin(), perm.end(), rng);
            for (auto& x : v) x = perm[z(rng)];
            break;
        }
        case 3: {  // random histogram, heavy-tailed counts, sparse over a wide range
            k = 2 + static_cast<int>(rng() % 200);
            std::vector<int> syms(k);
            for (auto& s : syms) s = static_cast<int>(rng() % (1u << 20));
            std::vector<double> w(k);
            for (auto& x : w) x = std::exp(12 * u01(rng));
            std::discrete_distribution<int> z(w.begin(), w.end());
            for (auto& x : v) x = syms[z(rng)];
            break;
        }
    }
    int maxv = 0;
    for (auto& x : v) {
        x += base;
        maxv = std::max(maxv, x);
    }
    stateNum = kind == 0 ? (rng() & 1 ? maxv + 1 : std::max(maxv + 1, 65536)) : 0;
    return v;
}

}  // namespace

TEST(SZ3_HuffmanEncoderV2, SingleSymbol) {
    for (int value : {0, 1, 7, 65535, INT_MAX - 1, INT_MAX}) {
        for (size_t n : {1u, 1000u}) {
            std::vector<int> bins(n, value);
            check_round_trip(bins, 0);
            if (value < INT_MAX) check_round_trip(bins, value + 1);
            if (value < 65536) check_round_trip(bins, 65536);
        }
    }
}

TEST(SZ3_HuffmanEncoderV2, TwoSymbolsEveryBitBoundary) {
    // One bit per symbol: the payload is exactly n bits, so n walks every residue mod 8.
    for (size_t n = 2; n <= 17; n++) {
        std::vector<int> bins(n);
        for (size_t i = 0; i < n; i++) bins[i] = (i * 5 + i / 3) % 2 == 0 ? 11 : 4;
        if (distinct(bins) < 2) bins[0] = bins[0] == 11 ? 4 : 11;
        check_round_trip(bins, 0);
        for (auto& b : bins) b = b == 11 ? 1 : 0;
        check_round_trip(bins, 2);
        check_round_trip(bins, 65536);
    }
}

TEST(SZ3_HuffmanEncoderV2, SmallAlphabetsEveryBitBoundary) {
    for (int k : {3, 4, 5, 8, 9, 255, 256, 257}) {
        for (size_t n = static_cast<size_t>(k); n <= static_cast<size_t>(k) + 16; n++) {
            std::vector<int> bins(n);
            for (size_t i = 0; i < n; i++) bins[i] = static_cast<int>((i * 7919) % k);
            check_round_trip(bins, 0);
            check_round_trip(bins, k);
        }
    }
}

TEST(SZ3_HuffmanEncoderV2, CodeLengthsAroundTheTableThreshold) {
    // The decoder walks the tree bit by bit when the longest code is <= 16 bits and switches to a 16-bit
    // lookup table above that. Fibonacci counts give a longest code of exactly k - 1 bits.
    std::mt19937_64 rng(1);
    for (int k = 2; k <= 25; k++) {
        std::vector<int> values = iota_values(k, 0);
        auto bins = fibonacci_stream(k, values);
        check_round_trip(bins, 0);
        check_round_trip(bins, k);
        check_round_trip(bins, 65536);
        // Same shape, symbols scattered over a wide range and the stream shuffled.
        for (auto& v : values) v = static_cast<int>(rng() % 2000001);
        std::sort(values.begin(), values.end());
        values.erase(std::unique(values.begin(), values.end()), values.end());
        if (static_cast<int>(values.size()) < k) continue;
        std::shuffle(values.begin(), values.end(), rng);
        bins = fibonacci_stream(k, values);
        std::shuffle(bins.begin(), bins.end(), rng);
        check_round_trip(bins, 0);
    }
}

TEST(SZ3_HuffmanEncoderV2, DeepCodesEveryPayloadBoundary) {
    // Longest code 17 and 20 bits (table path), payload ending at every residue mod 8.
    for (int k : {18, 21}) {
        auto base = fibonacci_stream(k, iota_values(k, 100));
        for (int extra = 0; extra < 8; extra++) {
            auto bins = base;
            bins.insert(bins.end(), extra, 100 + k - 1);  // more of the 1-bit symbol
            check_round_trip(bins, 0);
            check_round_trip(bins, 100 + k);
        }
    }
}

TEST(SZ3_HuffmanEncoderV2, CodeLongerThan32Bits) {
    // 34 symbols with Fibonacci counts (about 15M values) give a 33-bit code.
    auto bins = fibonacci_stream(34, iota_values(34, 0));
    int depth = 0;
    const uint64_t bits = optimal_bits(bins, &depth);
    ASSERT_EQ(depth, 33);
    Encoded e = encode_v2(bins, 0);
    EXPECT_EQ(payload_bits(e), bits);
    EXPECT_TRUE(decode_separate(e, bins.size()) == bins);
}

TEST(SZ3_HuffmanEncoderV2, TableMixesShortAndLongCodes) {
    // Many leaves at depths 10..16 alongside a chain below 16, so table entries hold both leaves and
    // internal nodes.
    std::vector<int> bins;
    for (int s = 0; s < 3000; s++) bins.insert(bins.end(), 1 + (s % 7), s);
    auto tail = fibonacci_stream(24, iota_values(24, 5000));
    bins.insert(bins.end(), tail.begin(), tail.end());
    check_round_trip(bins, 0);
    check_round_trip(bins, 65536);
}

TEST(SZ3_HuffmanEncoderV2, LargeAndSkewedAlphabets) {
    std::mt19937_64 rng(2);
    // Every bin of the default quantizer (radius 32768), in the dense and the map representation.
    for (size_t n : {50000u, 140000u}) {
        std::vector<int> bins(n);
        for (auto& b : bins) b = static_cast<int>(rng() % 65536);
        check_round_trip(bins, 65536);
        check_round_trip(bins, 0);
    }
    // A radius of 2^20 and 2^27 with a few thousand symbols (map mode), and one span just under 2^31.
    for (int span : {1 << 21, 1 << 28, INT_MAX - 1}) {
        std::vector<int> bins(20000);
        for (auto& b : bins) b = static_cast<int>(rng() % 3000) * (span / 3000);
        bins[0] = 0;
        bins[1] = span - 1;
        check_round_trip(bins, 0);
        if (span <= (1 << 28)) check_round_trip(bins, span);
    }
    // Large magnitudes, small span.
    for (int base : {INT_MAX - 40, 1, 1000000000}) {
        std::vector<int> bins(3000);
        for (auto& b : bins) b = base + static_cast<int>(rng() % 41);
        check_round_trip(bins, 0);
    }
    // 99.99% one symbol.
    std::vector<int> bins(200000, 32768);
    for (int i = 0; i < 20; i++) bins[rng() % bins.size()] = static_cast<int>(rng() % 65536);
    check_round_trip(bins, 65536);
    check_round_trip(bins, 0);
}

TEST(SZ3_HuffmanEncoderV2, RepresentationsAgree) {
    // Auto, forced map (0x40) and forced dense (0x80) choose the same tree, so the payload is identical;
    // the raw fixed-width mode (0x01) round-trips too.
    std::mt19937_64 rng(3);
    for (int trial = 0; trial < 40; trial++) {
        int stateNum = 0;
        auto bins = random_stream(rng, trial % 4, 1 + rng() % 3000, stateNum);
        if (stateNum > (1 << 20)) stateNum = 0;
        int lo = *std::min_element(bins.begin(), bins.end());
        int hi = *std::max_element(bins.begin(), bins.end());
        if (static_cast<int64_t>(hi) - lo > (1 << 22)) continue;
        Encoded a = check_round_trip(bins, stateNum);
        Encoded m = check_round_trip(bins, stateNum, 0x40);
        Encoded d = check_round_trip(bins, stateNum, 0x80);
        EXPECT_TRUE(a.data == m.data && m.data == d.data);
        check_round_trip(bins, stateNum, 0x01);
    }
}

TEST(SZ3_HuffmanEncoderV2, RawModeEveryWidth) {
    for (int width = 1; width <= 30; width++) {
        for (size_t n : {1u, 2u, 3u, 7u, 8u, 9u, 31u, 32u, 33u, 100u}) {
            std::vector<int> bins(n);
            for (size_t i = 0; i < n; i++) bins[i] = static_cast<int>((i * 2654435761u) % (1u << width)) + 5;
            bins[0] = 5;
            if (n > 1) bins[1] = (1 << width) + 4;
            check_round_trip(bins, 0, 0x01);
        }
    }
}

TEST(SZ3_HuffmanEncoderV2, RandomHistograms) {
    std::mt19937_64 rng(4);
    for (int trial = 0; trial < 400; trial++) {
        const size_t n = trial % 10 == 0 ? 1 + rng() % 8 : 1 + rng() % (trial % 3 == 0 ? 30000 : 2000);
        int stateNum = 0;
        auto bins = random_stream(rng, trial % 4, n, stateNum);
        check_round_trip(bins, stateNum);
    }
}

TEST(SZ3_HuffmanEncoderV2, EncoderObjectsAreReusable) {
    // One encoder through several inputs; one decoder through several streams.
    std::mt19937_64 rng(5);
    V2 enc, dec;
    for (int trial = 0; trial < 30; trial++) {
        int stateNum = 0;
        auto bins = random_stream(rng, trial % 4, 1 + rng() % 5000, stateNum);
        enc.preprocess_encode(bins, stateNum);
        std::vector<uint8_t> buf(enc.size_est() + sizeof(size_t) + 16 + 8 * bins.size());
        uint8_t* p = buf.data();
        enc.save(p);
        SZ3::write<size_t>(bins.size(), p);
        enc.encode(bins, p);
        const size_t written = p - buf.data();
        // Encoding the same input twice from one object gives the same payload.
        std::vector<uint8_t> buf2(16 + 8 * bins.size());
        uint8_t* p2 = buf2.data();
        enc.encode(bins, p2);
        std::vector<uint8_t> tail(buf.data() + written - (p2 - buf2.data()), buf.data() + written);
        EXPECT_TRUE(std::equal(tail.begin(), tail.end(), buf2.begin()));

        const uint8_t* q = buf.data();
        size_t rem = written;
        dec.load(q, rem);
        size_t count = 0;
        SZ3::read(count, q, rem);
        EXPECT_TRUE(dec.decode(q, count, rem) == bins);
        EXPECT_EQ(rem, 0u);
    }
}

TEST(SZ3_HuffmanEncoderV2, BioMDQuantizationStreams) {
    // The bins ALGO_BIOMD hands this encoder: a water-box trajectory through SZBioMDDecomposition.
    std::mt19937 rng(7);
    std::uniform_real_distribution<float> jitter(-0.02f, 0.02f);
    const size_t frames = 6, atoms = 999;
    std::vector<float> traj(frames * atoms * 3), pos(atoms * 3);
    for (size_t a = 0; a < atoms; a++)
        for (int d = 0; d < 3; d++)
            pos[a * 3 + d] = static_cast<float>((a / (d == 0 ? 1 : d == 1 ? 10 : 100)) % 10) * 0.31f;
    for (size_t f = 0; f < frames; f++)
        for (size_t i = 0; i < atoms * 3; i++) traj[f * atoms * 3 + i] = (pos[i] += jitter(rng));

    for (int radius : {32768, 1024, 1 << 20}) {
        for (double eb : {1e-1, 1e-3, 1e-5, 1e-7}) {
            for (int dims : {1, 2, 3}) {
                SZ3::Config conf = dims == 1   ? SZ3::Config(traj.size())
                                   : dims == 2 ? SZ3::Config(frames * atoms, size_t{3})
                                               : SZ3::Config(frames, atoms, size_t{3});
                std::vector<float> data(traj);
                std::vector<int> bins;
                SZ3::LinearQuantizer<float> q(eb, radius);
                if (dims == 1) bins = SZ3::make_decomposition_biomd<float, 1>(conf, q).compress(conf, data.data());
                if (dims == 2) bins = SZ3::make_decomposition_biomd<float, 2>(conf, q).compress(conf, data.data());
                if (dims == 3) bins = SZ3::make_decomposition_biomd<float, 3>(conf, q).compress(conf, data.data());
                check_round_trip(bins, 2 * radius);
            }
        }
    }
}

TEST(SZ3_HuffmanEncoderV2, StreamFormatIsStable) {
    // Pins the bytes of four streams (dense, map, table-path, raw).
    std::vector<int> a(1000), b(3000), c = fibonacci_stream(20, iota_values(20, 3)), d(77);
    for (int i = 0; i < 1000; i++) a[i] = (i * i + 3 * i) % 37 + 12;
    for (int i = 0; i < 3000; i++) b[i] = 32768 + ((i * 7919) % 97) - 48 + (i % 101 == 0 ? 70000 : 0);
    for (int i = 0; i < 77; i++) d[i] = (i * 31) % 1000;
    EXPECT_EQ(fnv1a(check_round_trip(a, 0)), 0x400a81b62de4f5eeull);
    EXPECT_EQ(fnv1a(check_round_trip(b, 65536)), 0xf29de062d2188424ull);
    EXPECT_EQ(fnv1a(check_round_trip(c, 0)), 0x73b66e4c02cbc398ull);
    EXPECT_EQ(fnv1a(check_round_trip(d, 0, 0x01)), 0x20e3cfd7450463adull);
}

TEST(SZ3_HuffmanEncoderV2, EmptyInputThrows) {
    V2 enc;
    EXPECT_THROW(enc.preprocess_encode(std::vector<int>{}, 0), std::invalid_argument);
}

TEST(SZ3_HuffmanEncoderV2, NegativeBinsThrow) {
    for (const auto& bins : {std::vector<int>{-1}, std::vector<int>{0, 5, -1}}) {
        V2 enc;
        EXPECT_THROW(enc.preprocess_encode(bins, 0), std::invalid_argument);
    }
}

TEST(SZ3_HuffmanEncoderV2, SpanWiderThanInt) {
    for (const auto& bins : {std::vector<int>{0, INT_MAX}, std::vector<int>{INT_MAX, 3, 0}}) {
        V2 enc;
        EXPECT_THROW(enc.preprocess_encode(bins, 0), std::invalid_argument);
    }
    check_round_trip({0, INT_MAX - 1}, 0);
    check_round_trip({1, INT_MAX, 1}, 0);
}

TEST(SZ3_HuffmanEncoderV2, LoadRejectsNegativeOffset) {
    Encoded e = check_round_trip({5, 6, 6, 7}, 0);
    // The tree header is a flag byte followed by the offset.
    memset(e.tree.data() + 1, 0xff, sizeof(int));
    V2 dec;
    const uint8_t* p = e.tree.data();
    size_t rem = e.tree.size();
    EXPECT_THROW(dec.load(p, rem), std::out_of_range);
}
