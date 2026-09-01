/**
 * @file test_composition.cpp
 * @brief Composition-level tests for modules that have no `ALGO_*` wiring.
 *
 * fz ships far more modules than algorithms, so several Quantizer / Encoder /
 * Decomposition modules are only ever exercised in isolation (quantize one
 * value, round-trip one bin vector). Isolated tests cannot see
 * composition-level problems: an output range the compressor rejects, a bin
 * width an encoder cannot represent, save/load ordering, a `size_est()` that
 * under-reports and lets the encoder run off the end of
 * `SZGenericCompressor`'s scratch buffer, or an encoder that assumes a zeroed
 * destination the compressor never provides.
 *
 * Each composition here is built with `make_compressor_sz_generic` and driven
 * through a full compress -> lossless -> decompress round trip.
 *
 * Modules covered: BitTruncationQuantizer, FixedPointQuantizer,
 * BitplaneEncoder, BitshuffleEncoder, RunlengthEncoder, ArithmeticEncoder,
 * TimeSeriesDecomposition.
 *
 * Three of them turned out to be defective rather than merely un-wired. Those
 * defects are pre-existing and are NOT fixed here; each has a dedicated test
 * that pins the root cause deterministically and says what to do if the test
 * ever starts failing (which would mean the bug was fixed):
 *
 *   - `BitshuffleEncoder::encode()` never clears its destination.
 *   - `ArithmeticEncoder::decode()` sign-extends the initial code word.
 *   - `RunlengthEncoder` can overrun the compressor's scratch buffer.
 *   - `TimeSeriesDecomposition` loses the timestep-0 reconstruction when
 *     `data_ts0 == nullptr`, which costs up to 2x the requested error bound.
 *
 * This file is deliberately test-only: none of these compositions are wired
 * into `Config::cmprAlgo` or `SZDispatcher`.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <random>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

#include "SZ3/compressor/SZGenericCompressor.hpp"
#include "SZ3/decomposition/BlockwiseDecomposition.hpp"
#include "SZ3/decomposition/Decomposition.hpp"
#include "SZ3/decomposition/NoPredictionDecomposition.hpp"
#include "SZ3/decomposition/TimeSeriesDecomposition.hpp"
#include "SZ3/encoder/ArithmeticEncoder.hpp"
#include "SZ3/encoder/BitplaneEncoder.hpp"
#include "SZ3/encoder/BitshuffleEncoder.hpp"
#include "SZ3/encoder/BypassEncoder.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/encoder/RunlengthEncoder.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/quantizer/BitTruncationQuantizer.hpp"
#include "SZ3/quantizer/FixedPointQuantizer.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "gtest/gtest.h"

namespace {

// ---------------------------------------------------------------------------
// Deterministic data generators (analytic, so runs are reproducible).
// ---------------------------------------------------------------------------
std::vector<float> smooth1Df(size_t n) {
    std::vector<float> v(n);
    for (size_t i = 0; i < n; i++) {
        const double x = static_cast<double>(i) / static_cast<double>(n);
        v[i] = static_cast<float>(std::sin(6.0 * x) + 0.25 * std::cos(31.0 * x));
    }
    return v;
}

std::vector<float> smooth3Df(size_t nz, size_t ny, size_t nx) {
    std::vector<float> v(nz * ny * nx);
    for (size_t z = 0; z < nz; z++) {
        for (size_t y = 0; y < ny; y++) {
            for (size_t x = 0; x < nx; x++) {
                const double a = static_cast<double>(x) * 0.11;
                const double b = static_cast<double>(y) * 0.07;
                const double c = static_cast<double>(z) * 0.13;
                v[(z * ny + y) * nx + x] = static_cast<float>(std::sin(a) * std::cos(b) + 0.3 * std::sin(a + b + c));
            }
        }
    }
    return v;
}

// Strictly positive field in [9, 31]. `BitTruncationQuantizer` puts the raw
// IEEE-754 pattern in the bin, so a sign change moves the bin by ~2^31.
std::vector<float> positive3Df(size_t nz, size_t ny, size_t nx) {
    std::vector<float> v(nz * ny * nx);
    for (size_t z = 0; z < nz; z++) {
        for (size_t y = 0; y < ny; y++) {
            for (size_t x = 0; x < nx; x++) {
                const double a = static_cast<double>(x) * 0.11;
                const double b = static_cast<double>(y) * 0.07;
                const double c = static_cast<double>(z) * 0.13;
                v[(z * ny + y) * nx + x] =
                    static_cast<float>(20.0 + 8.0 * std::sin(a) * std::cos(b) + 3.0 * std::sin(a + b + c));
            }
        }
    }
    return v;
}

// A single-binade band, [1, 1 + 2^-7). Every value shares one exponent, so the
// truncated bit patterns span only 2^16 - the only regime in which
// `HuffmanEncoder` can afford `BitTruncationQuantizer` bins at all (see
// `BitTruncationQuantizerIsUnaffordableForHuffman`).
std::vector<float> narrowBand1Df(size_t n) {
    std::vector<float> v(n);
    for (size_t i = 0; i < n; i++) {
        const double t = static_cast<double>(i) / static_cast<double>(n);
        v[i] = static_cast<float>(1.0 + 0.00390625 * (1.0 + std::sin(6.0 * t)));
    }
    return v;
}

std::vector<double> smooth3Dd(size_t nz, size_t ny, size_t nx) {
    const auto f = smooth3Df(nz, ny, nx);
    return std::vector<double>(f.begin(), f.end());
}

// Deterministic pseudo-random field in [-1, 1]. `std::mt19937` is specified by
// the standard, so this is reproducible across platforms. Used where a
// genuinely incompressible bin stream is needed.
std::vector<float> pseudoRandom1Df(size_t n, unsigned seed) {
    std::mt19937 g(seed);
    std::vector<float> v(n);
    for (auto &x : v) {
        x = static_cast<float>(static_cast<double>(g() % 2000001u) / 1000000.0 - 1.0);
    }
    return v;
}

// A 2D (time x space) field: smooth in space, slowly drifting in time.
std::vector<float> timeSeries2Df(size_t nt, size_t ns) {
    std::vector<float> v(nt * ns);
    for (size_t t = 0; t < nt; t++) {
        for (size_t s = 0; s < ns; s++) {
            const double x = static_cast<double>(s) * 0.05;
            const double u = static_cast<double>(t) * 0.02;
            v[t * ns + s] = static_cast<float>(std::sin(x + u) + 0.2 * std::cos(3.0 * x) + 0.01 * u);
        }
    }
    return v;
}

SZ3::Config absConf(const std::vector<size_t> &dims, double eb, int quantbinCnt = 65536) {
    SZ3::Config conf;
    auto d = dims;
    conf.setDims(d.begin(), d.end());
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = eb;
    conf.quantbinCnt = quantbinCnt;
    return conf;
}

// ---------------------------------------------------------------------------
// Round-trip driver. `make` is called twice: compression and decompression
// need independent module state (a Huffman tree built during encode must not
// leak into decode), exactly as the shipped `SZAlgo*.hpp` wirings do.
// ---------------------------------------------------------------------------
template <class T>
struct PipeResult {
    size_t cmpSize = 0;
    double maxAbsErr = 0.0;
    double maxRelErr = 0.0;
    std::vector<T> dec;
};

template <class T, class Factory>
PipeResult<T> runPipe(const SZ3::Config &conf, const std::vector<T> &orig, Factory make) {
    // Generous capacity: some of these encoders expand rather than compress,
    // and Lossless_zstd refuses any buffer below ZSTD_compressBound().
    const size_t cap = 32 * orig.size() * sizeof(T) + (1u << 16);
    std::vector<SZ3::uchar> cmp(cap, 0);

    std::vector<T> data = orig;
    PipeResult<T> r;
    r.cmpSize = make()->compress(conf, data.data(), cmp.data(), cap);

    r.dec.assign(conf.num, static_cast<T>(0));
    make()->decompress(conf, cmp.data(), r.cmpSize, r.dec.data());

    for (size_t i = 0; i < conf.num; i++) {
        const double e = std::fabs(static_cast<double>(orig[i]) - static_cast<double>(r.dec[i]));
        if (e > r.maxAbsErr) r.maxAbsErr = e;
        const double denom = std::fabs(static_cast<double>(orig[i]));
        if (denom > 0.0) {
            r.maxRelErr = std::max(r.maxRelErr, e / denom);
        }
    }
    return r;
}

template <class T>
void expectBitIdentical(const std::vector<T> &a, const std::vector<T> &b, const char *what) {
    ASSERT_EQ(a.size(), b.size()) << what;
    ASSERT_EQ(0, std::memcmp(a.data(), b.data(), a.size() * sizeof(T)))
        << what
        << ": reconstruction differs from the reference pipeline; a lossless encoder swap must not "
           "change the decompressed data";
}

template <class T>
void reportSize(const char *name, const SZ3::Config &conf, const PipeResult<T> &r) {
    const size_t raw = conf.num * sizeof(T);
    std::printf("[ COMPOSE  ] %-44s cmp=%8zu raw=%8zu ratio=%6.2fx maxAbsErr=%.6g\n", name, r.cmpSize, raw,
                static_cast<double>(raw) / static_cast<double>(r.cmpSize), r.maxAbsErr);
}

// ===========================================================================
// Section 1 - encoder swaps at the standard Q = int width.
//
// Reference pipeline is the shipped one: BlockwiseDecomposition +
// LorenzoPredictor + LinearQuantizer + HuffmanEncoder + Lossless_zstd. Each
// test replaces only the encoder. All of these encoders are lossless, so the
// decompressed data must match the Huffman reference bit for bit.
// ===========================================================================

template <class T, SZ3::uint N, class Encoder>
auto makeLorenzoPipe(const SZ3::Config &conf, Encoder encoder) {
    return SZ3::make_compressor_sz_generic<T, N>(
        SZ3::make_decomposition_blockwise<T, N>(conf, SZ3::LorenzoPredictor<T, N, 1>(conf.absErrorBound),
                                                SZ3::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
        encoder, SZ3::Lossless_zstd());
}

template <class T, SZ3::uint N, class Encoder>
auto makeNopredPipe(const SZ3::Config &conf, Encoder encoder) {
    return SZ3::make_compressor_sz_generic<T, N>(
        SZ3::make_decomposition_noprediction<T, N>(conf,
                                                   SZ3::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
        encoder, SZ3::Lossless_zstd());
}

template <class T, SZ3::uint N, class Encoder>
void checkLosslessEncoderSwap(const SZ3::Config &conf, const std::vector<T> &orig, Encoder encoder, const char *name) {
    const auto ref = runPipe<T>(conf, orig, [&] { return makeLorenzoPipe<T, N>(conf, SZ3::HuffmanEncoder<int>()); });
    const auto got = runPipe<T>(conf, orig, [&] { return makeLorenzoPipe<T, N>(conf, encoder); });

    reportSize("Lorenzo + Huffman (reference)", conf, ref);
    reportSize(name, conf, got);

    EXPECT_LE(ref.maxAbsErr, conf.absErrorBound) << "reference pipeline broke its own bound";
    EXPECT_LE(got.maxAbsErr, conf.absErrorBound) << name << " violated the error bound";
    expectBitIdentical(ref.dec, got.dec, name);
    EXPECT_LT(got.cmpSize, conf.num * sizeof(T)) << name << " did not compress below the raw input size";
}

TEST(SZ3_CompositionTest, BitplaneEncoder1D) {
    const size_t n = 4096;
    const auto conf = absConf({n}, 1e-3);
    checkLosslessEncoderSwap<float, 1>(conf, smooth1Df(n), SZ3::BitplaneEncoder<int>(), "Lorenzo + Bitplane 1D");
}

TEST(SZ3_CompositionTest, BitplaneEncoder3D) {
    const auto conf = absConf({16, 16, 16}, 1e-3);
    checkLosslessEncoderSwap<float, 3>(conf, smooth3Df(16, 16, 16), SZ3::BitplaneEncoder<int>(),
                                       "Lorenzo + Bitplane 3D");
}

// The no-predictor decomposition produces a far wider, higher-entropy bin
// stream than Lorenzo - a much harder case for a plane-oriented encoder.
TEST(SZ3_CompositionTest, BitplaneEncoderOnNoPredictionBins) {
    const size_t n = 4096;
    const auto conf = absConf({n}, 1e-3);
    const auto orig = smooth1Df(n);

    const auto ref =
        runPipe<float>(conf, orig, [&] { return makeNopredPipe<float, 1>(conf, SZ3::HuffmanEncoder<int>()); });
    const auto got =
        runPipe<float>(conf, orig, [&] { return makeNopredPipe<float, 1>(conf, SZ3::BitplaneEncoder<int>()); });

    reportSize("NoPred + Huffman (reference)", conf, ref);
    reportSize("NoPred + Bitplane", conf, got);

    EXPECT_LE(got.maxAbsErr, conf.absErrorBound);
    expectBitIdentical(ref.dec, got.dec, "NoPred + Bitplane");
    EXPECT_LT(got.cmpSize, conf.num * sizeof(float));
}

TEST(SZ3_CompositionTest, RunlengthEncoder1D) {
    // Lorenzo on smooth data puts nearly every sample in the same bin, so the
    // bin stream is a handful of very long runs - the case RLE is built for,
    // and the only case in which it is safe (see the overrun test below).
    const size_t n = 4096;
    const auto conf = absConf({n}, 1e-2);
    checkLosslessEncoderSwap<float, 1>(conf, smooth1Df(n), SZ3::RunlengthEncoder<int>(), "Lorenzo + Runlength 1D");
}

TEST(SZ3_CompositionTest, RunlengthEncoder3D) {
    const auto conf = absConf({16, 16, 16}, 1e-2);
    checkLosslessEncoderSwap<float, 3>(conf, smooth3Df(16, 16, 16), SZ3::RunlengthEncoder<int>(),
                                       "Lorenzo + Runlength 3D");
}

// BitshuffleEncoder behind the zero-fill adapter. Everything else about the
// module composes: types line up, save/load is a no-op, the byte count it
// advances matches what decode consumes, and the result is bit-exact.
TEST(SZ3_CompositionTest, BitshuffleEncoder1D) {
    const size_t n = 4096;
    const auto conf = absConf({n}, 1e-3);
    checkLosslessEncoderSwap<float, 1>(conf, smooth1Df(n), SZ3::BitshuffleEncoder<int>(1),
                                       "Lorenzo + Bitshuffle(bits=1) 1D [zero-filled]");
}

TEST(SZ3_CompositionTest, BitshuffleEncoder3D) {
    const auto conf = absConf({16, 16, 16}, 1e-3);
    checkLosslessEncoderSwap<float, 3>(conf, smooth3Df(16, 16, 16), SZ3::BitshuffleEncoder<int>(8),
                                       "Lorenzo + Bitshuffle(bits=8) 3D [zero-filled]");
}

// ---------------------------------------------------------------------------
// Regression test for a FIXED defect: BitshuffleEncoder::encode() used to only
// OR bits into its destination and never clear it, while SZGenericCompressor
// hands it malloc() memory. Stale heap bits then leaked into decoded values.
// ---------------------------------------------------------------------------
TEST(SZ3_CompositionTest, BitshuffleEncoderToleratesDirtyDestination) {
    const std::vector<int> bins = {0, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144};
    const size_t payload = bins.size() * sizeof(int);

    auto roundTrip = [&](SZ3::uchar fill) {
        std::vector<SZ3::uchar> buf(payload + 64, fill);
        SZ3::uchar *wp = buf.data();
        SZ3::BitshuffleEncoder<int> enc(1);
        enc.preprocess_encode(bins, 256);
        enc.encode(bins, wp);

        const SZ3::uchar *rp = buf.data();
        SZ3::BitshuffleEncoder<int> dec(1);
        return dec.decode(rp, bins.size());
    };

    EXPECT_EQ(bins, roundTrip(0x00)) << "exact round trip on a zeroed destination";
    EXPECT_EQ(bins, roundTrip(0xFF)) << "stale destination bits are leaking through encode()";
    EXPECT_EQ(bins, roundTrip(0xA5)) << "stale destination bits are leaking through encode()";
}

// ---------------------------------------------------------------------------
// Regression test for a FIXED defect: RunlengthEncoder had no size_est()
// override, so SZGenericCompressor under-allocated its scratch buffer. RLE
// emits sizeof(T)+sizeof(int) bytes per run, and an incompressible bin stream
// has one run per sample, which is exactly the whole buffer - while the
// decomposition header and the 8-byte bin count had already been written into
// it. This composes the worst case end to end.
// ---------------------------------------------------------------------------
TEST(SZ3_CompositionTest, RunlengthEncoderSurvivesOneRunPerSample) {
    const size_t n = 4096;
    const auto conf = absConf({n}, 2e-5);
    const auto orig = pseudoRandom1Df(n, 12345);

    // Confirm this input really is the worst case for RLE before relying on it.
    auto decomp = SZ3::make_decomposition_noprediction<float, 1>(
        conf, SZ3::LinearQuantizer<float>(conf.absErrorBound, conf.quantbinCnt / 2));
    std::vector<float> probe = orig;
    const std::vector<int> bins = decomp.compress(conf, probe.data());
    size_t runs = 1;
    for (size_t i = 1; i < bins.size(); i++) {
        if (bins[i] != bins[i - 1]) runs++;
    }
    EXPECT_EQ(runs, n) << "expected one run per sample on an incompressible bin stream";

    const auto ref =
        runPipe<float>(conf, orig, [&] { return makeNopredPipe<float, 1>(conf, SZ3::HuffmanEncoder<int>()); });
    const auto got =
        runPipe<float>(conf, orig, [&] { return makeNopredPipe<float, 1>(conf, SZ3::RunlengthEncoder<int>()); });

    reportSize("NoPred + Huffman (reference)", conf, ref);
    reportSize("NoPred + Runlength (worst case)", conf, got);

    EXPECT_LE(got.maxAbsErr, conf.absErrorBound) << "Runlength pipeline violated the error bound";
    expectBitIdentical(ref.dec, got.dec, "NoPred + Runlength");
}

// ---------------------------------------------------------------------------
// Regression test for a FIXED defect: ArithmeticEncoder::decode() primed its
// code word with a SIGNED right shift, so any stream whose first coded bit was
// 1 - roughly half of all inputs - sign-extended outside the 44-bit MAX_CODE
// window and desynced on the first symbol. A prior sweep over 60 random
// streams gave 30 correct, 23 wrong, 7 segfaults. The single input in
// test_encoder.cpp (i % 100) happens to start with a 0 bit, which is why it
// never caught this.
// ---------------------------------------------------------------------------
TEST(SZ3_CompositionTest, ArithmeticEncoderRoundTripsLeadingOneBitStreams) {
    auto roundTrip = [](const std::vector<int> &bins, int stateNum, unsigned *firstByte) {
        std::vector<SZ3::uchar> buf(1u << 20, 0);
        SZ3::uchar *p = buf.data();
        SZ3::ArithmeticEncoder<int> enc(false);
        enc.preprocess_encode(bins, stateNum);
        enc.save(p);
        const SZ3::uchar *dataStart = p;
        enc.encode(bins, p);
        enc.postprocess_encode();
        if (firstByte) *firstByte = *dataStart;

        const SZ3::uchar *rp = buf.data();
        SZ3::ArithmeticEncoder<int> dec(false);
        size_t remaining = buf.size();
        dec.load(rp, remaining);
        auto out = dec.decode(rp, bins.size());
        dec.postprocess_decode();
        return out;
    };

    // The specific input that used to fail: {2, 0} over a 4-symbol alphabet
    // starts in the upper half of the interval, so the first code bit is 1.
    unsigned fb = 0;
    const std::vector<int> tricky = {2, 0};
    EXPECT_EQ(tricky, roundTrip(tricky, 4, &fb));
    EXPECT_NE(0u, fb & 0x80u) << "sanity: this input must still produce a leading 1 bit";

    // Randomized sweep. Before the fix this was ~50% correct.
    std::mt19937 gen(20260831);
    size_t leadingOne = 0;
    for (int trial = 0; trial < 60; trial++) {
        const int stateNum = 4 + static_cast<int>(gen() % 61);
        const size_t n = 8 + gen() % 200;
        std::vector<int> bins(n);
        for (size_t i = 0; i < n; i++) bins[i] = static_cast<int>(gen() % static_cast<unsigned>(stateNum));
        unsigned f = 0;
        EXPECT_EQ(bins, roundTrip(bins, stateNum, &f)) << "trial " << trial << " stateNum=" << stateNum;
        if (f & 0x80u) leadingOne++;
    }
    std::printf("[ COMPOSE  ] %-44s %zu/60 streams had a leading 1 bit\n", "Arithmetic randomized sweep", leadingOne);
    EXPECT_GT(leadingOne, 0u) << "sweep never produced the failing case; it is not exercising the fix";
}

// ===========================================================================
// Section 2 - FixedPointQuantizer (int64_t bins) through the pipeline.
// ===========================================================================

template <class T, SZ3::uint N, class Encoder>
auto makeFixedPointPipe(const SZ3::FixedPointQuantizer<T> &q, Encoder encoder) {
    return SZ3::make_compressor_sz_generic<T, N>(SZ3::make_decomposition_noprediction<T, N>(SZ3::Config(), q), encoder,
                                                 SZ3::Lossless_zstd());
}

template <class T>
double maxAbs(const std::vector<T> &v) {
    double m = 0.0;
    for (const T &x : v) m = std::max(m, std::fabs(static_cast<double>(x)));
    return m;
}

TEST(SZ3_CompositionTest, FixedPointQuantizer3DDouble) {
    const auto conf = absConf({16, 16, 16}, 1e-3);
    const auto orig = smooth3Dd(16, 16, 16);

    SZ3::FixedPointQuantizer<double> q(16);
    q.calibrate(maxAbs(orig));
    const double bound = q.max_abs_error();

    const auto bypass =
        runPipe<double>(conf, orig, [&] { return makeFixedPointPipe<double, 3>(q, SZ3::BypassEncoder<int64_t>()); });
    const auto plane =
        runPipe<double>(conf, orig, [&] { return makeFixedPointPipe<double, 3>(q, SZ3::BitplaneEncoder<int64_t>()); });

    reportSize("FixedPoint(16b) + Bypass", conf, bypass);
    reportSize("FixedPoint(16b) + Bitplane", conf, plane);
    std::printf("[ COMPOSE  ] FixedPoint(16b) max_abs_error() = %.6g\n", bound);

    EXPECT_LE(bypass.maxAbsErr, bound) << "FixedPointQuantizer exceeded its own max_abs_error()";
    EXPECT_LE(plane.maxAbsErr, bound);
    expectBitIdentical(bypass.dec, plane.dec, "FixedPoint + Bitplane vs Bypass");
    EXPECT_LT(bypass.cmpSize, conf.num * sizeof(double));
    EXPECT_LT(plane.cmpSize, conf.num * sizeof(double));
}

TEST(SZ3_CompositionTest, FixedPointQuantizer1DFloat) {
    const size_t n = 4096;
    const auto conf = absConf({n}, 1e-3);
    const auto orig = smooth1Df(n);

    SZ3::FixedPointQuantizer<float> q(12);
    q.calibrate(maxAbs(orig));
    // The quantizer rounds in double and stores the reconstruction as float,
    // so allow a little slack on top of the theoretical bound.
    const double bound = q.max_abs_error() * (1.0 + 1e-5);

    const auto bypass =
        runPipe<float>(conf, orig, [&] { return makeFixedPointPipe<float, 1>(q, SZ3::BypassEncoder<int64_t>()); });
    const auto plane =
        runPipe<float>(conf, orig, [&] { return makeFixedPointPipe<float, 1>(q, SZ3::BitplaneEncoder<int64_t>()); });
    const auto huff =
        runPipe<float>(conf, orig, [&] { return makeFixedPointPipe<float, 1>(q, SZ3::HuffmanEncoder<int64_t>()); });

    reportSize("FixedPoint(12b) + Bypass", conf, bypass);
    reportSize("FixedPoint(12b) + Bitplane", conf, plane);
    reportSize("FixedPoint(12b) + Huffman", conf, huff);
    std::printf("[ COMPOSE  ] FixedPoint(12b) max_abs_error() = %.6g\n", q.max_abs_error());

    EXPECT_LE(bypass.maxAbsErr, bound);
    EXPECT_LE(plane.maxAbsErr, bound);
    EXPECT_LE(huff.maxAbsErr, bound);
    expectBitIdentical(bypass.dec, plane.dec, "FixedPoint + Bitplane vs Bypass");
    expectBitIdentical(bypass.dec, huff.dec, "FixedPoint + Huffman vs Bypass");

    const size_t raw = conf.num * sizeof(float);
    EXPECT_LT(bypass.cmpSize, raw);
    EXPECT_LT(plane.cmpSize, raw);
    // Honest negative result: HuffmanEncoder EXPANDS this pipeline past the raw
    // input. `FixedPointQuantizer` is a scalar quantizer with no predictor, so
    // its 12-bit bins are near-uniform over ~2500 distinct values; the
    // serialised Huffman tree (~13 bytes per node, 2*distinct-1 nodes) then
    // costs more than the whole payload. This is the general shape of the
    // problem: Huffman is the wrong partner for wide, high-entropy bins.
    EXPECT_GT(huff.cmpSize, raw)
        << "FixedPoint(12b) + Huffman now compresses below the raw input; if HuffmanEncoder's tree "
           "serialisation got cheaper, turn this into a plain EXPECT_LT";
}

// ===========================================================================
// Section 3 - BitTruncationQuantizer, paired with BitplaneEncoder.
//
// The bin is the raw IEEE-754 bit pattern as uint64_t, so it is never negative
// and the output-range contract holds for both float and double. T = float here.
// ===========================================================================

template <class T, SZ3::uint N, class Encoder>
auto makeBitTruncPipe(int keepBytes, Encoder encoder) {
    return SZ3::make_compressor_sz_generic<T, N>(
        SZ3::make_decomposition_noprediction<T, N>(SZ3::Config(), SZ3::BitTruncationQuantizer<T>(keepBytes)), encoder,
        SZ3::Lossless_zstd());
}

// Independent model of the quantizer: zero the low (4 - keepBytes) bytes.
float truncateFloat(float x, int keepBytes) {
    uint32_t bits;
    std::memcpy(&bits, &x, sizeof(bits));
    const int dropBits = (4 - keepBytes) * 8;
    bits &= (dropBits >= 32) ? 0u : ~((uint32_t{1} << dropBits) - 1);
    float out;
    std::memcpy(&out, &bits, sizeof(out));
    return out;
}

// The pairing BitplaneEncoder.hpp advertises, measured against the pairing it
// advertises first (BypassEncoder + Lossless_zstd) on a realistic field.
TEST(SZ3_CompositionTest, BitTruncationQuantizerBitplaneVsBypass) {
    const auto conf = absConf({16, 16, 16}, 1.0 /*unused by this pipeline*/);
    const auto orig = positive3Df(16, 16, 16);
    const int keepBytes = 2;  // sign + 8 exponent bits + 7 mantissa bits

    const auto bypass = runPipe<float>(
        conf, orig, [&] { return makeBitTruncPipe<float, 3>(keepBytes, SZ3::BypassEncoder<uint64_t>()); });
    const auto plane = runPipe<float>(
        conf, orig, [&] { return makeBitTruncPipe<float, 3>(keepBytes, SZ3::BitplaneEncoder<uint64_t>()); });

    reportSize("BitTrunc(keep=2) + Bypass", conf, bypass);
    reportSize("BitTrunc(keep=2) + Bitplane", conf, plane);
    std::printf("[ COMPOSE  ] BitTrunc Bitplane/Bypass = %.3f (<1 means Bitplane wins)\n",
                static_cast<double>(plane.cmpSize) / static_cast<double>(bypass.cmpSize));

    // Truncation matches an independent model of the transform, exactly.
    for (size_t i = 0; i < conf.num; i++) {
        ASSERT_EQ(truncateFloat(orig[i], keepBytes), plane.dec[i]) << "at " << i;
    }
    // 7 mantissa bits kept, truncated toward zero -> relative error < 2^-7.
    EXPECT_LT(plane.maxRelErr, std::ldexp(1.0, -7));

    expectBitIdentical(bypass.dec, plane.dec, "BitTrunc + Bitplane vs Bypass");

    const size_t raw = conf.num * sizeof(float);
    EXPECT_LT(plane.cmpSize, raw) << "BitTrunc + Bitplane did not beat the raw input";
    EXPECT_LT(bypass.cmpSize, raw) << "BitTrunc + Bypass did not beat the raw input";
}

// Mixed-sign data: the bit-pattern bins jump by ~2^31 across zero, which is
// the case most likely to break the output-range contract. For float it holds
// (the 32-bit pattern widens to a non-negative int64_t) and BitplaneEncoder
// just needs 32 planes. HuffmanEncoder is not viable here at all - see
// `BitTruncationQuantizerIsUnaffordableForHuffman`.
TEST(SZ3_CompositionTest, BitTruncationQuantizerMixedSign) {
    const size_t n = 4096;
    const auto conf = absConf({n}, 1.0 /*unused by this pipeline*/);
    const auto orig = smooth1Df(n);  // spans negative and positive
    const int keepBytes = 3;

    const auto bypass = runPipe<float>(
        conf, orig, [&] { return makeBitTruncPipe<float, 1>(keepBytes, SZ3::BypassEncoder<uint64_t>()); });
    const auto plane = runPipe<float>(
        conf, orig, [&] { return makeBitTruncPipe<float, 1>(keepBytes, SZ3::BitplaneEncoder<uint64_t>()); });

    reportSize("BitTrunc(keep=3) mixed-sign + Bypass", conf, bypass);
    reportSize("BitTrunc(keep=3) mixed-sign + Bitplane", conf, plane);

    for (size_t i = 0; i < conf.num; i++) {
        ASSERT_EQ(truncateFloat(orig[i], keepBytes), plane.dec[i]) << "at " << i;
    }
    EXPECT_LT(plane.maxRelErr, std::ldexp(1.0, -15));  // 15 mantissa bits kept
    expectBitIdentical(bypass.dec, plane.dec, "BitTrunc mixed-sign Bitplane vs Bypass");
    EXPECT_LT(plane.cmpSize, conf.num * sizeof(float));
    EXPECT_LT(bypass.cmpSize, conf.num * sizeof(float));
}

// The literal Bitplane-vs-Huffman comparison, on the only kind of field where
// Huffman can afford these bins: a single binade, so the truncated patterns
// span 2^16 rather than millions. See the next test for why.
TEST(SZ3_CompositionTest, BitTruncationQuantizerBitplaneVsHuffmanNarrowBand) {
    const size_t n = 4096;
    const auto conf = absConf({n}, 1.0 /*unused by this pipeline*/);
    const auto orig = narrowBand1Df(n);
    const int keepBytes = 3;

    const auto huff = runPipe<float>(
        conf, orig, [&] { return makeBitTruncPipe<float, 1>(keepBytes, SZ3::HuffmanEncoder<uint64_t>()); });
    const auto plane = runPipe<float>(
        conf, orig, [&] { return makeBitTruncPipe<float, 1>(keepBytes, SZ3::BitplaneEncoder<uint64_t>()); });

    reportSize("BitTrunc(keep=3) narrow + Huffman", conf, huff);
    reportSize("BitTrunc(keep=3) narrow + Bitplane", conf, plane);
    std::printf("[ COMPOSE  ] BitTrunc Bitplane/Huffman = %.3f (<1 means Bitplane wins)\n",
                static_cast<double>(plane.cmpSize) / static_cast<double>(huff.cmpSize));

    for (size_t i = 0; i < conf.num; i++) {
        ASSERT_EQ(truncateFloat(orig[i], keepBytes), plane.dec[i]) << "at " << i;
    }
    expectBitIdentical(huff.dec, plane.dec, "BitTrunc + Bitplane vs Huffman");
    EXPECT_LT(plane.cmpSize, conf.num * sizeof(float));
}

// ---------------------------------------------------------------------------
// FINDING: `HuffmanEncoder` is not a usable partner for
// `BitTruncationQuantizer` on realistic data, and is unsafe on signed data.
//
// `HuffmanEncoder::init()` sizes its state table by the bin *span*
// (`stateNum = max - offset + 2`), not by the number of distinct bins, and
// `createHuffmanTree()` then allocates roughly 137 bytes per state (pool, qqq,
// code, cout). BitTruncation bins are raw IEEE-754 patterns, so any variation
// in the exponent field spans millions of states, and mixed-sign data spans
// ~2^31 - which overflows the `int stateNum` outright.
//
// Asserted arithmetically; instantiating the tree is the thing to avoid.
// ---------------------------------------------------------------------------
TEST(SZ3_CompositionTest, BitTruncationQuantizerIsUnaffordableForHuffman) {
    auto span = [](const std::vector<float> &v, int keepBytes) {
        std::set<int64_t> distinct;
        int64_t lo = std::numeric_limits<int64_t>::max();
        int64_t hi = std::numeric_limits<int64_t>::min();
        SZ3::BitTruncationQuantizer<float> q(keepBytes);
        for (float x : v) {
            float tmp = x;
            const int64_t b = q.quantize_and_overwrite(tmp, 0.0f);
            distinct.insert(b);
            lo = std::min(lo, b);
            hi = std::max(hi, b);
        }
        return std::make_pair(hi - lo + 2, static_cast<int64_t>(distinct.size()));
    };

    const auto pos = span(positive3Df(16, 16, 16), 2);
    const auto mixed = span(smooth1Df(4096), 2);

    std::printf("[ COMPOSE  ] %-44s stateNum=%lld distinct=%lld huffTree~%.0f MB\n",
                "BitTrunc bins, positive field [9,31]", static_cast<long long>(pos.first),
                static_cast<long long>(pos.second), 137.0 * static_cast<double>(pos.first) / 1e6);
    std::printf("[ COMPOSE  ] %-44s stateNum=%lld distinct=%lld huffTree~%.0f MB\n",
                "BitTrunc bins, mixed-sign field [-1,1]", static_cast<long long>(mixed.first),
                static_cast<long long>(mixed.second), 137.0 * static_cast<double>(mixed.first) / 1e6);

    // A few hundred distinct symbols, but a state table millions of entries wide.
    EXPECT_LT(pos.second, 1000);
    EXPECT_GT(pos.first, 1000000) << "HuffmanEncoder's state table is no longer sized by the bin span; a "
                                     "direct BitTrunc + Huffman composition may now be affordable";
    // Mixed-sign bins do not even fit the `int stateNum` HuffmanEncoder uses.
    EXPECT_GT(mixed.first, static_cast<int64_t>(std::numeric_limits<int>::max()));
}

// ===========================================================================
// Section 4 - TimeSeriesDecomposition through SZGenericCompressor.
//
// Dimension semantics differ from every other decomposition: dims must be 2D,
// dims[0] is time and dims[1] is space, and the predictor/quantizer work on
// N-1 spatial dimensions. `data_ts0` is an optional reference frame for
// timestep 0; when null, timestep 0 is predicted spatially instead. This
// mirrors `make_sz_timebased2` in tools/mdz/include/mdz.hpp, the only place
// the module is reachable from today.
// ===========================================================================

template <class T>
auto makeTimeSeriesPipe(const SZ3::Config &conf, T *dataTs0) {
    return SZ3::make_compressor_sz_generic<T, 2>(
        SZ3::make_decomposition_timeseries<T, 2>(conf, SZ3::LorenzoPredictor<T, 1, 1>(conf.absErrorBound),
                                                 SZ3::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2),
                                                 dataTs0),
        SZ3::HuffmanEncoder<int>(), SZ3::Lossless_zstd());
}

TEST(SZ3_CompositionTest, TimeSeriesDecompositionWithReferenceFrame) {
    const size_t nt = 32, ns = 128;
    const auto conf = absConf({nt, ns}, 1e-3);
    const auto orig = timeSeries2Df(nt, ns);

    // Reference frame for timestep 0, identical on both sides as the caller
    // contract requires. Held outside runPipe so both factory calls see it.
    std::vector<float> ts0(ns);
    for (size_t s = 0; s < ns; s++) {
        ts0[s] = orig[s] * 0.98f;
    }

    const auto r = runPipe<float>(conf, orig, [&] { return makeTimeSeriesPipe<float>(conf, ts0.data()); });
    reportSize("TimeSeries (with data_ts0)", conf, r);

    EXPECT_LE(r.maxAbsErr, conf.absErrorBound) << "TimeSeriesDecomposition violated the error bound";
    EXPECT_LT(r.cmpSize, conf.num * sizeof(float));
}

// ---------------------------------------------------------------------------
// Regression test for a FIXED defect: with data_ts0 == nullptr, timestep 0 was
// compressed spatially through `block_data` constructed with data_valid = true,
// which works on an internal padded copy and never writes back. So compress()
// predicted timestep 1 from timestep 0's ORIGINAL values while decompress()
// predicted it from the RECONSTRUCTED ones, and the mismatch added a full
// quantization step on top of timestep 1's own error - up to 2x the bound
// (measured 1.97x). TimeSeriesDecomposition now mirrors the reconstruction back
// into `data` inside the quantize loop.
//
// This is the path mdz takes for current_method == 4 ("TS") - see
// MDZ_Compress() in tools/mdz/include/mdz.hpp.
// ---------------------------------------------------------------------------
TEST(SZ3_CompositionTest, TimeSeriesDecompositionNullReferenceFrameHonoursBound) {
    const size_t nt = 32, ns = 128;
    const auto conf = absConf({nt, ns}, 1e-3);
    const auto orig = timeSeries2Df(nt, ns);

    const auto r = runPipe<float>(conf, orig, [&] { return makeTimeSeriesPipe<float>(conf, nullptr); });
    reportSize("TimeSeries (data_ts0 = nullptr)", conf, r);
    std::printf("[ COMPOSE  ] TimeSeries nullptr-ts0 maxAbsErr / eb = %.4f\n", r.maxAbsErr / conf.absErrorBound);

    EXPECT_LE(r.maxAbsErr, conf.absErrorBound + 1e-9)
        << "data_ts0 == nullptr must honour the requested bound; the timestep-0 write-back is broken again";
    EXPECT_LT(r.cmpSize, conf.num * sizeof(float));
}

}  // namespace
