// zfp CODEC 5 modules: ZFPDecomposition + ZFPEncoder.
//
// The stronger check -- that the compressed stream is byte-identical to upstream
// zfp 1.0.1 -- needs libzfp to link against and so is not run here. These cover the
// properties the pair has to hold on its own: the round trip is exact through the
// seam, the error bound is honoured, partial blocks in every dimension work, and the
// two stages agree on the intermediate layout.

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <cmath>
#include <random>
#include <vector>

#include "SZ3/api/sz.hpp"
#include "SZ3/decomposition/ZFPDecomposition.hpp"
#include "SZ3/encoder/ZFPEncoder.hpp"
#include "gtest/gtest.h"

namespace {

SZ3::Config makeConf(size_t nx, size_t ny, size_t nz, double eb) {
    SZ3::Config conf;
    std::vector<size_t> d{nx, ny, nz};
    conf.setDims(d.begin(), d.end());
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = eb;
    return conf;
}

enum Content { Smooth, Noise, Sparse, Zero };

std::vector<float> makeField(size_t nx, size_t ny, size_t nz, Content c) {
    std::vector<float> d(nx * ny * nz);
    std::mt19937 g(1234);
    std::uniform_real_distribution<float> u(-1, 1);
    for (size_t i = 0; i < d.size(); i++) {
        const double x = double(i % nx) / nx, y = double((i / nx) % ny) / ny, z = double(i / (nx * ny)) / nz;
        switch (c) {
            case Smooth: d[i] = float(std::sin(6 * x) * std::cos(5 * y) + 0.3 * std::sin(9 * z)); break;
            case Noise:  d[i] = u(g); break;
            case Sparse: d[i] = (i % 97 == 0) ? 1e3f : 0.f; break;
            case Zero:   d[i] = 0.f; break;
        }
    }
    return d;
}

/// Drive the pair end to end and report the largest absolute error.
double roundTrip(size_t nx, size_t ny, size_t nz, double eb, Content content, size_t *cmp_bytes = nullptr) {
    const auto conf = makeConf(nx, ny, nz, eb);
    const auto original = makeField(nx, ny, nz, content);
    std::vector<float> work = original, out(original.size(), 0.f);
    std::vector<SZ3::uchar> buf(original.size() * sizeof(float) * 2 + (1u << 16));

    SZ3::ZFPDecomposition<float, int, 3> dec;
    SZ3::ZFPEncoder<int, 3> enc(conf);
    auto coeffs = dec.compress(conf, work.data());
    enc.preprocess_encode(coeffs, 0);
    SZ3::uchar *pos = buf.data();
    enc.encode(coeffs, pos);
    if (cmp_bytes) *cmp_bytes = static_cast<size_t>(pos - buf.data());

    SZ3::ZFPDecomposition<float, int, 3> dec2;
    SZ3::ZFPEncoder<int, 3> enc2(conf);
    const SZ3::uchar *rp = buf.data();
    size_t rlen = buf.size();
    auto back = enc2.decode(rp, coeffs.size(), rlen);
    dec2.decompress(conf, back, out.data());

    double worst = 0;
    for (size_t i = 0; i < original.size(); i++)
        worst = std::max(worst, std::fabs(double(original[i]) - out[i]));
    return worst;
}

TEST(SZ3_ZFP, HonoursTheErrorBound) {
    for (double eb : {1e-1, 1e-2, 1e-3, 1e-4}) {
        for (Content c : {Smooth, Noise, Sparse}) {
            EXPECT_LE(roundTrip(16, 16, 16, eb, c), eb) << "eb=" << eb << " content=" << int(c);
        }
    }
}

TEST(SZ3_ZFP, HandlesPartialBlocksInEveryDimension) {
    // Sizes chosen so each dimension in turn is not a multiple of 4, and one where none is.
    // Every extent must exceed 1: Config::setDims drops extents of 1, which would make conf.N
    // disagree with the N this module is instantiated with.
    const size_t shapes[][3] = {{13, 16, 16}, {16, 7, 16}, {16, 16, 11}, {13, 7, 11}, {5, 3, 2}, {2, 2, 2}};
    for (const auto &s : shapes) {
        EXPECT_LE(roundTrip(s[0], s[1], s[2], 1e-2, Smooth), 1e-2)
            << s[0] << "x" << s[1] << "x" << s[2];
    }
}

TEST(SZ3_ZFP, AllZeroFieldIsExactAndCheap) {
    size_t bytes = 0;
    EXPECT_EQ(roundTrip(16, 16, 16, 1e-2, Zero, &bytes), 0.0);
    // One flag bit per block, padded to the stream's word size.
    EXPECT_LE(bytes, 64u);
}

TEST(SZ3_ZFP, TighterBoundCostsMoreBytes) {
    size_t loose = 0, tight = 0;
    roundTrip(16, 16, 16, 1e-2, Smooth, &loose);
    roundTrip(16, 16, 16, 1e-5, Smooth, &tight);
    EXPECT_GT(tight, loose);
}

TEST(SZ3_ZFP, DecompositionCarriesNoBinRange) {
    // The bound lives in the encoder, so the transform stage declares no usable range.
    SZ3::ZFPDecomposition<float, int, 3> dec;
    EXPECT_EQ(dec.get_out_range().first, 0);
    EXPECT_EQ(dec.get_out_range().second, 0);
}

TEST(SZ3_ZFP, EncoderSizeEstimateCoversTheOutput) {
    const auto conf = makeConf(16, 16, 16, 1e-5);
    const auto original = makeField(16, 16, 16, Noise);
    std::vector<float> work = original;
    SZ3::ZFPDecomposition<float, int, 3> dec;
    SZ3::ZFPEncoder<int, 3> enc(conf);
    auto coeffs = dec.compress(conf, work.data());
    enc.preprocess_encode(coeffs, 0);
    std::vector<SZ3::uchar> buf(enc.size_est());
    SZ3::uchar *pos = buf.data();
    enc.encode(coeffs, pos);
    EXPECT_LE(static_cast<size_t>(pos - buf.data()), enc.size_est());
}

}  // namespace

// ---------------------------------------------------------------------------------------------
// Coverage the six tests above do not reach: the other scalar type, the other dimensionalities,
// the ALGO_ZFP wiring, anisotropic shapes, and streams that have been tampered with.
// ---------------------------------------------------------------------------------------------

namespace {

/// Round trip through whichever (scalar, dimensionality) pair is asked for.
template <class T, class Int, SZ3::uint N>
double roundTripND(std::vector<size_t> dims, double eb) {
    SZ3::Config conf;
    conf.setDims(dims.begin(), dims.end());
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = eb;
    std::vector<T> original(conf.num);
    for (size_t i = 0; i < original.size(); i++) original[i] = static_cast<T>(std::sin(i * 0.05) * 3.0);
    std::vector<T> work = original, out(original.size(), T{0});
    std::vector<SZ3::uchar> buf(original.size() * sizeof(T) * 2 + (1u << 16));

    SZ3::ZFPDecomposition<T, Int, N> dec;
    SZ3::ZFPEncoder<Int, N> enc(conf);
    auto coeffs = dec.compress(conf, work.data());
    enc.preprocess_encode(coeffs, 0);
    SZ3::uchar *pos = buf.data();
    enc.encode(coeffs, pos);

    SZ3::ZFPDecomposition<T, Int, N> dec2;
    SZ3::ZFPEncoder<Int, N> enc2(conf);
    const SZ3::uchar *rp = buf.data();
    size_t rlen = buf.size();
    auto back = enc2.decode(rp, coeffs.size(), rlen);
    dec2.decompress(conf, back, out.data());

    double worst = 0;
    for (size_t i = 0; i < original.size(); i++)
        worst = std::max(worst, std::fabs(static_cast<double>(original[i]) - static_cast<double>(out[i])));
    return worst;
}

}  // namespace

TEST(SZ3_ZFP, CoversBothScalarTypesAndAllDimensionalities) {
    for (double eb : {1e-2, 1e-4}) {
        EXPECT_LE((roundTripND<float, int32_t, 1>({4096}, eb)), eb) << "float 1D eb=" << eb;
        EXPECT_LE((roundTripND<float, int32_t, 2>({64, 67}, eb)), eb) << "float 2D eb=" << eb;
        EXPECT_LE((roundTripND<float, int32_t, 3>({17, 33, 21}, eb)), eb) << "float 3D eb=" << eb;
        EXPECT_LE((roundTripND<double, int64_t, 1>({4096}, eb)), eb) << "double 1D eb=" << eb;
        EXPECT_LE((roundTripND<double, int64_t, 2>({64, 67}, eb)), eb) << "double 2D eb=" << eb;
        EXPECT_LE((roundTripND<double, int64_t, 3>({17, 33, 21}, eb)), eb) << "double 3D eb=" << eb;
    }
}

TEST(SZ3_ZFP, BlocksFollowTheArrayLayout) {
    // SZ3 is row-major, so a smooth isotropic field should compress about equally well however the
    // extents are ordered. If the module's axes were reversed relative to SZ3's, the 4^N blocks
    // would stop being neighbourhoods and anisotropic shapes would cost far more bytes.
    const double eb = 1e-3;
    std::vector<size_t> sizes;
    for (auto dims : std::vector<std::vector<size_t>>{{16, 32, 64}, {64, 32, 16}, {8, 16, 256}}) {
        SZ3::Config conf;
        conf.setDims(dims.begin(), dims.end());
        conf.errorBoundMode = SZ3::EB_ABS;
        conf.absErrorBound = eb;
        std::vector<float> original(conf.num);
        size_t i = 0;
        for (size_t z = 0; z < dims[0]; z++)
            for (size_t y = 0; y < dims[1]; y++)
                for (size_t x = 0; x < dims[2]; x++)
                    original[i++] = static_cast<float>(std::sin(0.05 * x) + std::cos(0.05 * y) + std::sin(0.05 * z));
        std::vector<float> work = original;
        std::vector<SZ3::uchar> buf(conf.num * sizeof(float) * 2 + (1u << 16));
        SZ3::ZFPDecomposition<float, int, 3> dec;
        SZ3::ZFPEncoder<int, 3> enc(conf);
        auto coeffs = dec.compress(conf, work.data());
        enc.preprocess_encode(coeffs, 0);
        SZ3::uchar *pos = buf.data();
        enc.encode(coeffs, pos);
        sizes.push_back(static_cast<size_t>(pos - buf.data()));
    }
    const size_t lo = *std::min_element(sizes.begin(), sizes.end());
    const size_t hi = *std::max_element(sizes.begin(), sizes.end());
    EXPECT_LT(static_cast<double>(hi) / static_cast<double>(lo), 1.25)
        << "same field, permuted extents: " << sizes[0] << " / " << sizes[1] << " / " << sizes[2];
}

TEST(SZ3_ZFP, RejectsTamperedStreams) {
    SZ3::Config conf(16, 16, 16);
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = 1e-3;
    const auto original = makeField(16, 16, 16, Smooth);
    std::vector<float> work = original;
    std::vector<SZ3::uchar> buf(1u << 20);

    SZ3::ZFPDecomposition<float, int, 3> dec;
    SZ3::ZFPEncoder<int, 3> enc(conf);
    auto coeffs = dec.compress(conf, work.data());
    enc.preprocess_encode(coeffs, 0);
    SZ3::uchar *pos = buf.data();
    enc.encode(coeffs, pos);

    // A coefficient count that is not a whole number of blocks is not a stream this encoder wrote.
    {
        SZ3::ZFPEncoder<int, 3> e(conf);
        const SZ3::uchar *rp = buf.data();
        size_t rlen = buf.size();
        EXPECT_THROW(e.decode(rp, coeffs.size() + 1, rlen), std::out_of_range);
    }
    // Nor is one claiming far more blocks than the coded bytes could hold.
    {
        SZ3::ZFPEncoder<int, 3> e(conf);
        const SZ3::uchar *rp = buf.data();
        size_t rlen = buf.size();
        EXPECT_THROW(e.decode(rp, coeffs.size() * 10, rlen), std::out_of_range);
    }
    // A coded length larger than the block count can produce must not be believed.
    {
        std::vector<SZ3::uchar> tampered(buf.begin(), buf.begin() + (pos - buf.data()));
        const uint64_t huge = 1ull << 30;
        std::memcpy(tampered.data(), &huge, sizeof(huge));
        SZ3::ZFPEncoder<int, 3> e(conf);
        const SZ3::uchar *rp = tampered.data();
        size_t rlen = tampered.size();
        EXPECT_THROW(e.decode(rp, coeffs.size(), rlen), std::out_of_range);
    }
    // A count that stays block-aligned but claims more blocks than the stream holds: the
    // capacity check alone scales with it, so the stream's own block count must be consulted.
    {
        SZ3::ZFPEncoder<int, 3> e(conf);
        const SZ3::uchar *rp = buf.data();
        const size_t blocks = (coeffs.size() - 1) / 65;
        size_t rlen = buf.size();
        EXPECT_THROW(e.decode(rp, 1 + (blocks * 4) * 65, rlen), std::out_of_range);
    }
    // A block count that disagrees with the configuration is caught by the decomposition.
    {
        auto bad = coeffs;
        bad[0] = static_cast<int>(coeffs[0]) + 1;
        std::vector<float> out(conf.num);
        SZ3::ZFPDecomposition<float, int, 3> d;
        EXPECT_THROW(d.decompress(conf, bad, out.data()), std::out_of_range);
    }
}

TEST(SZ3_ZFP, RejectsAConfigurationOfTheWrongRank) {
    // Config::setDims drops extents of 1, so a 3D module handed such a Config would otherwise index
    // past the end of conf.dims.
    SZ3::Config conf(16, 1, 16);  // becomes 2D
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = 1e-3;
    std::vector<float> data(conf.num, 1.f);
    SZ3::ZFPDecomposition<float, int, 3> dec;
    EXPECT_THROW(dec.compress(conf, data.data()), std::invalid_argument);
}

TEST(SZ3_ZFP, AlgoZfpRoundTripsThroughTheApi) {
    // The modules are wired together by SZAlgoZFP; nothing above exercises that wiring, or checks
    // that the error bound survives being written into the stream and read back.
    for (int nd = 1; nd <= 3; nd++) {
        std::vector<size_t> dims = nd == 1   ? std::vector<size_t>{4096}
                                   : nd == 2 ? std::vector<size_t>{64, 67}
                                             : std::vector<size_t>{17, 33, 21};
        size_t n = 1;
        for (auto d : dims) n *= d;
        std::vector<float> original(n);
        for (size_t i = 0; i < n; i++) original[i] = static_cast<float>(std::sin(i * 0.01) * 3.0);

        for (double eb : {1e-2, 1e-4}) {
            SZ3::Config conf;
            conf.setDims(dims.begin(), dims.end());
            conf.cmprAlgo = SZ3::ALGO_ZFP;
            conf.errorBoundMode = SZ3::EB_ABS;
            conf.absErrorBound = eb;

            size_t cmpSize = 0;
            char *cmp = SZ_compress(conf, original.data(), cmpSize);
            ASSERT_GT(cmpSize, 0u);

            // A fresh Config, as a reader of the file would have.
            std::vector<float> out(n, 0.f);
            float *op = out.data();
            SZ3::Config readBack;
            SZ_decompress(readBack, cmp, cmpSize, op);
            delete[] cmp;

            double worst = 0;
            for (size_t i = 0; i < n; i++)
                worst = std::max(worst, std::fabs(static_cast<double>(original[i]) - out[i]));
            EXPECT_LE(worst, eb) << nd << "D eb=" << eb;
        }
    }
}

/// The NOTICE records a floor no bound can push past, because a block is coded against one
/// exponent. Pin it both ways: a bound well under the floor is missed, and the error still
/// stops at the fraction of the block maximum the NOTICE quotes.
TEST(SZ3_ZFP, HasAnErrorFloorATighterBoundCannotCross) {
    float peak = 0;
    for (float v : makeField(16, 16, 16, Smooth)) peak = std::max(peak, std::fabs(v));
    const double eb = 1e-10;
    const double worst = roundTrip(16, 16, 16, eb, Smooth);
    EXPECT_GT(worst, eb) << "the floor would have to be gone";
    EXPECT_LE(worst, std::ldexp(1.0, -24) * peak) << "and it sits where the NOTICE says";
}

TEST(SZ3_ZFP, AlgoZfpRejectsFourDimensions) {
    SZ3::Config conf(8, 8, 8, 8);
    conf.cmprAlgo = SZ3::ALGO_ZFP;
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = 1e-3;
    std::vector<float> data(conf.num, 1.f);
    size_t cmpSize = 0;
    EXPECT_THROW(SZ_compress(conf, data.data(), cmpSize), std::invalid_argument);
}
