#include <cmath>
#include <cstdint>
#include <functional>
#include <vector>

#include "SZ3/decomposition/BlockwiseDecomposition.hpp"
#include "SZ3/decomposition/InterpolationDecomposition.hpp"
#include "SZ3/decomposition/MGARDFusedDecomposition.hpp"
#include "SZ3/decomposition/NoPredictionDecomposition.hpp"
#include "SZ3/decomposition/SPERRFusedDecomposition.hpp"
#include "SZ3/decomposition/SVDDecomposition.hpp"
#include "SZ3/decomposition/SZBioMDDecomposition.hpp"
#include "SZ3/decomposition/SZBioMDXtcDecomposition.hpp"
#include "SZ3/decomposition/TimeSeriesDecomposition.hpp"
#include "SZ3/decomposition/ZFPDecomposition.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "gtest/gtest.h"

// ----- Shared helpers ------------------------------------------------------

// Single roundtrip helper: compress -> save -> load -> decompress -> verify.
// The factory lambda fully owns Decomp construction so each test handles its
// quirks (multiple quantizers, no-arg ctors, etc.) without bespoke boilerplate.
template <typename Decomp, typename T>
static void runRoundtrip(SZ3::Config& conf, const std::vector<T>& original, double tol,
                         std::function<Decomp()> make) {
    std::vector<T> data = original;
    auto decomp = make();
    auto bins = decomp.compress(conf, data.data());
    EXPECT_EQ(bins.size(), conf.num);

    std::vector<unsigned char> buf(1u << 22);  // 4 MiB; covers SVD-class state
    unsigned char* sp = buf.data();
    decomp.save(sp);
    size_t saved = sp - buf.data();

    auto decomp2 = make();
    const unsigned char* lp = buf.data();
    size_t remaining = saved;
    decomp2.load(lp, remaining);

    std::vector<T> dec(conf.num);
    decomp2.decompress(conf, bins, dec.data());

    T max_err = 0;
    for (size_t i = 0; i < conf.num; i++) {
        T e = std::fabs(original[i] - dec[i]);
        if (e > max_err) max_err = e;
    }
    EXPECT_LE(max_err, tol) << "max abs err " << max_err << " exceeds " << tol;
}

// Convenience: convert a (z,y,x) lambda into a flat row-major vector.
template <typename T, typename F>
static std::vector<T> make1D(size_t dim, F f) {
    std::vector<T> v(dim);
    for (size_t i = 0; i < dim; i++) v[i] = static_cast<T>(f(i));
    return v;
}
template <typename T, typename F>
static std::vector<T> make2D(size_t H, size_t W, F f) {
    std::vector<T> v(H * W);
    for (size_t y = 0; y < H; y++)
        for (size_t x = 0; x < W; x++) v[y * W + x] = static_cast<T>(f(y, x));
    return v;
}
template <typename T, typename F>
static std::vector<T> make3D(size_t Z, size_t Y, size_t X, F f) {
    std::vector<T> v(Z * Y * X);
    for (size_t z = 0; z < Z; z++)
        for (size_t y = 0; y < Y; y++)
            for (size_t x = 0; x < X; x++) v[(z * Y + y) * X + x] = static_cast<T>(f(z, y, x));
    return v;
}

// Boilerplate to set ABS-mode error bound on Config.
static void setAbsBound(SZ3::Config& conf, double eb) {
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = eb;
}

// ----- MGARDFusedDecomposition (1D / 2D / 3D) -----------------------------------

TEST(SZ3_DecompositionTest, MGARDFusedDecomposition1D) {
    constexpr SZ3::uint N = 1;
    constexpr size_t dim = 257;
    const double eb = 1e-3;
    SZ3::Config conf(dim);
    setAbsBound(conf, eb);
    auto original = make1D<float>(dim, [](size_t i) { return std::sin(0.05f * i) + 0.3f * std::cos(0.13f * i); });

    using Decomp = SZ3::MGARDFusedDecomposition<float, N>;
    runRoundtrip<Decomp, float>(conf, original, eb, [&] { return Decomp(eb, 32768); });
}

TEST(SZ3_DecompositionTest, MGARDFusedDecomposition2D) {
    constexpr SZ3::uint N = 2;
    constexpr size_t dim = 33;
    const double eb = 1e-3;
    SZ3::Config conf(dim, dim);
    setAbsBound(conf, eb);
    auto original = make2D<float>(dim, dim, [](size_t y, size_t x) {
        return std::sin(0.05f * x) * std::cos(0.03f * y);
    });

    using Decomp = SZ3::MGARDFusedDecomposition<float, N>;
    runRoundtrip<Decomp, float>(conf, original, eb, [&] { return Decomp(eb, 32768); });
}

TEST(SZ3_DecompositionTest, MGARDFusedDecomposition3D) {
    constexpr SZ3::uint N = 3;
    constexpr size_t dim = 17;  // floor(log2(17)) = 4 multigrid levels
    const double eb = 1e-2;
    SZ3::Config conf(dim, dim, dim);
    setAbsBound(conf, eb);
    auto original = make3D<float>(dim, dim, dim, [](size_t z, size_t y, size_t x) {
        return std::sin(0.1f * z) + 0.5f * std::cos(0.07f * y) + 0.25f * std::sin(0.03f * x);
    });

    using Decomp = SZ3::MGARDFusedDecomposition<float, N>;
    runRoundtrip<Decomp, float>(conf, original, eb, [&] { return Decomp(eb, 32768); });
}

// ----- NoPredictionDecomposition -------------------------------------------

TEST(SZ3_DecompositionTest, NoPredictionDecomposition1D) {
    constexpr SZ3::uint N = 1;
    constexpr size_t dim = 100;
    const double eb = 1e-2;
    SZ3::Config conf(dim);
    auto original = make1D<float>(dim, [](size_t i) { return std::sin(0.07f * i); });

    using Quant = SZ3::LinearQuantizer<float>;
    using Decomp = SZ3::NoPredictionDecomposition<float, N, Quant>;
    runRoundtrip<Decomp, float>(conf, original, eb,
                                [&] { return Decomp(conf, Quant(eb, conf.quantbinCnt / 2)); });
}

// ----- InterpolationDecomposition ------------------------------------------

TEST(SZ3_DecompositionTest, InterpolationDecomposition2D) {
    constexpr SZ3::uint N = 2;
    constexpr size_t dim = 33;
    const double eb = 1e-2;
    SZ3::Config conf(dim, dim);
    auto original = make2D<float>(dim, dim, [](size_t y, size_t x) {
        return std::sin(0.1f * x) * std::cos(0.07f * y);
    });

    using Quant = SZ3::LinearQuantizer<float>;
    using Decomp = SZ3::InterpolationDecomposition<float, N, Quant>;
    runRoundtrip<Decomp, float>(conf, original, eb,
                                [&] { return Decomp(conf, Quant(eb, conf.quantbinCnt / 2)); });
}

// ----- BlockwiseDecomposition + LorenzoPredictor (SZ2 classic) -------------

TEST(SZ3_DecompositionTest, BlockwiseDecomposition2D_Lorenzo) {
    constexpr SZ3::uint N = 2;
    constexpr size_t dim = 32;
    const double eb = 1e-2;
    SZ3::Config conf(dim, dim);
    setAbsBound(conf, eb);
    conf.blockSize = 8;
    auto original = make2D<float>(dim, dim, [](size_t y, size_t x) {
        return std::sin(0.05f * x) + std::cos(0.03f * y);
    });

    using Pred = SZ3::LorenzoPredictor<float, N, 1>;
    using Quant = SZ3::LinearQuantizer<float>;
    using Decomp = SZ3::BlockwiseDecomposition<float, N, Pred, Quant>;
    runRoundtrip<Decomp, float>(conf, original, eb,
                                [&] { return Decomp(conf, Pred(eb), Quant(eb, conf.quantbinCnt / 2)); });
}

// ----- TimeSeriesDecomposition (per-timestep error accumulates ~T*eb) ------

TEST(SZ3_DecompositionTest, TimeSeriesDecomposition2D) {
    constexpr SZ3::uint N = 2;
    constexpr size_t T_dim = 4, S_dim = 32;
    const double eb = 1e-2;
    SZ3::Config conf(T_dim, S_dim);
    setAbsBound(conf, eb);
    conf.blockSize = 8;
    auto original = make2D<float>(T_dim, S_dim,
                                  [](size_t t, size_t s) { return std::sin(0.07f * s) + 0.001f * t; });

    using Pred = SZ3::LorenzoPredictor<float, N - 1, 1>;
    using Quant = SZ3::LinearQuantizer<float>;
    using Decomp = SZ3::TimeSeriesDecomposition<float, N, Pred, Quant>;
    // Per-timestep prediction means worst-case error grows ~T_dim * eb.
    runRoundtrip<Decomp, float>(conf, original, eb * T_dim, [&] {
        return Decomp(conf, Pred(eb), Quant(eb, conf.quantbinCnt / 2), nullptr);
    });
}

// ----- SZBioMDDecomposition (1D delta quantization) ------------------------

TEST(SZ3_DecompositionTest, SZBioMDDecomposition1D) {
    constexpr SZ3::uint N = 1;
    constexpr size_t dim = 64;
    const double eb = 1e-2;
    SZ3::Config conf(dim);
    setAbsBound(conf, eb);
    auto original = make1D<float>(dim, [](size_t i) { return 1.0f + 0.001f * i; });

    using Quant = SZ3::LinearQuantizer<float>;
    using Decomp = SZ3::SZBioMDDecomposition<float, N, Quant>;
    runRoundtrip<Decomp, float>(conf, original, eb,
                                [&] { return Decomp(conf, Quant(eb, conf.quantbinCnt / 2)); });
}

// ----- SZBioMDXtcDecomposition (2D single-frame, XTC offset) ---------------

TEST(SZ3_DecompositionTest, SZBioMDXtcDecomposition2D) {
    constexpr SZ3::uint N = 2;
    constexpr size_t H = 4, W = 12;
    const double eb = 1e-2;
    SZ3::Config conf(H, W);
    setAbsBound(conf, eb);
    auto original = make2D<float>(H, W, [](size_t y, size_t x) { return 1.0f + 0.001f * (y + x); });

    using Quant = SZ3::LinearQuantizer<float>;
    using Decomp = SZ3::SZBioMDXtcDecomposition<float, N, Quant>;
    runRoundtrip<Decomp, float>(conf, original, eb,
                                [&] { return Decomp(conf, Quant(eb, SZ3::XTC_radius, false)); });
}

// ----- SVDDecomposition (ST-HOSVD + residual quantization) -----------------

TEST(SZ3_DecompositionTest, SVDDecomposition2D) {
    constexpr SZ3::uint N = 2;
    constexpr size_t dim = 16;
    const double eb = 1e-2;
    SZ3::Config conf(dim, dim);
    setAbsBound(conf, eb);
    conf.svd_energy_threshold = 0.95;
    auto original = make2D<float>(dim, dim, [](size_t y, size_t x) {
        return std::sin(0.07f * x) * std::cos(0.05f * y);
    });

    using Quant = SZ3::LinearQuantizer<float>;
    using Decomp = SZ3::SVDDecomposition<float, N, Quant>;
    runRoundtrip<Decomp, float>(conf, original, eb, [&] {
        return Decomp(conf, Quant(eb * 0.1, conf.quantbinCnt / 2), Quant(eb, conf.quantbinCnt / 2));
    });
}

// ----- SPERRFusedDecomposition (3D-only wavelet+SPECK, PWE mode) ----------------

TEST(SZ3_DecompositionTest, SPERRFusedDecomposition3D) {
    constexpr SZ3::uint N = 3;
    constexpr size_t dim = 8;
    const double eb = 1e-2;
    SZ3::Config conf(dim, dim, dim);
    setAbsBound(conf, eb);
    auto original = make3D<float>(dim, dim, dim, [](size_t z, size_t y, size_t x) {
        return std::sin(0.1f * z) + 0.5f * std::cos(0.07f * y) + 0.25f * std::sin(0.03f * x);
    });

    using Decomp = SZ3::SPERRFusedDecomposition<float, N>;
    runRoundtrip<Decomp, float>(conf, original, eb, [&] { return Decomp(); });
}

// ----- ZFPDecomposition (truly different: empty save/load, no eb param) ----
// ZFP's bin count is 1 + n_blocks*(1 + block_size), not conf.num, and its
// `decomposition.save/load` are no-ops -- so it doesn't fit the shared helper.

TEST(SZ3_DecompositionTest, ZFPDecomposition2D) {
    constexpr SZ3::uint N = 2;
    constexpr size_t dim = 16;
    SZ3::Config conf(dim, dim);
    conf.absErrorBound = 1e-2;
    auto original = make2D<float>(dim, dim, [](size_t y, size_t x) {
        return std::sin(0.05f * x) + std::cos(0.03f * y);
    });

    std::vector<float> data = original;
    auto decomp = SZ3::make_decomposition_zfp<float, int, N>();
    auto bins = decomp.compress(conf, data.data());
    EXPECT_GT(bins.size(), 0u);

    auto decomp2 = SZ3::make_decomposition_zfp<float, int, N>();
    std::vector<float> dec(conf.num, 0.0f);
    decomp2.decompress(conf, bins, dec.data());

    float max_err = 0;
    for (size_t i = 0; i < conf.num; i++) {
        float e = std::fabs(original[i] - dec[i]);
        if (e > max_err) max_err = e;
    }
    EXPECT_LE(max_err, 0.1f) << "ZFP 2D transform roundtrip error too large";
}
