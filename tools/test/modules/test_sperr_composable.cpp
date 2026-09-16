// Tests for the *split* SPERR modules -- the decomposed counterpart of the fused
// SZ3::SPERRFusedDecomposition:
//
//   SPERRTransform      (preprocessor)  conditioner + CDF9/7 wavelet
//   ScalarQuantizer     (quantizer)     mid-tread quantization of coefficients   [pre-existing]
//   OutlierQuantizer    (quantizer)     sparse corrections, LinearQuantizer::unpred style
//   SPERREncoder        (encoder)       SPECK bitstream                          [pre-existing]
//
// The headline result these tests are here to establish is *numerical equivalence*
// with SZ3::SPERRFusedDecomposition: same input, same Config, bit-identical reconstruction.
// The fused module stays the supported path; nothing here modifies it.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <vector>

#include "SZ3/compressor/SZGenericCompressor.hpp"
#include "SZ3/decomposition/SPERRDecomposition.hpp"
#include "SZ3/decomposition/SPERRFusedDecomposition.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/encoder/SPERREncoder.hpp"
#include "SZ3/lossless/Lossless_bypass.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/preprocessor/SPERRTransform.hpp"
#include "SZ3/quantizer/OutlierQuantizer.hpp"
#include "SZ3/quantizer/ScalarQuantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "gtest/gtest.h"

namespace {

constexpr SZ3::uint N3 = 3;

// ----- Fixtures ------------------------------------------------------------

template <typename T, typename F>
std::vector<T> make3D(size_t Z, size_t Y, size_t X, F f) {
    std::vector<T> v(Z * Y * X);
    for (size_t z = 0; z < Z; z++)
        for (size_t y = 0; y < Y; y++)
            for (size_t x = 0; x < X; x++) v[(z * Y + y) * X + x] = static_cast<T>(f(z, y, x));
    return v;
}

SZ3::Config makeConf(size_t dim, double eb) {
    SZ3::Config conf(dim, dim, dim);
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = eb;
    return conf;
}

// Smooth, well-behaved field: no outliers expected at a loose bound.
std::vector<float> smoothField(size_t dim) {
    return make3D<float>(dim, dim, dim, [](size_t z, size_t y, size_t x) {
        return std::sin(0.1 * z) + 0.5 * std::cos(0.07 * y) + 0.25 * std::sin(0.03 * x);
    });
}

// Smooth base plus deterministic sharp spikes -- the wavelet cannot follow these,
// so PWE outlier corrections are guaranteed.
std::vector<float> spikyField(size_t dim) {
    std::vector<float> v = smoothField(dim);
    uint32_t state = 12345u;
    for (size_t i = 0; i < v.size(); i++) {
        state = state * 1664525u + 1013904223u;
        if ((state >> 24) < 24) {  // ~9% of samples
            v[i] += ((state & 1u) ? 1.0f : -1.0f) * (3.0f + static_cast<float>((state >> 8) & 0x3f) * 0.1f);
        }
    }
    return v;
}

SZ3::SPERR::dims_type sperrDims(const SZ3::Config &conf) {
    return SZ3::SPERR::dims_type{conf.dims[2], conf.dims[1], conf.dims[0]};
}

// Compress/save/load/decompress a decomposition, mirroring what SZGenericCompressor
// does around the encoder stage (minus the encoder itself).
template <typename Decomp, typename T>
std::vector<T> decompositionRoundtrip(const SZ3::Config &conf, const std::vector<T> &original,
                                      std::vector<int64_t> *bins_out = nullptr,
                                      std::vector<unsigned char> *state_out = nullptr) {
    std::vector<T> data = original;
    Decomp decomp;
    std::vector<int64_t> bins = decomp.compress(conf, data.data());

    std::vector<unsigned char> buf(1u << 24);
    unsigned char *sp = buf.data();
    decomp.save(sp);
    const size_t saved = static_cast<size_t>(sp - buf.data());

    Decomp decomp2;
    const unsigned char *lp = buf.data();
    size_t remaining = saved;
    decomp2.load(lp, remaining);
    EXPECT_EQ(static_cast<size_t>(lp - buf.data()), saved) << "load() must consume exactly what save() wrote";

    std::vector<T> dec(conf.num, T{0});
    decomp2.decompress(conf, bins, dec.data());

    if (bins_out != nullptr) *bins_out = bins;
    if (state_out != nullptr) state_out->assign(buf.data(), buf.data() + saved);
    return dec;
}

template <typename T>
double maxAbsErr(const std::vector<T> &a, const std::vector<T> &b) {
    double m = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
        m = std::max(m, std::fabs(static_cast<double>(a[i]) - static_cast<double>(b[i])));
    }
    return m;
}

// Full pipeline roundtrip through SZGenericCompressor with a caller-chosen encoder.
template <typename T, typename Encoder>
std::vector<T> pipelineRoundtrip(const SZ3::Config &conf, const std::vector<T> &original, Encoder encoder,
                                 size_t *cmp_size_out = nullptr) {
    std::vector<T> data = original;
    std::vector<unsigned char> cmp(conf.num * sizeof(T) * 4 + (1u << 16));

    auto compressor =
        SZ3::make_compressor_sz_generic<T, N3>(SZ3::SPERRDecomposition<T, N3>(), encoder, SZ3::Lossless_zstd());
    const size_t cmp_size = compressor->compress(conf, data.data(), cmp.data(), cmp.size());
    if (cmp_size_out != nullptr) *cmp_size_out = cmp_size;

    auto decompressor =
        SZ3::make_compressor_sz_generic<T, N3>(SZ3::SPERRDecomposition<T, N3>(), encoder, SZ3::Lossless_zstd());
    std::vector<T> dec(conf.num, T{0});
    decompressor->decompress(conf, cmp.data(), cmp_size, dec.data());
    return dec;
}

}  // namespace

// ----- 1. Transform round-trip --------------------------------------------

TEST(SZ3_SPERRComposable, TransformForwardInverseRoundTrip) {
    constexpr size_t dim = 16;
    SZ3::Config conf = makeConf(dim, 1e-2);
    const std::vector<float> original = smoothField(dim);

    SZ3::SPERRTransform<float, N3> transform(sperrDims(conf));
    std::vector<double> coeffs = transform.preprocess(original.data(), conf.num);
    ASSERT_EQ(coeffs.size(), conf.num);
    EXPECT_FALSE(transform.constant_field());

    std::vector<float> dec(conf.num, 0.0f);
    transform.postprocess(std::move(coeffs), dec.data(), conf.num);

    // Lossless up to double->float rounding; the transform itself adds no quantization.
    EXPECT_LT(maxAbsErr(original, dec), 1e-4) << "conditioner+CDF9/7 roundtrip should be near-lossless";
}

TEST(SZ3_SPERRComposable, TransformStagesMatchFusedTransform) {
    constexpr size_t dim = 16;
    SZ3::Config conf = makeConf(dim, 1e-2);
    const std::vector<float> original = smoothField(dim);

    // Split: condition, then wavelet, as two separate calls.
    SZ3::SPERRTransform<float, N3> transform(sperrDims(conf));
    std::vector<double> conditioned = transform.condition(original.data(), conf.num);
    std::vector<double> staged_coeffs = transform.forward_wavelet(conditioned);

    // Fused: SPERRFusedDecomposition::forward() does both at once.
    SZ3::SPERRFusedDecomposition<float, N3> fused;
    std::vector<float> data = original;
    const SZ3::SPERR3DLayout layout = fused.prepare(conf);
    const SZ3::SPERRFrame frame = fused.forward(conf, layout, data.data());

    ASSERT_EQ(staged_coeffs.size(), frame.wavelet_coeffs.size());
    for (size_t i = 0; i < staged_coeffs.size(); i++) {
        ASSERT_EQ(staged_coeffs[i], frame.wavelet_coeffs[i]) << "wavelet coefficient " << i << " differs";
    }
    EXPECT_EQ(transform.header(), frame.conditioner_header) << "conditioner headers must be identical";
    EXPECT_EQ(transform.constant_field(), frame.constant_field);
    ASSERT_EQ(conditioned.size(), frame.conditioned_values.size());
    for (size_t i = 0; i < conditioned.size(); i++) {
        ASSERT_EQ(conditioned[i], frame.conditioned_values[i]) << "conditioned value " << i << " differs";
    }
}

TEST(SZ3_SPERRComposable, TransformStateSaveLoad) {
    constexpr size_t dim = 8;
    SZ3::Config conf = makeConf(dim, 1e-2);
    const std::vector<float> original = smoothField(dim);

    SZ3::SPERRTransform<float, N3> transform(sperrDims(conf));
    (void)transform.preprocess(original.data(), conf.num);

    std::vector<unsigned char> buf(256);
    unsigned char *sp = buf.data();
    transform.save(sp);
    const size_t saved = static_cast<size_t>(sp - buf.data());
    EXPECT_EQ(saved, transform.size_est());

    SZ3::SPERRTransform<float, N3> loaded(sperrDims(conf));
    const unsigned char *lp = buf.data();
    size_t remaining = saved;
    loaded.load(lp, remaining);
    EXPECT_EQ(remaining, 0u);
    EXPECT_EQ(loaded.header(), transform.header());
    EXPECT_EQ(loaded.constant_field(), transform.constant_field());
}

// ----- 2. OutlierQuantizer as a standalone module ---------------------------

TEST(SZ3_SPERRComposable, OutlierQuantizerDeadzoneAndCorrection) {
    const double tol = 0.1;
    SZ3::OutlierQuantizer<double, int64_t> oq(tol);

    // Inside the tolerance band -> bin 0, value snaps to the approximation.
    double v = 1.05;
    EXPECT_EQ(oq.quantize_and_overwrite(v, 1.0), 0);
    EXPECT_DOUBLE_EQ(v, 1.0);

    // Outside -> non-zero bin, and the corrected value lands within tolerance.
    double w = 1.6;
    const int64_t bin = oq.quantize_and_overwrite(w, 1.0);
    EXPECT_NE(bin, 0);
    EXPECT_LE(std::fabs(w - 1.6), tol);
    EXPECT_DOUBLE_EQ(w, oq.recover(1.0, bin));
}

TEST(SZ3_SPERRComposable, OutlierQuantizerCollectApplySaveLoad) {
    const double tol = 0.05;
    std::vector<double> truth(200);
    std::vector<double> approx(200);
    for (size_t i = 0; i < truth.size(); i++) {
        truth[i] = std::sin(0.05 * static_cast<double>(i));
        approx[i] = truth[i] + ((i % 17 == 0) ? 0.9 : 0.01);  // every 17th sample is an outlier
    }

    SZ3::OutlierQuantizer<double, int64_t> oq(tol);
    const size_t collected = oq.collect(truth, approx);
    EXPECT_EQ(collected, oq.size());
    EXPECT_EQ(collected, (truth.size() + 16) / 17);

    std::vector<unsigned char> buf(1u << 16);
    unsigned char *sp = buf.data();
    oq.save(sp);
    const size_t saved = static_cast<size_t>(sp - buf.data());
    EXPECT_EQ(saved, oq.size_est());

    SZ3::OutlierQuantizer<double, int64_t> loaded(1.0);
    const unsigned char *lp = buf.data();
    size_t remaining = saved;
    loaded.load(lp, remaining);
    EXPECT_EQ(remaining, 0u);
    EXPECT_DOUBLE_EQ(loaded.get_tolerance(), tol);
    ASSERT_EQ(loaded.size(), oq.size());

    std::vector<double> corrected = approx;
    loaded.apply(corrected);
    for (size_t i = 0; i < truth.size(); i++) {
        EXPECT_LE(std::fabs(truth[i] - corrected[i]), tol) << "index " << i;
    }

    // Dense <-> sparse interop (the shape SPERRFusedDecomposition ships outliers in).
    const std::vector<int64_t> dense = oq.to_dense(truth.size());
    SZ3::OutlierQuantizer<double, int64_t> from_dense(tol);
    from_dense.from_dense(dense);
    ASSERT_EQ(from_dense.size(), oq.size());
    for (size_t i = 0; i < oq.size(); i++) {
        EXPECT_EQ(from_dense.get_corrections()[i].pos, oq.get_corrections()[i].pos);
        EXPECT_EQ(from_dense.get_corrections()[i].bin, oq.get_corrections()[i].bin);
    }
}

// ----- 3. Four pieces composed vs. the fused module -------------------------

TEST(SZ3_SPERRComposable, SplitDecompositionMatchesFusedBitForBit) {
    constexpr size_t dim = 16;
    const double eb = 1e-2;
    SZ3::Config conf = makeConf(dim, eb);
    const std::vector<float> original = smoothField(dim);

    const std::vector<float> fused =
        decompositionRoundtrip<SZ3::SPERRFusedDecomposition<float, N3>, float>(conf, original);
    const std::vector<float> split = decompositionRoundtrip<SZ3::SPERRDecomposition<float, N3>, float>(conf, original);

    ASSERT_EQ(fused.size(), split.size());
    for (size_t i = 0; i < fused.size(); i++) {
        ASSERT_EQ(fused[i], split[i]) << "reconstruction differs at " << i;
    }
    EXPECT_LE(maxAbsErr(original, split), eb);
}

TEST(SZ3_SPERRComposable, SplitInternalsMatchFusedInternals) {
    constexpr size_t dim = 16;
    const double eb = 1e-2;
    SZ3::Config conf = makeConf(dim, eb);
    const std::vector<float> original = spikyField(dim);

    // Fused pipeline, driven stage by stage through its public API.
    SZ3::SPERRFusedDecomposition<float, N3> fused;
    std::vector<float> fused_data = original;
    const SZ3::SPERR3DLayout layout = fused.prepare(conf);
    SZ3::SPERRFrame frame = fused.forward(conf, layout, fused_data.data());
    fused.quantize_and_collect_outliers(frame);

    // Split pipeline.
    SZ3::SPERRDecomposition<float, N3> split;
    std::vector<float> split_data = original;
    const std::vector<int64_t> split_bins = split.compress(conf, split_data.data());

    // (a) quantization step chosen by the rate-control policy
    EXPECT_DOUBLE_EQ(split.get_q(), frame.q);

    // (b) quantized wavelet coefficients -- what the encoder stage sees
    ASSERT_EQ(split_bins.size(), frame.coeff_bins.size());
    for (size_t i = 0; i < split_bins.size(); i++) {
        ASSERT_EQ(split_bins[i], frame.coeff_bins[i]) << "coefficient bin " << i << " differs";
    }

    // (c) outliers: same positions, same quantized corrections
    const auto &corrections = split.get_outlier_quantizer().get_corrections();
    ASSERT_GT(frame.outliers.size(), 0u) << "the spiky fixture must actually produce outliers";
    ASSERT_EQ(corrections.size(), frame.outliers.size());
    SZ3::ScalarQuantizer<double, int64_t> reference(eb, 1.1, -0.25);
    for (size_t i = 0; i < corrections.size(); i++) {
        EXPECT_EQ(corrections[i].pos, frame.outliers[i].pos) << "outlier " << i << " position differs";
        double err = frame.outliers[i].err;
        EXPECT_EQ(corrections[i].bin, reference.quantize_and_overwrite(err, 0.0))
            << "outlier " << i << " correction bin differs";
    }
}

TEST(SZ3_SPERRComposable, SplitCoefficientStreamMatchesFusedSpeckPayload) {
    constexpr size_t dim = 16;
    const double eb = 1e-2;
    SZ3::Config conf = makeConf(dim, eb);
    const std::vector<float> original = spikyField(dim);

    // Fused: the SPECK stream is stashed inside the decomposition and emitted by save().
    // Layout of save(): [size_t stream_len][17B conditioner header][SPECK payload]
    std::vector<unsigned char> fused_state;
    (void)decompositionRoundtrip<SZ3::SPERRFusedDecomposition<float, N3>, float>(conf, original, nullptr, &fused_state);
    ASSERT_GT(fused_state.size(), sizeof(size_t) + 17u);
    size_t stream_len = 0;
    std::memcpy(&stream_len, fused_state.data(), sizeof(size_t));
    ASSERT_EQ(stream_len + sizeof(size_t), fused_state.size());
    const unsigned char *fused_header = fused_state.data() + sizeof(size_t);
    const unsigned char *fused_speck = fused_header + 17;
    const size_t fused_speck_len = stream_len - 17;
    ASSERT_GT(fused_speck_len, 512u) << "the comparison must be against a substantive payload";

    // Split: the same SPECK payload, but produced by the encoder module outside the
    // decomposition, from the bins the decomposition returned.
    SZ3::SPERRDecomposition<float, N3> split;
    std::vector<float> split_data = original;
    const std::vector<int64_t> split_bins = split.compress(conf, split_data.data());
    SZ3::SPERREncoder<int64_t, N3> encoder(sperrDims(conf));
    const std::vector<unsigned char> split_speck = encoder.encode_stream(split_bins);

    ASSERT_EQ(split_speck.size(), fused_speck_len) << "SPECK payload length differs";
    EXPECT_EQ(0, std::memcmp(split_speck.data(), fused_speck, fused_speck_len)) << "SPECK payload differs";

    // The conditioner header agrees on meta byte + mean; bytes 9..16 hold `q` in the
    // fused module (Conditioner::save_q), whereas the split path keeps `q` in the
    // ScalarQuantizer's own serialized state. Both must name the same step.
    EXPECT_EQ(0, std::memcmp(split.get_transform().header().data(), fused_header, 9u))
        << "conditioner meta/mean differ";
    const SZ3::SPERR::Conditioner conditioner;
    SZ3::SPERR::condi_type fused_condi{};
    std::memcpy(fused_condi.data(), fused_header, fused_condi.size());
    EXPECT_DOUBLE_EQ(conditioner.retrieve_q(fused_condi), split.get_q());
}

// ----- 4. Error bound ------------------------------------------------------

TEST(SZ3_SPERRComposable, ErrorBoundHoldsAcrossBounds) {
    constexpr size_t dim = 16;
    const std::vector<float> original = spikyField(dim);
    const double bounds[] = {1e-1, 1e-2, 1e-3};

    for (double eb : bounds) {
        SZ3::Config conf = makeConf(dim, eb);
        const std::vector<float> split =
            decompositionRoundtrip<SZ3::SPERRDecomposition<float, N3>, float>(conf, original);
        const std::vector<float> fused =
            decompositionRoundtrip<SZ3::SPERRFusedDecomposition<float, N3>, float>(conf, original);
        EXPECT_LE(maxAbsErr(original, split), eb) << "eb = " << eb;
        for (size_t i = 0; i < split.size(); i++) {
            ASSERT_EQ(split[i], fused[i]) << "eb = " << eb << ", index " << i;
        }
    }
}

// ----- 5. Constant field ---------------------------------------------------

TEST(SZ3_SPERRComposable, ConstantFieldMatchesFused) {
    constexpr size_t dim = 16;
    const double eb = 1e-3;
    SZ3::Config conf = makeConf(dim, eb);
    const std::vector<float> original(conf.num, 3.5f);

    SZ3::SPERRDecomposition<float, N3> probe;
    std::vector<float> probe_data = original;
    const std::vector<int64_t> probe_bins = probe.compress(conf, probe_data.data());
    EXPECT_TRUE(probe.get_transform().constant_field());
    EXPECT_EQ(probe_bins.size(), conf.num) << "a constant field still emits a full-length (zero) bin stream";
    EXPECT_TRUE(std::all_of(probe_bins.begin(), probe_bins.end(), [](int64_t b) { return b == 0; }));
    EXPECT_EQ(probe.get_outlier_quantizer().size(), 0u) << "a constant field has no outliers";

    const std::vector<float> split = decompositionRoundtrip<SZ3::SPERRDecomposition<float, N3>, float>(conf, original);
    const std::vector<float> fused =
        decompositionRoundtrip<SZ3::SPERRFusedDecomposition<float, N3>, float>(conf, original);
    for (size_t i = 0; i < split.size(); i++) {
        ASSERT_EQ(split[i], 3.5f) << "constant field must reconstruct exactly, index " << i;
        ASSERT_EQ(split[i], fused[i]) << "index " << i;
    }
}

TEST(SZ3_SPERRComposable, ConstantFieldThroughFullPipeline) {
    constexpr size_t dim = 8;
    SZ3::Config conf = makeConf(dim, 1e-3);
    const std::vector<float> original(conf.num, -7.25f);

    const std::vector<float> dec =
        pipelineRoundtrip<float>(conf, original, SZ3::SPERREncoder<int64_t, N3>(sperrDims(conf)));
    for (size_t i = 0; i < dec.size(); i++) {
        ASSERT_EQ(dec[i], -7.25f) << "index " << i;
    }
}

// ----- 6. Genuine outliers -------------------------------------------------

TEST(SZ3_SPERRComposable, FieldWithOutliersIsCorrected) {
    constexpr size_t dim = 16;
    const double eb = 1e-3;
    SZ3::Config conf = makeConf(dim, eb);
    const std::vector<float> original = spikyField(dim);

    SZ3::SPERRDecomposition<float, N3> split;
    std::vector<float> data = original;
    (void)split.compress(conf, data.data());
    const size_t n_outliers = split.get_outlier_quantizer().size();
    EXPECT_GT(n_outliers, 0u) << "expected genuine outliers at eb = " << eb;
    EXPECT_LT(n_outliers, conf.num) << "outliers should stay sparse";

    const std::vector<float> dec = decompositionRoundtrip<SZ3::SPERRDecomposition<float, N3>, float>(conf, original);
    EXPECT_LE(maxAbsErr(original, dec), eb);

    // Without the outlier stage the bound would be violated -- confirm the stage is
    // load-bearing rather than decorative.
    SZ3::SPERRDecomposition<float, N3> no_outliers;
    std::vector<float> data2 = original;
    std::vector<int64_t> bins = no_outliers.compress(conf, data2.data());
    no_outliers.get_outlier_quantizer().clear();
    std::vector<float> uncorrected(conf.num, 0.0f);
    no_outliers.decompress(conf, bins, uncorrected.data());
    EXPECT_GT(maxAbsErr(original, uncorrected), eb) << "outlier corrections should be what enforces the bound";
}

// ----- 7. Composition through SZGenericCompressor, SPECK and non-SPECK ------

TEST(SZ3_SPERRComposable, PipelineWithSpeckEncoder) {
    constexpr size_t dim = 16;
    const double eb = 1e-2;
    SZ3::Config conf = makeConf(dim, eb);
    const std::vector<float> original = spikyField(dim);

    size_t cmp_size = 0;
    const std::vector<float> dec =
        pipelineRoundtrip<float>(conf, original, SZ3::SPERREncoder<int64_t, N3>(sperrDims(conf)), &cmp_size);
    EXPECT_GT(cmp_size, 0u);
    EXPECT_LT(cmp_size, conf.num * sizeof(float)) << "pipeline should actually compress";
    EXPECT_LE(maxAbsErr(original, dec), eb);

    // Same reconstruction as the fused module.
    const std::vector<float> fused =
        decompositionRoundtrip<SZ3::SPERRFusedDecomposition<float, N3>, float>(conf, original);
    for (size_t i = 0; i < dec.size(); i++) {
        ASSERT_EQ(dec[i], fused[i]) << "index " << i;
    }
}

// The point of splitting the transform out: it composes with encoders that have
// nothing to do with SPECK.
TEST(SZ3_SPERRComposable, PipelineWithHuffmanEncoder) {
    constexpr size_t dim = 16;
    const double eb = 1e-2;
    SZ3::Config conf = makeConf(dim, eb);
    const std::vector<float> original = spikyField(dim);

    size_t huff_size = 0;
    const std::vector<float> huff =
        pipelineRoundtrip<float>(conf, original, SZ3::HuffmanEncoder<int64_t>(), &huff_size);
    EXPECT_GT(huff_size, 0u);
    EXPECT_LE(maxAbsErr(original, huff), eb);

    size_t speck_size = 0;
    const std::vector<float> speck =
        pipelineRoundtrip<float>(conf, original, SZ3::SPERREncoder<int64_t, N3>(sperrDims(conf)), &speck_size);

    // The encoder is lossless, so swapping it must not move a single reconstructed value.
    ASSERT_EQ(huff.size(), speck.size());
    for (size_t i = 0; i < huff.size(); i++) {
        ASSERT_EQ(huff[i], speck[i]) << "encoder choice must not change the reconstruction, index " << i;
    }
    RecordProperty("speck_bytes", static_cast<int>(speck_size));
    RecordProperty("huffman_bytes", static_cast<int>(huff_size));
}

// Rate comparison against the fused ALGO_SPERR wiring (same encoder, same lossless
// stage). The two ship the same information through different slots: the fused module
// hides the SPECK stream in decomposition state and sends the dense outlier map to the
// encoder, the split module does the opposite.
TEST(SZ3_SPERRComposable, PipelineRateComparableToFusedWiring) {
    constexpr size_t dim = 16;
    const double eb = 1e-2;
    SZ3::Config conf = makeConf(dim, eb);
    const std::vector<float> original = spikyField(dim);

    std::vector<unsigned char> cmp_fused(conf.num * sizeof(float) * 4 + (1u << 16));
    std::vector<float> fused_data = original;
    auto fused_pipeline = SZ3::make_compressor_sz_generic<float, N3>(SZ3::SPERRFusedDecomposition<float, N3>(),
                                                                     SZ3::SPERREncoder<int64_t, N3>(sperrDims(conf)),
                                                                     SZ3::Lossless_bypass());
    const size_t fused_size = fused_pipeline->compress(conf, fused_data.data(), cmp_fused.data(), cmp_fused.size());

    std::vector<unsigned char> cmp_split(conf.num * sizeof(float) * 4 + (1u << 16));
    std::vector<float> split_data = original;
    auto split_pipeline = SZ3::make_compressor_sz_generic<float, N3>(
        SZ3::SPERRDecomposition<float, N3>(), SZ3::SPERREncoder<int64_t, N3>(sperrDims(conf)), SZ3::Lossless_bypass());
    const size_t split_size = split_pipeline->compress(conf, split_data.data(), cmp_split.data(), cmp_split.size());

    auto fused_dec_pipeline = SZ3::make_compressor_sz_generic<float, N3>(
        SZ3::SPERRFusedDecomposition<float, N3>(), SZ3::SPERREncoder<int64_t, N3>(sperrDims(conf)),
        SZ3::Lossless_bypass());
    std::vector<float> fused_dec(conf.num, 0.0f);
    fused_dec_pipeline->decompress(conf, cmp_fused.data(), fused_size, fused_dec.data());

    auto split_dec_pipeline = SZ3::make_compressor_sz_generic<float, N3>(
        SZ3::SPERRDecomposition<float, N3>(), SZ3::SPERREncoder<int64_t, N3>(sperrDims(conf)), SZ3::Lossless_bypass());
    std::vector<float> split_dec(conf.num, 0.0f);
    split_dec_pipeline->decompress(conf, cmp_split.data(), split_size, split_dec.data());

    for (size_t i = 0; i < split_dec.size(); i++) {
        ASSERT_EQ(split_dec[i], fused_dec[i]) << "index " << i;
    }
    EXPECT_LE(maxAbsErr(original, split_dec), eb);
    EXPECT_LT(split_size, 2 * fused_size) << "split wiring should stay in the same ballpark";

    // Same comparison with Zstd, which is what actually squeezes the sparse outlier
    // list the split module keeps in its own serialized state.
    std::vector<unsigned char> cmp_fz(conf.num * sizeof(float) * 4 + (1u << 16));
    std::vector<float> fz_data = original;
    auto fused_zstd = SZ3::make_compressor_sz_generic<float, N3>(SZ3::SPERRFusedDecomposition<float, N3>(),
                                                                 SZ3::SPERREncoder<int64_t, N3>(sperrDims(conf)),
                                                                 SZ3::Lossless_zstd());
    const size_t fused_zstd_size = fused_zstd->compress(conf, fz_data.data(), cmp_fz.data(), cmp_fz.size());

    size_t split_zstd_size = 0;
    (void)pipelineRoundtrip<float>(conf, original, SZ3::SPERREncoder<int64_t, N3>(sperrDims(conf)), &split_zstd_size);

    SZ3::SPERRDecomposition<float, N3> counter;
    std::vector<float> counter_data = original;
    (void)counter.compress(conf, counter_data.data());

    RecordProperty("fused_wiring_bytes_bypass", static_cast<int>(fused_size));
    RecordProperty("split_wiring_bytes_bypass", static_cast<int>(split_size));
    RecordProperty("fused_wiring_bytes_zstd", static_cast<int>(fused_zstd_size));
    RecordProperty("split_wiring_bytes_zstd", static_cast<int>(split_zstd_size));
    RecordProperty("outlier_count", static_cast<int>(counter.get_outlier_quantizer().size()));
}

// ----- 8. PSNR mode --------------------------------------------------------

TEST(SZ3_SPERRComposable, PsnrModeMatchesFused) {
    constexpr size_t dim = 16;
    SZ3::Config conf(dim, dim, dim);
    conf.errorBoundMode = SZ3::EB_PSNR;
    conf.psnrErrorBound = 60.0;
    const std::vector<float> original = smoothField(dim);

    const std::vector<float> split = decompositionRoundtrip<SZ3::SPERRDecomposition<float, N3>, float>(conf, original);
    const std::vector<float> fused =
        decompositionRoundtrip<SZ3::SPERRFusedDecomposition<float, N3>, float>(conf, original);

    ASSERT_EQ(split.size(), fused.size());
    for (size_t i = 0; i < split.size(); i++) {
        ASSERT_EQ(split[i], fused[i]) << "PSNR-mode reconstruction differs at " << i;
    }

    SZ3::SPERRDecomposition<float, N3> probe;
    std::vector<float> probe_data = original;
    (void)probe.compress(conf, probe_data.data());
    EXPECT_EQ(probe.get_outlier_quantizer().size(), 0u) << "PSNR mode carries no outlier corrections";
}

// ----- 9. Guard rails ------------------------------------------------------

TEST(SZ3_SPERRComposable, RejectsNon3DConfiguration) {
    SZ3::Config conf2d(16, 16);
    conf2d.errorBoundMode = SZ3::EB_ABS;
    conf2d.absErrorBound = 1e-2;
    std::vector<float> data(conf2d.num, 1.0f);

    SZ3::SPERRDecomposition<float, 2> split;
    EXPECT_THROW((void)split.compress(conf2d, data.data()), std::invalid_argument);
}

TEST(SZ3_SPERRComposable, OutRangeStartsAtZero) {
    SZ3::SPERRDecomposition<float, N3> split;
    EXPECT_EQ(split.get_out_range().first, 0);
}
