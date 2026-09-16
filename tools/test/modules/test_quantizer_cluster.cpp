/**
 * Tests for ClusterQuantizer.
 *
 * It is a pure scalar quantizer (it ignores `pred`), so every test drives it
 * directly rather than through a compressor.
 */

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <random>
#include <vector>

#include "SZ3/quantizer/ClusterQuantizer.hpp"
#include "gtest/gtest.h"

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Slack for the comparison itself: the bound is computed with a different rounding
/// path than the quantizer uses, and ties-to-even can land exactly on the bound.

// ---------------------------------------------------------------------------
// GranularBitRoundQuantizer
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// ClusterQuantizer
// ---------------------------------------------------------------------------

/// Build MD-like data: a uniform lattice of levels plus small noise.
static std::vector<float> makeClusteredData(size_t n, float start, float offset, int nlevels, float noise_sigma,
                                            uint32_t seed) {
    std::mt19937 gen(seed);
    std::normal_distribution<float> noise(0.f, noise_sigma);
    std::uniform_int_distribution<int> lvl(0, nlevels - 1);
    std::vector<float> data(n);
    for (auto &v : data) {
        v = start + offset * static_cast<float>(lvl(gen)) + noise(gen);
    }
    return data;
}

TEST(SZ3_ClusterQuantizerTest, ClusterOutRangeStartsAtZero) {
    SZ3::ClusterQuantizer<float> q(10.0f, 2.5f, 12);
    EXPECT_EQ(q.get_out_range().first, 0);
    EXPECT_EQ(q.get_out_range().second, 13);  // bin 0 (unpred) + 12 levels
}

TEST(SZ3_ClusterQuantizerTest, ClusterRejectsBadLevels) {
    EXPECT_THROW(SZ3::ClusterQuantizer<float>(0.f, 0.f, 4), std::invalid_argument);
    EXPECT_THROW(SZ3::ClusterQuantizer<float>(0.f, -1.f, 4), std::invalid_argument);
    EXPECT_THROW(SZ3::ClusterQuantizer<float>(0.f, 1.f, 0), std::invalid_argument);
    EXPECT_THROW(SZ3::ClusterQuantizer<float>(SZ3::ClusterLevels{}), std::invalid_argument);
    EXPECT_NO_THROW(SZ3::ClusterQuantizer<float>(0.f, 1.f, 4));
}

/// End-to-end: derive levels from a sample, quantize, verify the bound and the codebook hit rate.
TEST(SZ3_ClusterQuantizerTest, ClusterDeriveLevelsAndRoundTrip) {
    const float start = 10.0f, offset = 2.5f;
    const int nlevels = 12;
    auto data = makeClusteredData(4000, start, offset, nlevels, 0.05f, 1234u);

    const SZ3::ClusterLevels levels = SZ3::derive_cluster_levels(data.data(), data.size(), 2000);
    ASSERT_TRUE(levels.valid()) << "k-means failed to find the planted lattice";
    EXPECT_NEAR(levels.level_offset, offset, 0.05f);
    EXPECT_GE(levels.level_num, nlevels);

    SZ3::ClusterQuantizer<float> q(levels);
    const double eb = q.get_eb();
    EXPECT_NEAR(eb, 0.5 * levels.level_offset, 1e-6);

    std::vector<int> bins(data.size());
    std::vector<float> recon(data.size());
    for (size_t i = 0; i < data.size(); i++) {
        float scratch = data[i];
        bins[i] = q.quantize_and_overwrite(scratch, 0.f);
        recon[i] = scratch;
        EXPECT_LE(std::fabs(static_cast<double>(scratch) - data[i]), eb) << "i=" << i;
        EXPECT_GE(bins[i], q.get_out_range().first);
        EXPECT_LT(bins[i], q.get_out_range().second);
    }
    // With noise this small essentially everything should land on the codebook.
    EXPECT_LT(q.get_unpred_num(), data.size() / 100) << "too many samples missed the codebook";

    SZ3::ClusterQuantizer<float> q2(levels);
    // Feed the same unpredictable values back for decoding.
    std::vector<unsigned char> buf(64 + q.get_unpred_num() * sizeof(float) + 64);
    unsigned char *sp = buf.data();
    q.save(sp);
    const unsigned char *lp = buf.data();
    size_t remaining = static_cast<size_t>(sp - buf.data());
    q2.load(lp, remaining);
    q2.predecompress_data();
    for (size_t i = 0; i < data.size(); i++) {
        EXPECT_EQ(q2.recover(0.f, bins[i]), recon[i]) << "i=" << i;
    }
}

TEST(SZ3_ClusterQuantizerTest, ClusterDeriveLevelsRejectsUnclusterableData) {
    std::mt19937 gen(5u);
    std::uniform_real_distribution<float> u(0.f, 1.f);
    std::vector<float> uniform(4000);
    for (auto &v : uniform) v = u(gen);
    EXPECT_FALSE(SZ3::derive_cluster_levels(uniform.data(), uniform.size(), 2000).valid());

    std::vector<float> constant(1000, 7.5f);
    EXPECT_FALSE(SZ3::derive_cluster_levels(constant.data(), constant.size()).valid());

    std::vector<float> tiny(3, 1.f);
    EXPECT_FALSE(SZ3::derive_cluster_levels(tiny.data(), tiny.size()).valid());
    EXPECT_FALSE(SZ3::derive_cluster_levels<float>(nullptr, 0).valid());
}

/// `derive_cluster_levels` must also work on double input (KmeansUtil itself only
/// instantiates for float, so the wrapper funnels the sample through float).
TEST(SZ3_ClusterQuantizerTest, ClusterDeriveLevelsDoubleInput) {
    auto f = makeClusteredData(3000, -3.0f, 1.75f, 6, 0.02f, 777u);
    std::vector<double> d(f.begin(), f.end());
    const SZ3::ClusterLevels levels = SZ3::derive_cluster_levels(d.data(), d.size(), 1500);
    ASSERT_TRUE(levels.valid());
    EXPECT_NEAR(levels.level_offset, 1.75f, 0.05f);

    SZ3::ClusterQuantizer<double> q(levels);
    for (double v : d) {
        double scratch = v;
        q.quantize_and_overwrite(scratch, 0.0);
        EXPECT_LE(std::fabs(scratch - v), q.get_eb());
    }
}

/// The derived lattice must span the full data range, not just the sampled subset,
/// and values outside that range must still round-trip exactly via the unpred path.
TEST(SZ3_ClusterQuantizerTest, ClusterOutsideSampledRange) {
    auto data = makeClusteredData(3000, 100.0f, 8.0f, 5, 0.02f, 31337u);
    const SZ3::ClusterLevels levels = SZ3::derive_cluster_levels(data.data(), data.size(), 1500);
    ASSERT_TRUE(levels.valid());

    SZ3::ClusterQuantizer<float> q(levels);
    const double eb = q.get_eb();

    // Values far outside [level(0), level(level_num-1)] -- including zero, negatives and
    // extreme magnitudes -- must be stored verbatim and recovered bit-exactly.
    const std::vector<float> outsiders = {0.f,
                                          -0.f,
                                          -1e6f,
                                          1e6f,
                                          -273.15f,
                                          std::numeric_limits<float>::max(),
                                          std::numeric_limits<float>::lowest(),
                                          1e-30f,
                                          std::numeric_limits<float>::infinity(),
                                          -std::numeric_limits<float>::infinity(),
                                          std::numeric_limits<float>::quiet_NaN()};
    std::vector<int> bins;
    for (float v : outsiders) {
        float scratch = v;
        const int bin = q.quantize_and_overwrite(scratch, 0.f);
        EXPECT_EQ(bin, 0) << "value " << v << " should not have hit the codebook";
        if (std::isnan(v)) {
            EXPECT_TRUE(std::isnan(scratch)) << "data must be left untouched on the unpred path";
        } else {
            EXPECT_EQ(scratch, v) << "data must be left untouched on the unpred path";
        }
        bins.push_back(bin);
    }
    EXPECT_EQ(q.get_unpred_num(), outsiders.size());

    q.predecompress_data();
    for (size_t i = 0; i < outsiders.size(); i++) {
        const float rec = q.recover(0.f, bins[i]);
        if (std::isnan(outsiders[i])) {
            EXPECT_TRUE(std::isnan(rec));
        } else {
            EXPECT_EQ(rec, outsiders[i]) << "i=" << i;
        }
    }

    // A value just past the top level but within eb of it is still allowed on the codebook
    // path -- the guarantee is about the error, not about being inside the range.
    const float top = levels.level_start + levels.level_offset * static_cast<float>(levels.level_num - 1);
    float near_top = top + static_cast<float>(eb) * 0.5f;
    const float near_top_ori = near_top;
    const int bin = q.quantize_and_overwrite(near_top, 0.f);
    EXPECT_GT(bin, 0);
    EXPECT_LE(std::fabs(static_cast<double>(near_top) - near_top_ori), eb);
}

TEST(SZ3_ClusterQuantizerTest, ClusterSaveLoad) {
    const float start = -3.0f, offset = 1.25f;
    const int nlevels = 8;
    SZ3::ClusterQuantizer<float> q(start, offset, nlevels, /*eb=*/0.2);

    std::mt19937 gen(2468u);
    std::uniform_int_distribution<int> lvl(0, nlevels - 1);
    std::normal_distribution<float> noise(0.f, 0.4f);  // large enough to force some unpreds
    const int N = 500;
    std::vector<float> originals(N), recon(N);
    std::vector<int> bins(N);
    for (int i = 0; i < N; i++) {
        originals[i] = start + offset * static_cast<float>(lvl(gen)) + noise(gen);
        float scratch = originals[i];
        bins[i] = q.quantize_and_overwrite(scratch, 0.f);
        recon[i] = scratch;
    }
    ASSERT_GT(q.get_unpred_num(), 0u) << "test needs some unpredictable values to be meaningful";
    ASSERT_LT(q.get_unpred_num(), static_cast<size_t>(N)) << "test needs some codebook hits too";

    std::vector<unsigned char> buffer(1024 + q.get_unpred_num() * sizeof(float));
    unsigned char *sp = buffer.data();
    q.save(sp);
    const size_t saved = static_cast<size_t>(sp - buffer.data());
    EXPECT_GT(saved, 0u);

    SZ3::ClusterQuantizer<float> q2;  // default-constructed: everything comes from the stream
    const unsigned char *lp = buffer.data();
    size_t remaining = saved;
    q2.load(lp, remaining);
    EXPECT_EQ(remaining, 0u) << "load() must consume exactly what save() wrote";
    EXPECT_EQ(q2.get_level_start(), start);
    EXPECT_EQ(q2.get_level_offset(), offset);
    EXPECT_EQ(q2.get_level_num(), nlevels);
    EXPECT_EQ(q2.get_eb(), 0.2);
    EXPECT_EQ(q2.get_unpred_num(), q.get_unpred_num());
    EXPECT_EQ(q2.get_out_range(), q.get_out_range());

    for (int i = 0; i < N; i++) {
        const float rec = q2.recover(0.f, bins[i]);
        EXPECT_EQ(rec, recon[i]) << "i=" << i;
        EXPECT_LE(std::fabs(static_cast<double>(rec) - originals[i]), q2.get_eb()) << "i=" << i;
    }
}

TEST(SZ3_ClusterQuantizerTest, ClusterLoadRejectsForeignStream) {
    SZ3::ClusterQuantizer<float> q;
    unsigned char buffer[64] = {0};
    buffer[0] = 0xAB;  // wrong uid
    const unsigned char *p = buffer;
    size_t remaining = sizeof(buffer);
    EXPECT_THROW(q.load(p, remaining), std::invalid_argument);
}

/// The core claim: whatever the data and whatever the codebook, the reconstruction error
/// never exceeds eb, and the unpredictable path is bit-exact.
TEST(SZ3_ClusterQuantizerTest, ClusterRandomizedGuarantee) {
    std::mt19937 gen(864213u);
    std::uniform_real_distribution<float> start_dist(-1000.f, 1000.f);
    std::uniform_real_distribution<float> offset_dist(1e-3f, 50.f);
    std::uniform_int_distribution<int> nlevel_dist(1, 64);
    std::uniform_int_distribution<int> ebmode(0, 2);

    for (int trial = 0; trial < 200; trial++) {
        const float start = start_dist(gen);
        const float offset = offset_dist(gen);
        const int nlevels = nlevel_dist(gen);
        const int mode = ebmode(gen);
        const double eb_arg = (mode == 0) ? 0.0 : (mode == 1 ? 0.1 * offset : 0.9 * offset);
        SZ3::ClusterQuantizer<float> q(start, offset, nlevels, eb_arg);
        const double eb = q.get_eb();
        EXPECT_GT(eb, 0.0);

        // Data spanning far more than the codebook range, plus specials.
        std::uniform_real_distribution<float> data_dist(start - 5.f * offset * nlevels, start + 5.f * offset * nlevels);
        std::vector<float> originals;
        std::vector<int> bins;
        std::vector<float> recon;
        for (int i = 0; i < 400; i++) {
            float v;
            switch (i % 40) {
                case 0:
                    v = 0.f;
                    break;
                case 1:
                    v = std::numeric_limits<float>::infinity();
                    break;
                case 2:
                    v = -std::numeric_limits<float>::infinity();
                    break;
                case 3:
                    v = std::numeric_limits<float>::quiet_NaN();
                    break;
                case 4:
                    v = std::numeric_limits<float>::max();
                    break;
                case 5:
                    v = std::numeric_limits<float>::lowest();
                    break;
                case 6:
                    v = 1e-38f;
                    break;
                default:
                    v = data_dist(gen);
                    break;
            }
            originals.push_back(v);
            float scratch = v;
            const int bin = q.quantize_and_overwrite(scratch, 12345.f);  // pred ignored
            bins.push_back(bin);
            recon.push_back(scratch);

            EXPECT_GE(bin, q.get_out_range().first) << "trial=" << trial << " i=" << i;
            EXPECT_LT(bin, q.get_out_range().second) << "trial=" << trial << " i=" << i;
            if (bin == 0) {
                // Exact path: value untouched.
                if (std::isnan(v)) {
                    EXPECT_TRUE(std::isnan(scratch));
                } else {
                    EXPECT_EQ(scratch, v) << "trial=" << trial << " i=" << i;
                }
            } else {
                ASSERT_FALSE(std::isnan(v)) << "NaN must never hit the codebook";
                EXPECT_LE(std::fabs(static_cast<double>(scratch) - static_cast<double>(v)), eb)
                    << "trial=" << trial << " i=" << i << " v=" << v;
            }
        }

        // Decode in encode order and confirm the reconstruction matches exactly.
        q.predecompress_data();
        for (size_t i = 0; i < bins.size(); i++) {
            const float rec = q.recover(0.f, bins[i]);
            if (std::isnan(recon[i])) {
                EXPECT_TRUE(std::isnan(rec)) << "trial=" << trial << " i=" << i;
            } else {
                EXPECT_EQ(rec, recon[i]) << "trial=" << trial << " i=" << i;
            }
        }
    }
}

TEST(SZ3_ClusterQuantizerTest, ClusterForceSaveUnpredIsExact) {
    SZ3::ClusterQuantizer<double> q(0.0f, 1.0f, 4);
    const double vals[] = {1.0 / 3.0, -1e300, 1e-300, 0.0};
    std::vector<int> bins;
    for (double v : vals) {
        const int bin = q.force_save_unpred(v);
        EXPECT_EQ(bin, 0);
        bins.push_back(bin);
    }
    q.predecompress_data();
    for (size_t i = 0; i < bins.size(); i++) {
        EXPECT_EQ(q.recover(0.0, bins[i]), vals[i]) << "i=" << i;
    }
}
