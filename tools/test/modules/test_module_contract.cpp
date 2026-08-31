// Group-level acceptance checks. Every module in the library must pass the contract for its
// group; see docs/MODULE_ACCEPTANCE.md.

#include "ModuleContract.hpp"
#include "SZ3/decomposition/BlockwiseDecomposition.hpp"
#include "SZ3/decomposition/InterpolationDecomposition.hpp"
#include "SZ3/decomposition/MGARDFusedDecomposition.hpp"
#include "SZ3/decomposition/MultiLevelDecomposition.hpp"
#include "SZ3/decomposition/NoPredictionDecomposition.hpp"
#include "SZ3/encoder/ArithmeticEncoder.hpp"
#include "SZ3/encoder/BitshuffleEncoder.hpp"
#include "SZ3/encoder/BypassEncoder.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/encoder/HuffmanEncoderV2.hpp"
#include "SZ3/encoder/RunlengthEncoder.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/quantizer/BitTruncationQuantizer.hpp"
#include "SZ3/quantizer/ClusterQuantizer.hpp"
#include "SZ3/quantizer/FixedPointQuantizer.hpp"
#include "SZ3/quantizer/GranularBitRoundQuantizer.hpp"
#include "SZ3/quantizer/LevelQuantizer.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
#include "SZ3/quantizer/LogDomainQuantizer.hpp"
#include "gtest/gtest.h"

// ----- Quantizer group -----------------------------------------------------

TEST(SZ3_ModuleContract, LinearQuantizer) {
    auto make = [] { return SZ3::LinearQuantizer<float>(1e-3); };
    SZ3_test::expectQuantizerContract<float>("LinearQuantizer", make);
    SZ3_test::expectAbsErrorBound<float>("LinearQuantizer", make, 1e-3);
}

TEST(SZ3_ModuleContract, LinearQuantizerDouble) {
    auto make = [] { return SZ3::LinearQuantizer<double>(1e-6); };
    SZ3_test::expectQuantizerContract<double>("LinearQuantizer<double>", make);
    SZ3_test::expectAbsErrorBound<double>("LinearQuantizer<double>", make, 1e-6);
}

TEST(SZ3_ModuleContract, LevelQuantizerQuadratic) {
    auto make = [] { return SZ3::LevelQuantizer<float>(1e-3, 32768, SZ3::LevelCurve::Quadratic); };
    SZ3_test::expectQuantizerContract<float>("LevelQuantizer/Quadratic", make);
    SZ3_test::expectAbsErrorBound<float>("LevelQuantizer/Quadratic", make, 1e-3);
}

TEST(SZ3_ModuleContract, LevelQuantizerLog) {
    auto make = [] { return SZ3::LevelQuantizer<float>(1e-3, 32768, SZ3::LevelCurve::Log); };
    SZ3_test::expectQuantizerContract<float>("LevelQuantizer/Log", make);
    SZ3_test::expectAbsErrorBound<float>("LevelQuantizer/Log", make, 1e-3);
}

TEST(SZ3_ModuleContract, LogDomainQuantizer) {
    SZ3_test::expectQuantizerContract<float>("LogDomainQuantizer", [] { return SZ3::LogDomainQuantizer<float>(1e-3); });
}

TEST(SZ3_ModuleContract, ClusterQuantizer) {
    auto make = [] { return SZ3::ClusterQuantizer<float>(0.0f, 1.0f, 8, 0.5); };
    SZ3_test::expectQuantizerContract<float>("ClusterQuantizer", make);
    SZ3_test::expectAbsErrorBound<float>("ClusterQuantizer", make, 0.5);
}

TEST(SZ3_ModuleContract, BitTruncationQuantizerFloat) {
    SZ3_test::expectQuantizerContract<float>("BitTruncationQuantizer<float>",
                                             [] { return SZ3::BitTruncationQuantizer<float>(2); });
}

TEST(SZ3_ModuleContract, BitTruncationQuantizerDouble) {
    SZ3_test::expectQuantizerContract<double>("BitTruncationQuantizer<double>",
                                              [] { return SZ3::BitTruncationQuantizer<double>(4); });
}

TEST(SZ3_ModuleContract, GranularBitRoundQuantizer) {
    SZ3_test::expectQuantizerContract<float>("GranularBitRoundQuantizer",
                                             [] { return SZ3::GranularBitRoundQuantizer<float>(4); });
}

TEST(SZ3_ModuleContract, FixedPointQuantizer) {
    // Needs its range before use, so the factory hands back a calibrated instance.
    SZ3_test::expectQuantizerContract<float>("FixedPointQuantizer", [] {
        SZ3::FixedPointQuantizer<float> q(16);
        q.calibrate(1.0f);
        return q;
    });
}

// ----- Encoder group ------------------------------------------------------
//
// `stateNum` is what the compressor would pass: `get_out_range().second`, or 0 where the
// decomposition has no bin range. Encoders that cannot derive their own range need a real one.

TEST(SZ3_ModuleContract, HuffmanEncoder) {
    SZ3_test::expectEncoderContract<int>("HuffmanEncoder", [] { return SZ3::HuffmanEncoder<int>(); });
}

TEST(SZ3_ModuleContract, HuffmanEncoderV2) {
    SZ3_test::expectEncoderContract<int>("HuffmanEncoderV2", [] { return SZ3::HuffmanEncoderV2<int>(); });
}

TEST(SZ3_ModuleContract, BypassEncoder) {
    SZ3_test::expectEncoderContract<int>("BypassEncoder", [] { return SZ3::BypassEncoder<int>(); });
}

TEST(SZ3_ModuleContract, RunlengthEncoder) {
    SZ3_test::expectEncoderContract<int>("RunlengthEncoder", [] { return SZ3::RunlengthEncoder<int>(); });
}

TEST(SZ3_ModuleContract, BitshuffleEncoder) {
    SZ3_test::expectEncoderContract<int>("BitshuffleEncoder", [] { return SZ3::BitshuffleEncoder<int>(32); });
}

TEST(SZ3_ModuleContract, ArithmeticEncoder) {
    // ArithmeticEncoder sizes its frequency model from `stateNum` and cannot derive it, so it is
    // exercised with a real range rather than the 0 the other encoders tolerate.
    SZ3_test::expectEncoderContract<int>(
        "ArithmeticEncoder", [] { return SZ3::ArithmeticEncoder<int>(); }, {}, /*stateNum=*/128);
}

// ----- Decomposition group ------------------------------------------------

namespace {

SZ3::Config conf3D(size_t d0, size_t d1, size_t d2, double eb) {
    SZ3::Config c(d0, d1, d2);
    c.errorBoundMode = SZ3::EB_ABS;
    c.absErrorBound = eb;
    return c;
}

}  // namespace

TEST(SZ3_ModuleContract, NoPredictionDecomposition) {
    const double eb = 1e-2;
    auto conf = conf3D(8, 12, 12, eb);
    SZ3_test::expectDecompositionContract<float>(
        "NoPredictionDecomposition", conf,
        [&] {
            return SZ3::NoPredictionDecomposition<float, 3, SZ3::LinearQuantizer<float>>(
                conf, SZ3::LinearQuantizer<float>(eb));
        },
        eb);
}

TEST(SZ3_ModuleContract, InterpolationDecomposition) {
    const double eb = 1e-2;
    auto conf = conf3D(9, 13, 13, eb);
    SZ3_test::expectDecompositionContract<float>(
        "InterpolationDecomposition", conf,
        [&] {
            return SZ3::InterpolationDecomposition<float, 3, SZ3::LinearQuantizer<float>>(
                conf, SZ3::LinearQuantizer<float>(eb));
        },
        eb);
}

TEST(SZ3_ModuleContract, BlockwiseDecompositionLorenzo) {
    const double eb = 1e-2;
    auto conf = conf3D(8, 12, 12, eb);
    SZ3_test::expectDecompositionContract<float>(
        "BlockwiseDecomposition/Lorenzo", conf,
        [&] {
            return SZ3::BlockwiseDecomposition<float, 3, SZ3::LorenzoPredictor<float, 3, 1>,
                                               SZ3::LinearQuantizer<float>>(
                conf, SZ3::LorenzoPredictor<float, 3, 1>(eb), SZ3::LinearQuantizer<float>(eb));
        },
        eb);
}

// MGARD quantizes multigrid coefficients per level and never checks the reconstruction against
// the bound, so the transform's own round-off passes through. That round-off scales with the
// field's magnitude, and the bound only holds while `max|x| * epsilon<T> << eb`: for float and
// eb = 1e-2, a magnitude of 1e6 already reaches 1.17x eb and 1e7 reaches 6.25x. The fields below
// stay inside the regime the module honours; see MGARDFusedDecomposition.hpp.
static std::vector<std::vector<float>> mgardFields(size_t num) {
    std::vector<std::vector<float>> fields;
    fields.push_back(std::vector<float>(num, 3.5f));
    std::vector<float> ramp(num), smooth(num);
    for (size_t i = 0; i < num; i++) {
        ramp[i] = static_cast<float>(i) * 1e-3f;
        smooth[i] = static_cast<float>(std::sin(static_cast<double>(i) * 0.013) * 3.0 +
                                       std::cos(static_cast<double>(i) * 0.007));
    }
    fields.push_back(ramp);
    fields.push_back(smooth);
    return fields;
}

TEST(SZ3_ModuleContract, MGARDFusedDecomposition) {
    const double eb = 1e-2;
    auto conf = conf3D(9, 9, 9, eb);
    SZ3_test::expectDecompositionContract<float>(
        "MGARDFusedDecomposition", conf, [&] { return SZ3::MGARDFusedDecomposition<float, 3>(eb); }, eb,
        mgardFields(conf.num));
}

TEST(SZ3_ModuleContract, MGARDDecompositionComposable) {
    const double eb = 1e-2;
    const int radius = 32768;
    auto conf = conf3D(9, 9, 9, eb);
    SZ3_test::expectDecompositionContract<float>(
        "MGARDDecomposition", conf,
        [&] {
            return SZ3::MGARDDecomposition<float, 3>(
                eb, [radius](double level_eb) { return SZ3::LinearQuantizer<float>(level_eb, radius); });
        },
        eb, mgardFields(conf.num));
}
