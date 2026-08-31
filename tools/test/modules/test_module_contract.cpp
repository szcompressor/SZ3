// Group-level acceptance checks. Every module in the library must pass the contract for its
// group; see docs/MODULE_ACCEPTANCE.md.

#include "ModuleContract.hpp"
#include "SZ3/encoder/ArithmeticEncoder.hpp"
#include "SZ3/encoder/BitshuffleEncoder.hpp"
#include "SZ3/encoder/BypassEncoder.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/encoder/HuffmanEncoderV2.hpp"
#include "SZ3/encoder/RunlengthEncoder.hpp"
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
