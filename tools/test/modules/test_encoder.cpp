#include <cmath>
#include <cstdint>

#include "SZ3/def.hpp"
#include "SZ3/utils/Config.hpp"  // must precede XtcBasedEncoder; XtcBased uses Config without including it
#include "SZ3/encoder/ArithmeticEncoder.hpp"
#include "SZ3/encoder/BitplaneEncoder.hpp"
#include "SZ3/encoder/BitplaneRLEEncoder.hpp"
#include "SZ3/encoder/BypassEncoder.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/encoder/HuffmanEncoderV2.hpp"
#include "SZ3/encoder/RunlengthEncoder.hpp"
#include "SZ3/encoder/BitshuffleEncoder.hpp"
#include "SZ3/encoder/SPERREncoder.hpp"
#include "SZ3/encoder/XtcBasedEncoder.hpp"
// NOTE: SZ3/encoder/ZFPEncoder.hpp is currently UNCOMPILABLE in isolation:
// it references `ZFP::IntCodec04<...>` but does not include
// `intcodec04.h`. Including it here breaks this whole TU. Skipped until
// the upstream header is fixed; documented in coverage matrix.
#include "gtest/gtest.h"

template <typename Encoder, typename T>
void runFunctionalTest() {
    int N = 1000;
    std::vector<SZ3::uchar> buffer_data(N * sizeof(T) * 2);
    std::vector<SZ3::uchar> buffer_conf(N * sizeof(T) * 2);
    std::vector<T> data(N);
    for (int i = 0; i < N; i++) {
        data[i] = i % 100;
    }

    size_t data_len = 0, conf_len = 0;
    {
        SZ3::uchar *buffer_data_pos = buffer_data.data();
        SZ3::uchar *buffer_conf_pos = buffer_conf.data();
        Encoder coder;
        coder.preprocess_encode(data, 100);
        coder.encode(data, buffer_data_pos);
        coder.save(buffer_conf_pos);
        data_len = buffer_data_pos - buffer_data.data();
        conf_len = buffer_conf_pos - buffer_conf.data();
    }
    {
        const SZ3::uchar *buffer_data_pos = buffer_data.data();
        const SZ3::uchar *buffer_conf_pos = buffer_conf.data();
        Encoder coder;
        coder.load(buffer_conf_pos, conf_len);
        auto dataDecoded = coder.decode(buffer_data_pos, N);
        for (int i = 0; i < N; i++) {
            EXPECT_EQ(data[i], dataDecoded[i]);
        }
    }
}

template <typename Encoder, typename T>
void runAllTest() {
    runFunctionalTest<Encoder, T>();
}

TEST(SZ3_EncoderTest, HuffmanEncoder) { runAllTest<SZ3::HuffmanEncoder<int>, int>(); }

TEST(SZ3_EncoderTest, RunlengthEncoder) { runAllTest<SZ3::RunlengthEncoder<int>, int>(); }

TEST(SZ3_EncoderTest, ArithmeticEncoder) { runAllTest<SZ3::ArithmeticEncoder<int>, int>(); }

TEST(SZ3_EncoderTest, BypassEncoder) { runAllTest<SZ3::BypassEncoder<int>, int>(); }


template<class T, SZ3::uint BITS, class BaseEncoder>
class BitshuffleEncoderWrapper : public BaseEncoder {
public:
    BitshuffleEncoderWrapper() : BaseEncoder(BITS) {}
};

TEST(SZ3_EncoderTest, BitshuffleEncoder) {
    runAllTest<BitshuffleEncoderWrapper<int, 4, SZ3::BitshuffleEncoder<int>>, int>();
    runAllTest<BitshuffleEncoderWrapper<int, 8, SZ3::BitshuffleEncoder<int>>, int>();
    runAllTest<BitshuffleEncoderWrapper<float, 4, SZ3::BitshuffleEncoder<float>>, float>();
    runAllTest<BitshuffleEncoderWrapper<float, 8, SZ3::BitshuffleEncoder<float>>, float>();
}

TEST(SZ3_EncoderTest, BitplaneEncoder) {
    runAllTest<SZ3::BitplaneEncoder<int>, int>();
}

TEST(SZ3_EncoderTest, BitplaneRLEEncoder) {
    runAllTest<SZ3::BitplaneRLEEncoder<int>, int>();
}

// Edge cases the generic harness doesn't cover.
template <typename Encoder>
static void roundtripBins(const std::vector<int>& bins) {
    std::vector<SZ3::uchar> data_buf(bins.size() * sizeof(int) * 4 + 64);
    std::vector<SZ3::uchar> conf_buf(64);
    SZ3::uchar* dp = data_buf.data();
    SZ3::uchar* cp = conf_buf.data();
    Encoder e;
    e.preprocess_encode(bins, 0);
    e.encode(bins, dp);
    e.save(cp);
    size_t conf_len = cp - conf_buf.data();
    const SZ3::uchar* cp_in = conf_buf.data();
    const SZ3::uchar* dp_in = data_buf.data();
    Encoder e2;
    e2.load(cp_in, conf_len);
    auto out = e2.decode(dp_in, bins.size());
    ASSERT_EQ(out.size(), bins.size());
    for (size_t i = 0; i < bins.size(); i++) EXPECT_EQ(out[i], bins[i]) << "i=" << i;
}

TEST(SZ3_EncoderTest, BitplaneEncoderEdgeCases) {
    roundtripBins<SZ3::BitplaneEncoder<int>>({});                       // empty
    roundtripBins<SZ3::BitplaneEncoder<int>>({0, 0, 0, 0, 0});          // all-zero
    roundtripBins<SZ3::BitplaneEncoder<int>>({-3, 0, 7, -1, 5, 2, 0});  // mixed signs
    roundtripBins<SZ3::BitplaneEncoder<int>>({65535, 0, 1});            // wide range
}

TEST(SZ3_EncoderTest, BitplaneRLEEncoderEdgeCases) {
    roundtripBins<SZ3::BitplaneRLEEncoder<int>>({});                       // empty
    roundtripBins<SZ3::BitplaneRLEEncoder<int>>({0, 0, 0, 0, 0});          // all-zero (sparse top planes)
    roundtripBins<SZ3::BitplaneRLEEncoder<int>>({-3, 0, 7, -1, 5, 2, 0});  // mixed signs
    // Force a mode mix by giving sparse top planes (favors RLE) + dense low planes (favors RAW).
    std::vector<int> mixed(2048);
    for (int i = 0; i < 2048; i++) mixed[i] = (i < 16) ? 65535 : (i & 1);
    roundtripBins<SZ3::BitplaneRLEEncoder<int>>(mixed);
}

// HuffmanEncoderV2 mirrors HuffmanEncoder's interface (preprocess_encode takes
// stateNum) but its preprocess_encode is non-virtual; default-construct + run.
TEST(SZ3_EncoderTest, HuffmanEncoderV2) { runAllTest<SZ3::HuffmanEncoderV2<int>, int>(); }

// SPERREncoder requires non-empty bins. 1D path uses SPECK1D_INT.
// NOTE: SPERREncoder<int, 2> and <int, 3> currently fail to link --
// `SPECK2D_INT::m_initialize_lists`, `m_clean_LIS`, `m_sorting_pass`, and the
// process_{S,P,I} member functions are not defined in the vendored
// sperr_core_impl.hpp (only SPECK1D_INT and SPECK3D_INT have impls in there).
// Documented in coverage matrix; 1D is the only roundtrip we can run today.
TEST(SZ3_EncoderTest, SPERREncoder1D) {
    std::vector<int> bins(512);
    for (size_t i = 0; i < bins.size(); i++) bins[i] = static_cast<int>(i % 17) - 8;
    roundtripBins<SZ3::SPERREncoder<int, 1>>(bins);
}

// XtcBasedEncoder is GROMACS-derived and operates on triplets (XYZ) of
// quantized integers. Run a tiny 4-triplet round trip.
TEST(SZ3_EncoderTest, XtcBasedEncoder) {
    roundtripBins<SZ3::XtcBasedEncoder<int>>({
        0, 0, 0,
        1, 0, 0,
        1, 1, 0,
        2, 1, 1,
    });
}
