// The two halves of the save/load contract, for the decompositions the MD algorithms use.
//
// save() has to write the same bytes for the same input, or a compressed file cannot be
// checksummed and two runs of the same pipeline produce different output.
//
// A member that save() writes but some compress() path never assigns takes whatever the memory
// held before, so the check has to control what that was: the decomposition is constructed over
// a buffer filled with two different patterns, and the two headers compared. Compressing twice
// on the heap instead would only catch it when the allocator happens to hand back dirty memory.

#include <cmath>
#include <cstring>
#include <memory>
#include <new>
#include <random>
#include <vector>

#include "SZ3/decomposition/BlockwiseDecomposition.hpp"
#include "SZ3/decomposition/SZBioMDDecomposition.hpp"
#include "SZ3/predictor/ComposedPredictor.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/predictor/RegressionPredictor.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "gtest/gtest.h"

namespace {

std::vector<float> noise(size_t count, uint32_t seed = 11) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> spread(-5.0f, 5.0f);
    std::vector<float> data(count);
    for (auto &value : data) {
        value = spread(rng);
    }
    return data;
}

/// A ramp is what a regression predictor fits best, so it is also what makes it store the most.
std::vector<float> ramp(size_t count, uint32_t seed = 3) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> jitter(-1.0f, 1.0f);
    std::vector<float> data(count);
    for (size_t i = 0; i < count; i++) {
        data[i] = static_cast<float>(i) * 0.5f + jitter(rng);
    }
    return data;
}

/// The compressor allocates from size_est() and save() writes without checking, so check the
/// bound directly: save with room to spare, so an overrun is this assertion, not a crash.
template <class Decomposition>
void expect_save_stays_within_size_est(Decomposition &decomposition, const SZ3::Config &conf, std::vector<float> data,
                                       const char *what) {
    SZ3::Config work = conf;
    decomposition.compress(work, data.data());

    const size_t declared = decomposition.size_est();
    std::vector<SZ3::uchar> buffer(declared + (1u << 20));
    SZ3::uchar *cursor = buffer.data();
    decomposition.save(cursor);

    const size_t written = static_cast<size_t>(cursor - buffer.data());
    EXPECT_LE(written, declared) << what << ": save() wrote " << written << " bytes where size_est() promised "
                                 << declared;
}

/// Compress over storage holding `poison`, and return what save() writes.
template <class Decomposition, class MakeQuantizer>
std::vector<SZ3::uchar> header_after_compress(const std::vector<size_t> &dims, unsigned char poison,
                                              MakeQuantizer make_quantizer) {
    SZ3::Config conf;
    conf.setDims(dims.begin(), dims.end());
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = 1e-3;

    std::vector<unsigned char> storage(sizeof(Decomposition) + alignof(Decomposition));
    std::memset(storage.data(), poison, storage.size());
    void *place = storage.data();
    size_t room = storage.size();
    place = std::align(alignof(Decomposition), sizeof(Decomposition), place, room);

    auto data = noise(conf.num);
    auto *decomposition = new (place) Decomposition(conf, make_quantizer(conf));
    decomposition->compress(conf, data.data());

    std::vector<SZ3::uchar> header(decomposition->size_est() + 4096);
    SZ3::uchar *cursor = header.data();
    decomposition->save(cursor);
    header.resize(static_cast<size_t>(cursor - header.data()));
    decomposition->~Decomposition();
    return header;
}

template <class Decomposition, class MakeQuantizer>
void expect_header_depends_only_on_input(const std::vector<size_t> &dims, MakeQuantizer make_quantizer) {
    const auto over_zeros = header_after_compress<Decomposition>(dims, 0x00, make_quantizer);
    const auto over_ones = header_after_compress<Decomposition>(dims, 0xcd, make_quantizer);
    ASSERT_EQ(over_zeros.size(), over_ones.size());
    EXPECT_EQ(over_zeros, over_ones) << "save() wrote bytes that came from the memory it was built over";
}

/// load() must charge remaining_length for exactly what it read, which is what keeps every
/// later parse inside the buffer.
template <class Decomposition, class MakeQuantizer>
void expect_load_charges_what_it_reads(const std::vector<size_t> &dims, MakeQuantizer make_quantizer) {
    SZ3::Config conf;
    conf.setDims(dims.begin(), dims.end());
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = 1e-3;

    Decomposition writer(conf, make_quantizer(conf));
    auto data = noise(conf.num);
    writer.compress(conf, data.data());

    std::vector<SZ3::uchar> stream(writer.size_est() + 4096);
    SZ3::uchar *write_cursor = stream.data();
    writer.save(write_cursor);
    const size_t written = static_cast<size_t>(write_cursor - stream.data());

    Decomposition reader(conf, make_quantizer(conf));
    const SZ3::uchar *read_cursor = stream.data();
    size_t remaining = stream.size();
    reader.load(read_cursor, remaining);

    const size_t advanced = static_cast<size_t>(read_cursor - stream.data());
    EXPECT_EQ(advanced, written) << "load() did not consume what save() wrote";
    EXPECT_EQ(stream.size() - remaining, advanced)
        << "load() advanced " << advanced << " bytes but charged " << stream.size() - remaining;
}

// SZBioMDDecomposition takes no quantizer: this one takes the helpers' second argument and drops it
template <SZ3::uint N>
struct BioMD : SZ3::SZBioMDDecomposition<float, N> {
    BioMD(const SZ3::Config &conf, int) : SZ3::SZBioMDDecomposition<float, N>(conf) {}
};
auto no_quantizer = [](const SZ3::Config &) { return 0; };

}  // namespace

TEST(SZ3_DecompositionSaveLoad, BioMDTwoDimensions) {
    expect_header_depends_only_on_input<BioMD<2>>({777, 3}, no_quantizer);
}

TEST(SZ3_DecompositionSaveLoad, BioMDThreeDimensions) {
    expect_header_depends_only_on_input<BioMD<3>>({5, 777, 3}, no_quantizer);
}

TEST(SZ3_DecompositionSaveLoad, BioMDChargesWhatItReads) {
    expect_load_charges_what_it_reads<BioMD<2>>({777, 3}, no_quantizer);
    expect_load_charges_what_it_reads<BioMD<3>>({5, 777, 3}, no_quantizer);
}

// The mode bytes and the box flag hold 0 or 1; any other value is a corrupt stream, not a value to truncate.
TEST(SZ3_DecompositionSaveLoad, BioMDRefusesUnknownModeAndBoxFlag) {
    const std::vector<size_t> dims = {5, 777, 3};
    const auto header = header_after_compress<BioMD<3>>(dims, 0x00, no_quantizer);
    // coded frames (8), fill value (4), step (8), water distances (8 + 8), water count (8), then 4 modes and the flag
    const size_t modes = 44, box = modes + SZ3::biomd::NUM_GROUPS;
    SZ3::Config conf;
    conf.setDims(dims.begin(), dims.end());
    for (size_t at : {modes, modes + 3, box}) {
        auto corrupt = header;
        ASSERT_LE(corrupt[at], 1) << "byte " << at << " is not where the mode bytes and the box flag are";
        corrupt[at] = 2;
        BioMD<3> reader(conf, 0);
        const SZ3::uchar *cursor = corrupt.data();
        size_t remaining = corrupt.size();
        EXPECT_THROW(reader.load(cursor, remaining), std::runtime_error) << "byte " << at;
    }
}

// BlockwiseDecomposition declared nothing while save() serialised the regression coefficients,
// which grow as n/blockSize. Below the default block size that ran off the compressor's buffer.
TEST(SZ3_DecompositionSaveLoad, BlockwiseSizeEstBoundsSaveAtSmallBlockSizes) {
    for (int blockSize : {2, 3, 4, 6, 8, 128}) {
        for (double eb : {1e-1, 1e-2, 1e-3, 1e-4}) {
            SZ3::Config conf(16384);
            conf.errorBoundMode = SZ3::EB_ABS;
            conf.absErrorBound = eb;
            conf.blockSize = blockSize;

            std::vector<std::shared_ptr<SZ3::concepts::PredictorInterface<float, 1>>> predictors;
            predictors.push_back(std::make_shared<SZ3::LorenzoPredictor<float, 1, 1>>(eb));
            predictors.push_back(
                std::make_shared<SZ3::RegressionPredictor<float, 1>>(static_cast<SZ3::uint>(blockSize), eb));

            auto decomposition =
                SZ3::make_decomposition_blockwise<float, 1>(conf, SZ3::ComposedPredictor<float, 1>(predictors),
                                                            SZ3::LinearQuantizer<float>(eb, conf.quantbinCnt / 2));

            SCOPED_TRACE("BlockSize=" + std::to_string(blockSize) + " eb=" + std::to_string(eb));
            expect_save_stays_within_size_est(decomposition, conf, ramp(conf.num), "BlockwiseDecomposition");
        }
    }
}

// A block RegressionPredictor declines falls back to Lorenzo, which reads a neighbour outside the block: with only
// the regression predictor's padding it read past the data in 3D and 4D.
template <SZ3::uint N>
void expect_regression_only_round_trip(const std::vector<size_t> &dims) {
    SZ3::Config conf;
    conf.setDims(dims.begin(), dims.end());
    const double eb = 1e-3;
    auto make = [&] {
        return SZ3::make_decomposition_blockwise<float, N>(conf, SZ3::RegressionPredictor<float, N>(conf.blockSize, eb),
                                                           SZ3::LinearQuantizer<float>(eb, conf.quantbinCnt / 2));
    };
    const auto input = noise(conf.num);
    auto data = input;
    auto writer = make();
    auto bins = writer.compress(conf, data.data());
    std::vector<SZ3::uchar> stream(writer.size_est());
    SZ3::uchar *write_cursor = stream.data();
    writer.save(write_cursor);
    auto reader = make();
    const SZ3::uchar *read_cursor = stream.data();
    size_t remaining = write_cursor - stream.data();
    reader.load(read_cursor, remaining);
    std::vector<float> output(conf.num);
    reader.decompress(conf, bins, output.data());
    for (size_t i = 0; i < conf.num; i++) {
        ASSERT_LE(std::fabs(output[i] - input[i]), eb) << dims.size() << "D, value " << i;
    }
}

TEST(SZ3_DecompositionSaveLoad, BlockwiseRegressionOnlyFallback) {
    expect_regression_only_round_trip<3>({40, 41, 43});
    expect_regression_only_round_trip<4>({12, 13, 14, 15});
}

TEST(SZ3_DecompositionSaveLoad, BioMDSizeEstBoundsSave) {
    for (double eb : {1e-1, 1e-2, 1e-3, 1e-4}) {
        SZ3::Config conf(5, 777, 3);
        conf.errorBoundMode = SZ3::EB_ABS;
        conf.absErrorBound = eb;
        SCOPED_TRACE("eb=" + std::to_string(eb));

        BioMD<3> biomd(conf, 0);
        auto coordinates = ramp(conf.num);  // up to 12 nm, within BIOMD's lattice (64 nm at 1e-4)
        for (auto &v : coordinates) v *= 0.002f;
        expect_save_stays_within_size_est(biomd, conf, coordinates, "SZBioMDDecomposition");
    }
}
