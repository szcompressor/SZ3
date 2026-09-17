// The two molecular-dynamics algorithms, through the public API a caller reaches them by.
//
// The integration suite runs these on the SDRBench fields, all of which are 1D or 2D, so the
// trajectory layout SZBioMDXtcDecomposition was written for -- dims {frames, atoms, xyz} -- is
// only ever reached through 3D climate fields. Those never have a frame whose values are all
// equal, so the fill-frame path that skips trailing frames and restores them on decompression
// does not run there at all. A trajectory buffer that holds fewer frames than it has room for
// does, which is the case a simulation writing through the HDF5 filter produces.

#include <cmath>
#include <cstdint>
#include <cstring>
#include <random>
#include <vector>

#include "SZ3/api/sz.hpp"
#include "gtest/gtest.h"

namespace {

constexpr SZ3::ALGO kAlgos[] = {SZ3::ALGO_BIOMD, SZ3::ALGO_BIOMDXTC};

// Reconstruction rounds to float, so the bound holds to within an ulp of the reconstructed value
// rather than exactly: 1.0014x at eb = 1e-4 here, and 1.0071x for ALGO_BIOMDXTC through the HDF5
// filter, where each chunk quantizes on its own.
constexpr double kBoundSlack = 1.01;

const char *algo_name(SZ3::ALGO algo) { return algo == SZ3::ALGO_BIOMD ? "ALGO_BIOMD" : "ALGO_BIOMDXTC"; }

/// A water-box trajectory: atoms placed on a jittered lattice, then walking a little each frame.
std::vector<float> make_trajectory(size_t frames, size_t atoms, uint32_t seed = 7) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> jitter(-0.02f, 0.02f);
    std::vector<float> positions(atoms * 3);
    const auto side = static_cast<size_t>(std::cbrt(static_cast<double>(atoms)) + 1);
    for (size_t a = 0; a < atoms; a++) {
        positions[a * 3 + 0] = static_cast<float>(a % side) * 0.31f;
        positions[a * 3 + 1] = static_cast<float>((a / side) % side) * 0.31f;
        positions[a * 3 + 2] = static_cast<float>(a / (side * side)) * 0.31f;
    }

    std::vector<float> traj(frames * atoms * 3);
    for (size_t f = 0; f < frames; f++) {
        for (size_t i = 0; i < atoms * 3; i++) {
            positions[i] += jitter(rng);
            traj[f * atoms * 3 + i] = positions[i];
        }
    }
    return traj;
}

/// Compress and decompress through the public API, returning the compressed bytes.
std::vector<char> round_trip(SZ3::ALGO algo, double eb, const std::vector<size_t> &dims,
                             const std::vector<float> &input, std::vector<float> &output) {
    SZ3::Config conf;
    conf.setDims(dims.begin(), dims.end());
    conf.cmprAlgo = algo;
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = eb;

    std::vector<char> compressed(SZ3::SZ_compress_size_bound<float>(conf));
    const size_t size = SZ_compress(conf, input.data(), compressed.data(), compressed.size());
    compressed.resize(size);

    output.assign(input.size(), 0.0f);
    float *out = output.data();
    SZ3::Config read_conf;
    SZ_decompress(read_conf, compressed.data(), compressed.size(), out);
    return compressed;
}

float max_abs_error(const std::vector<float> &a, const std::vector<float> &b) {
    float worst = 0;
    for (size_t i = 0; i < a.size(); i++) {
        worst = std::max(worst, std::fabs(a[i] - b[i]));
    }
    return worst;
}

}  // namespace

/// dims {frames, atoms, xyz} is the layout a trajectory arrives in.
TEST(SZ3_BioMD, TrajectoryRoundTripHonoursTheBound) {
    const std::vector<size_t> dims = {21, 2652, 3};
    const auto input = make_trajectory(dims[0], dims[1]);

    for (SZ3::ALGO algo : kAlgos) {
        for (double eb : {1e-1, 1e-2, 1e-3, 1e-4}) {
            std::vector<float> output;
            const auto compressed = round_trip(algo, eb, dims, input, output);
            ASSERT_EQ(output.size(), input.size()) << algo_name(algo) << " eb=" << eb;
            EXPECT_LE(max_abs_error(input, output), eb * kBoundSlack)
                << algo_name(algo) << " eb=" << eb;
            EXPECT_LT(compressed.size(), input.size() * sizeof(float))
                << algo_name(algo) << " eb=" << eb << " did not compress";
        }
    }
}

/// Trailing frames that hold one value are skipped on compression and written back on decompression.
TEST(SZ3_BioMD, TrailingFilledFramesAreRestored) {
    const std::vector<size_t> dims = {16, 500, 3};
    const size_t frame = dims[1] * dims[2];
    const size_t written_frames = 9;
    const float fill = 1.25f;

    auto input = make_trajectory(dims[0], dims[1]);
    for (size_t f = written_frames; f < dims[0]; f++) {
        std::fill(input.begin() + f * frame, input.begin() + (f + 1) * frame, fill);
    }

    for (SZ3::ALGO algo : kAlgos) {
        std::vector<float> output;
        const auto compressed = round_trip(algo, 1e-3, dims, input, output);
        ASSERT_EQ(output.size(), input.size()) << algo_name(algo);
        EXPECT_LE(max_abs_error(input, output), 1e-3 * kBoundSlack) << algo_name(algo);
        // The filled tail is written back verbatim, not quantized.
        for (size_t i = written_frames * frame; i < input.size(); i++) {
            ASSERT_EQ(output[i], fill) << algo_name(algo) << " at element " << i;
        }
        EXPECT_LT(compressed.size(), input.size() * sizeof(float)) << algo_name(algo);
    }
}

/// Everything after the first frame is fill, which is the largest the skipped range can be.
TEST(SZ3_BioMD, EveryFrameAfterTheFirstIsFill) {
    const std::vector<size_t> dims = {8, 300, 3};
    const size_t frame = dims[1] * dims[2];
    const float fill = -3.5f;

    auto input = make_trajectory(dims[0], dims[1]);
    std::fill(input.begin() + frame, input.end(), fill);

    for (SZ3::ALGO algo : kAlgos) {
        std::vector<float> output;
        round_trip(algo, 1e-3, dims, input, output);
        ASSERT_EQ(output.size(), input.size()) << algo_name(algo);
        EXPECT_LE(max_abs_error(input, output), 1e-3 * kBoundSlack) << algo_name(algo);
        for (size_t i = frame; i < input.size(); i++) {
            ASSERT_EQ(output[i], fill) << algo_name(algo) << " at element " << i;
        }
    }
}

/// A single frame is the shape a filter writes when it chunks a trajectory one frame at a time.
TEST(SZ3_BioMD, SingleFrameTrajectory) {
    const std::vector<size_t> dims = {1, 1024, 3};
    const auto input = make_trajectory(dims[0], dims[1]);

    for (SZ3::ALGO algo : kAlgos) {
        std::vector<float> output;
        round_trip(algo, 1e-3, dims, input, output);
        ASSERT_EQ(output.size(), input.size()) << algo_name(algo);
        EXPECT_LE(max_abs_error(input, output), 1e-3 * kBoundSlack) << algo_name(algo);
    }
}

/// The single-frame paths, which 1D and 2D data take.
TEST(SZ3_BioMD, OneAndTwoDimensionalInput) {
    const auto atoms = make_trajectory(1, 4096);

    for (SZ3::ALGO algo : kAlgos) {
        for (const std::vector<size_t> &dims :
             {std::vector<size_t>{atoms.size()}, std::vector<size_t>{4096, 3}}) {
            std::vector<float> output;
            round_trip(algo, 1e-3, dims, atoms, output);
            ASSERT_EQ(output.size(), atoms.size()) << algo_name(algo) << " " << dims.size() << "D";
            EXPECT_LE(max_abs_error(atoms, output), 1e-3 * kBoundSlack)
                << algo_name(algo) << " " << dims.size() << "D";
        }
    }
}

/// XtcBasedEncoder walks its magicInts table to the end whenever no entry fits, which every
/// input shorter than two atoms does.
TEST(SZ3_BioMD, InputsTooShortForOneTriplet) {
    for (SZ3::ALGO algo : kAlgos) {
        for (size_t n : {size_t{1}, size_t{2}, size_t{3}, size_t{5}, size_t{6}}) {
            std::vector<float> input(n);
            for (size_t i = 0; i < n; i++) {
                input[i] = 0.5f * static_cast<float>(i) - 1.0f;
            }
            std::vector<float> output;
            round_trip(algo, 1e-3, {n}, input, output);
            ASSERT_EQ(output.size(), n) << algo_name(algo) << " n=" << n;
            EXPECT_LE(max_abs_error(input, output), 1e-3 * kBoundSlack) << algo_name(algo) << " n=" << n;
        }
    }
}

/// The same input has to produce the same file, or a trajectory cannot be checksummed.
TEST(SZ3_BioMD, CompressionIsDeterministic) {
    const std::vector<size_t> dims = {5, 777, 3};
    const auto input = make_trajectory(dims[0], dims[1]);

    for (SZ3::ALGO algo : kAlgos) {
        std::vector<float> discard;
        const auto first = round_trip(algo, 1e-3, dims, input, discard);
        const auto second = round_trip(algo, 1e-3, dims, input, discard);
        ASSERT_EQ(first.size(), second.size()) << algo_name(algo);
        EXPECT_EQ(0, std::memcmp(first.data(), second.data(), first.size())) << algo_name(algo);
    }
}
