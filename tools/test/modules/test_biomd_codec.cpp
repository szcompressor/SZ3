// ALGO_BIOMD on synthetic molecular-dynamics systems: rigid 3- and 4-site water, bonded chains and ions, as GROMACS
// lays them out (molecules in order, water in one block): the bound, the input it refuses, recompression, and corrupt
// or truncated streams, which are refused without reading or writing out of bounds.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

#include "SZ3/api/sz.hpp"
#include "gtest/gtest.h"

namespace {

struct SystemSpec {
    size_t waters = 800;
    bool four_site = false;
    size_t chains = 2, chain_len = 60;  // heavy atoms per chain, each with one hydrogen
    size_t ions = 10;
    size_t frames = 1;
    uint32_t seed = 1;
};

struct Vec {
    double x, y, z;
};
Vec operator+(Vec a, Vec b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
Vec operator*(double s, Vec a) { return {s * a.x, s * a.y, s * a.z}; }
double dot(Vec a, Vec b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
Vec unit(Vec a) { return (1.0 / std::sqrt(dot(a, a))) * a; }
Vec cross(Vec a, Vec b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }

// frames x atoms x 3 floats, in nm
std::vector<float> make_system(const SystemSpec &s, size_t *atoms) {
    std::mt19937 rng(s.seed);
    std::uniform_real_distribution<double> U(-1, 1);
    std::normal_distribution<double> G(0, 1);
    auto rvec = [&]() { return unit(Vec{G(rng), G(rng), G(rng)}); };
    const size_t per_water = s.four_site ? 4 : 3;
    const size_t n = s.chains * s.chain_len * 2 + s.waters * per_water + s.ions;
    *atoms = n;
    const double box = std::cbrt(double(n) * 0.01);
    const double r = 0.09572, theta = 104.52 * 3.14159265358979323846 / 180, a_vs = 0.128;  // TIP4P geometry
    const double bonds[3] = {0.1529, 0.1335, 0.1471}, bh = 0.109;
    std::vector<Vec> O(s.waters), u(s.waters), v(s.waters), chain_dir(s.chains * s.chain_len), ion(s.ions);
    std::vector<Vec> chain0(s.chains);
    for (size_t w = 0; w < s.waters; w++) {
        O[w] = {box * (U(rng) + 1) / 2, box * (U(rng) + 1) / 2, box * (U(rng) + 1) / 2};
        u[w] = rvec();
        v[w] = unit(cross(u[w], rvec()));
    }
    for (auto &d : chain_dir) d = rvec();
    for (auto &c : chain0) c = {box * (U(rng) + 1) / 2, box * (U(rng) + 1) / 2, box * (U(rng) + 1) / 2};
    for (auto &i : ion) i = {box * (U(rng) + 1) / 2, box * (U(rng) + 1) / 2, box * (U(rng) + 1) / 2};
    std::vector<float> x(s.frames * n * 3);
    for (size_t f = 0; f < s.frames; f++) {
        float *X = &x[f * n * 3];
        size_t k = 0;
        auto put = [&](Vec p) {
            X[3 * k] = float(p.x);
            X[3 * k + 1] = float(p.y);
            X[3 * k + 2] = float(p.z);
            k++;
        };
        for (size_t c = 0; c < s.chains; c++) {
            Vec p = chain0[c];
            for (size_t a = 0; a < s.chain_len; a++) {
                const Vec d = chain_dir[c * s.chain_len + a];
                p = p + bonds[a % 3] * d;
                put(p);
                put(p + bh * unit(cross(d, rvec())));
            }
        }
        for (size_t w = 0; w < s.waters; w++) {
            const Vec h1 = O[w] + r * u[w], h2 = O[w] + r * (std::cos(theta) * u[w] + std::sin(theta) * v[w]);
            put(O[w]);
            put(h1);
            put(h2);
            if (s.four_site) put(O[w] + a_vs * ((h1 + (-1.0) * O[w]) + (h2 + (-1.0) * O[w])));
        }
        for (auto &p : ion) put(p);
        // the next frame: molecules move and turn a little
        for (size_t w = 0; w < s.waters; w++) {
            O[w] = O[w] + 0.003 * rvec();
            u[w] = unit(u[w] + 0.05 * rvec());
            v[w] = unit(cross(u[w], cross(v[w], u[w])));
        }
        for (auto &d : chain_dir) d = unit(d + 0.03 * rvec());
        for (auto &p : ion) p = p + 0.005 * rvec();
    }
    return x;
}

template <class T>
struct RoundTrip {
    std::vector<char> bytes;
    std::vector<T> out;
    SZ3::ALGO algo;
    double max_err;
};

template <class T>
RoundTrip<T> round_trip(const std::vector<T> &x, size_t frames, size_t atoms, double eb) {
    SZ3::Config conf;
    if (frames == 1)
        conf = SZ3::Config(atoms, 3);
    else
        conf = SZ3::Config(frames, atoms, 3);
    conf.cmprAlgo = SZ3::ALGO_BIOMD;
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = eb;
    size_t size = 0;
    char *c = SZ_compress(conf, x.data(), size);
    RoundTrip<T> r;
    r.bytes.assign(c, c + size);
    delete[] c;
    SZ3::Config dconf;
    r.out.resize(x.size());
    T *o = r.out.data();
    SZ_decompress(dconf, r.bytes.data(), r.bytes.size(), o);
    r.algo = static_cast<SZ3::ALGO>(dconf.cmprAlgo);  // the algorithm the stream holds
    r.max_err = 0;
    for (size_t i = 0; i < x.size(); i++) r.max_err = std::max(r.max_err, std::fabs(double(r.out[i]) - double(x[i])));
    return r;
}

TEST(BioMD, WithinBoundAndBeatsFourBytesPerValue) {
    for (bool four : {false, true})
        for (size_t frames : {size_t(1), size_t(5), size_t(20)})
            for (double eb : {5e-4, 1e-3, 1e-5, 0.1}) {
                SystemSpec s;
                s.four_site = four;
                s.frames = frames;
                size_t n;
                const auto x = make_system(s, &n);
                const auto r = round_trip(x, frames, n, eb);
                SCOPED_TRACE(testing::Message() << "four_site=" << four << " frames=" << frames << " eb=" << eb);
                EXPECT_EQ(r.algo, SZ3::ALGO_BIOMD);
                EXPECT_LE(r.max_err, eb);
                if (eb == 5e-4) {
                    EXPECT_GT(double(x.size() * 4) / double(r.bytes.size()), 3.0);
                }
            }
}

TEST(BioMD, LooseBoundsStayWithinTheBound) {
    SystemSpec s;
    s.frames = 20;
    size_t n;
    auto x = make_system(s, &n);
    for (size_t t = 1; t < s.frames; t++)  // frames that do not move
        std::copy(x.begin(), x.begin() + 3 * n, x.begin() + t * 3 * n);
    for (double eb : {0.02, 0.1}) {
        const auto r = round_trip(x, s.frames, n, eb);
        EXPECT_EQ(r.algo, SZ3::ALGO_BIOMD);
        EXPECT_LE(r.max_err, eb);
    }
}

TEST(BioMD, DoubleInput) {
    SystemSpec s;
    s.frames = 4;
    size_t n;
    const auto xf = make_system(s, &n);
    const std::vector<double> x(xf.begin(), xf.end());
    for (double eb : {5e-4, 1e-7}) {
        const auto r = round_trip(x, s.frames, n, eb);
        EXPECT_LE(r.max_err, eb) << "eb=" << eb;
    }
}

TEST(BioMD, EdgeCasesStayWithinBound) {
    std::mt19937 rng(3);
    std::uniform_real_distribution<float> U(0.f, 5.f);
    std::vector<std::pair<std::string, std::vector<float>>> cases;
    for (size_t n : {1, 2, 3, 4, 5, 7}) {
        std::vector<float> x(3 * n);
        for (auto &v : x) v = U(rng);
        cases.push_back({"tiny N=" + std::to_string(n), x});
    }
    cases.push_back({"constant", std::vector<float>(3000, 1.234f)});
    cases.push_back({"zeros", std::vector<float>(3000, 0.f)});
    SystemSpec s;
    size_t n;
    const auto base = make_system(s, &n);
    {
        auto x = base;  // waters split over a periodic boundary
        for (size_t i = s.chains * s.chain_len * 2; i + 2 < n; i += 30) x[3 * (i + 1)] += 3.1f;
        cases.push_back({"broken waters", x});
    }
    {
        auto x = base;  // atoms in no particular order
        std::vector<size_t> perm(n);
        std::iota(perm.begin(), perm.end(), 0);
        std::shuffle(perm.begin(), perm.end(), rng);
        for (size_t i = 0; i < n; i++)
            for (int c = 0; c < 3; c++) x[3 * i + c] = base[3 * perm[i] + c];
        cases.push_back({"shuffled", x});
    }
    {
        auto x = base;
        for (auto &v : x) v += 1e5f;
        cases.push_back({"offset 1e5 nm", x});
    }
    {
        auto x = base;
        for (auto &v : x) v *= 10.f;  // Angstrom
        cases.push_back({"Angstrom", x});
    }
    for (const auto &c : cases) {
        const size_t atoms = c.second.size() / 3;
        const auto r = round_trip(c.second, 1, atoms, 5e-4);
        EXPECT_EQ(r.algo, SZ3::ALGO_BIOMD) << c.first;
        EXPECT_LE(r.max_err, 5e-4) << c.first;
    }
}

// Input BIOMD does not code goes to LORENZO_REG, within the bound: shapes other than {frames, atoms, 3}, coordinates
// beyond the lattice the bound allows, and NaN or Inf outside trailing fill, which come back exactly. The input stays
// as it was.
template <class T>
std::pair<SZ3::ALGO, std::vector<T>> compress_decompress(const std::vector<T> &x, const std::vector<size_t> &dims,
                                                         double eb) {
    SZ3::Config conf;
    conf.setDims(dims.begin(), dims.end());
    conf.cmprAlgo = SZ3::ALGO_BIOMD;
    conf.errorBoundMode = SZ3::EB_ABS;
    conf.absErrorBound = eb;
    size_t size = 0;
    const std::vector<T> before = x;
    std::unique_ptr<char[]> cmp(SZ_compress(conf, x.data(), size));
    EXPECT_EQ(memcmp(before.data(), x.data(), x.size() * sizeof(T)), 0) << "the input changed";
    SZ3::Config dec;
    std::unique_ptr<T[]> y(SZ_decompress<T>(dec, cmp.get(), size));
    return {static_cast<SZ3::ALGO>(dec.cmprAlgo), std::vector<T>(y.get(), y.get() + x.size())};
}

template <class T>
void expect_within(const std::vector<T> &x, const std::vector<T> &y, double eb) {
    for (size_t i = 0; i < x.size(); i++) ASSERT_LE(std::fabs(double(y[i]) - double(x[i])), eb) << "at " << i;
}

TEST(BioMD, FallsBackOnWhatItDoesNotCode) {
    SystemSpec s;
    size_t n;
    const auto x = make_system(s, &n);
    EXPECT_EQ(compress_decompress(x, {n, 3}, 5e-4).first, SZ3::ALGO_BIOMD);
    for (const std::vector<size_t> &dims : {std::vector<size_t>{x.size()}, std::vector<size_t>{3, n},
                                            std::vector<size_t>{2, 5, n / 10, 3}}) {  // n = 2650
        const auto r = compress_decompress(x, dims, 5e-4);
        EXPECT_EQ(r.first, SZ3::ALGO_LORENZO_REG);
        expect_within(x, r.second, 5e-4);
    }
    {
        auto y = x;
        for (auto &v : y) v += 1e6f;  // beyond the lattice of 5e-4
        const auto r = compress_decompress(y, {n, 3}, 5e-4);
        EXPECT_EQ(r.first, SZ3::ALGO_LORENZO_REG);
        expect_within(y, r.second, 5e-4);
    }
    for (float bad : {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::infinity()}) {
        auto y = x;
        y[100] = bad;
        auto r = compress_decompress(y, {n, 3}, 5e-4);
        EXPECT_EQ(r.first, SZ3::ALGO_LORENZO_REG);
        EXPECT_EQ(memcmp(&r.second[100], &bad, sizeof(float)), 0);
        y[100] = r.second[100] = 0;
        expect_within(y, r.second, 5e-4);
    }
    // a bound so large that the lattice would pass the largest value of the type
    for (const auto &c : {std::make_pair(3e38f, 1e38), std::make_pair(1.0f, 1e30)}) {
        const std::vector<float> y = {c.first, 0, 0, 0, c.first, 0};
        const auto r = compress_decompress(y, {2, 3}, c.second);
        for (float v : r.second) EXPECT_TRUE(std::isfinite(v));
        expect_within(y, r.second, c.second);
    }
    const std::vector<double> z = {1.5e308, 0, 0, 0, 1.5e308, 0};
    const auto r = compress_decompress(z, {2, 3}, 5e307);
    for (double v : r.second) EXPECT_TRUE(std::isfinite(v));
    expect_within(z, r.second, 5e307);
}

// Config drops dimensions of 1: one atom arrives as {frames, 3}, one atom of one frame as {3}.
TEST(BioMD, OneAtom) {
    for (size_t frames : {size_t(1), size_t(5)}) {
        std::vector<float> x(frames * 3);
        for (size_t i = 0; i < x.size(); i++) x[i] = 1.0f + 0.01f * float(i);
        SZ3::Config conf(frames, 1, 3);
        conf.cmprAlgo = SZ3::ALGO_BIOMD;
        conf.errorBoundMode = SZ3::EB_ABS;
        conf.absErrorBound = 5e-4;
        size_t size = 0;
        char *c = SZ_compress(conf, x.data(), size);
        std::vector<float> y(x.size());
        float *o = y.data();
        SZ3::Config dconf;
        SZ_decompress(dconf, c, size, o);
        delete[] c;
        EXPECT_EQ(dconf.cmprAlgo, SZ3::ALGO_BIOMD);
        for (size_t i = 0; i < x.size(); i++) EXPECT_LE(std::fabs(y[i] - x[i]), 5e-4) << frames << " frames";
    }
}

TEST(BioMD, TrailingFillFramesComeBackExactly) {
    SystemSpec s;
    s.frames = 8;
    size_t n;
    auto x = make_system(s, &n);
    std::fill(x.begin() + 5 * n * 3, x.end(), 0.f);  // a chunk with 5 frames written
    const auto r = round_trip(x, s.frames, n, 5e-4);
    EXPECT_LE(r.max_err, 5e-4);
    for (size_t i = 5 * n * 3; i < x.size(); i++) ASSERT_EQ(r.out[i], 0.f);
}

// GROMACS rewrites a chunk after restarting from a checkpoint: decompress, append, compress again.
TEST(BioMD, RecompressionChangesNothing) {
    for (size_t frames : {size_t(1), size_t(6)})
        for (double eb : {5e-4, 1e-3}) {
            SystemSpec s;
            s.frames = frames;
            size_t n;
            const auto x = make_system(s, &n);
            const auto a = round_trip(x, frames, n, eb);
            const auto b = round_trip(a.out, frames, n, eb);
            EXPECT_EQ(0, memcmp(a.out.data(), b.out.data(), a.out.size() * sizeof(float))) << frames << " " << eb;
            EXPECT_LE(b.max_err, eb);
        }
}

// The lattice depends on the bound alone: a decoded value that crosses a power of two, or a frame appended to a chunk
// that raises its largest |x|, lands on the same lattice points when compressed again.
TEST(BioMD, RecompressionAcrossAPowerOfTwoStaysWithinTheBound) {
    const double eb = 5e-4;
    SystemSpec s;
    s.frames = 4;
    size_t n;
    auto x = make_system(s, &n);
    size_t k = 0;
    for (size_t i = 0; i < x.size(); i++)
        if (std::fabs(x[i]) > std::fabs(x[k])) k = i;
    const double p2 = std::exp2(std::ceil(std::log2(std::fabs(double(x[k])))));
    x[k] = float(std::copysign(p2 - 1e-4, double(x[k])));
    auto cur = x;
    for (int round = 0; round < 3; round++) {
        cur = round_trip(cur, s.frames, n, eb).out;
        for (size_t i = 0; i < x.size(); i++) ASSERT_LE(std::fabs(double(cur[i]) - double(x[i])), eb) << round;
    }
    // two frames and two of fill, then a third frame whose largest |x| is past the power of two
    std::vector<float> chunk(x.begin(), x.begin() + 2 * n * 3);
    chunk.resize(4 * n * 3, 0.f);
    auto dec = round_trip(chunk, 4, n, eb).out;
    for (size_t i = 0; i < n * 3; i++) dec[2 * n * 3 + i] = x[2 * n * 3 + i];
    dec[2 * n * 3] = float(p2 + 0.5);
    const auto again = round_trip(dec, 4, n, eb).out;
    for (size_t i = 0; i < 2 * n * 3; i++) ASSERT_LE(std::fabs(double(again[i]) - double(x[i])), eb) << i;
}

TEST(BioMD, CorruptStreamsAreRefusedOrDecodeInBounds) {
    SystemSpec s;
    s.frames = 3;
    s.four_site = true;
    size_t n;
    const auto x = make_system(s, &n);
    const auto r = round_trip(x, s.frames, n, 5e-4);
    // SZ3's stream: 16-byte header (magic, version, payload size), the payload, then the config
    uint64_t payload = 0;
    memcpy(&payload, r.bytes.data() + 8, 8);
    ASSERT_LT(16 + payload, r.bytes.size());
    std::vector<float> out(x.size());
    std::mt19937 rng(9);
    auto attempt = [&](std::vector<char> b) {
        try {
            SZ3::Config dconf;
            float *o = out.data();
            SZ_decompress(dconf, b.data(), b.size(), o);
        } catch (std::exception &) {
        }
    };
    for (size_t cut = 16; cut < 16 + payload; cut += std::max<size_t>(1, payload / 300)) {  // payload truncated
        std::vector<char> b(r.bytes.begin(), r.bytes.begin() + cut);
        b.insert(b.end(), r.bytes.begin() + 16 + payload, r.bytes.end());
        const uint64_t p = cut - 16;
        memcpy(b.data() + 8, &p, 8);
        attempt(b);
    }
    for (int k = 0; k < 300; k++) {  // payload bytes changed
        auto b = r.bytes;
        for (int m = 0; m < 4; m++) b[16 + rng() % payload] ^= char(1 + rng() % 255);
        attempt(b);
    }
    SUCCEED();
}

// Found by an audit: the rounding of x / step to the lattice and of q step back to double took more than the one ulp
// of slack the step left.
TEST(BioMD, DoubleNearHalfLatticePointsStaysWithinBound) {
    const std::vector<double> x = {510.39700369221237, 0, 0, 467.45074901980968, 0, 0};
    const double eb = 5.8345702187110699e-06;
    const auto r = round_trip(x, 1, 2, eb);
    EXPECT_EQ(r.algo, SZ3::ALGO_BIOMD);
    EXPECT_LE(r.max_err, eb);
}

// A water whose O and H2 swap ends of the coordinate range between frames (2^28 lattice steps apart): the H2 symbol
// against the previous frame must still fit 32 bits.
template <class T>
void water_jumps_across_the_range() {
    const double eb = std::ldexp(1.0, -10), M = std::ldexp(1.0, 18);  // M / step = 2^28
    const size_t waters = 40, atoms = 3 * waters + 2, frames = 3;
    std::vector<T> x(frames * atoms * 3);
    std::mt19937 rng(1);
    std::uniform_real_distribution<double> U(0, 3);
    const double r = 0.09572, th = 104.52 * std::acos(-1.0) / 180;
    std::vector<double> O(3 * waters);
    for (auto &v : O) v = U(rng);
    for (size_t f = 0; f < frames; f++) {
        T *X = &x[f * atoms * 3];
        for (size_t w = 0; w < waters; w++) {
            double o[3] = {O[3 * w] + 0.001 * f, O[3 * w + 1], O[3 * w + 2]};
            double h1[3] = {o[0], o[1] + r, o[2]}, h2[3] = {o[0], o[1] + r * std::cos(th), o[2] + r * std::sin(th)};
            if (w == 7 && f > 0) o[0] = h1[0] = (f == 1 ? M : -M), h2[0] = (f == 1 ? -M : M);
            for (int c = 0; c < 3; c++)
                X[9 * w + c] = T(o[c]), X[9 * w + 3 + c] = T(h1[c]), X[9 * w + 6 + c] = T(h2[c]);
        }
        X[9 * waters] = T(M), X[9 * waters + 3] = T(-M);  // two ions fix max|x|
    }
    const auto rt = round_trip(x, frames, atoms, eb);
    EXPECT_EQ(rt.algo, SZ3::ALGO_BIOMD);
    EXPECT_LE(rt.max_err, eb);
}
TEST(BioMD, WaterJumpsAcrossTheCoordinateRange) {
    water_jumps_across_the_range<float>();
    water_jumps_across_the_range<double>();
}

// The O-H distance at 16384 - 1e-5 lattice steps: its square rounds to 2^28, which load() refuses, so the codec must
// not code the water on it.
TEST(BioMD, WaterGeometryAtTheLatticeLimitDecodes) {
    SystemSpec s;
    size_t n;
    const auto x = make_system(s, &n);
    const auto r = round_trip(x, s.frames, n, 3.3979795989402915e-06);
    EXPECT_LE(r.max_err, 3.3979795989402915e-06);
}

// An unwritten chunk is all fill, NaN included.
TEST(BioMD, ChunkOfFillOnly) {
    const size_t frames = 3, atoms = 10;
    const std::vector<float> x(frames * atoms * 3, std::numeric_limits<float>::quiet_NaN());
    const auto r = round_trip(x, frames, atoms, 5e-4);
    EXPECT_EQ(r.algo, SZ3::ALGO_BIOMD);
    for (float v : r.out) EXPECT_TRUE(std::isnan(v));
}

// Few atoms and many frames: what each frame stores beyond its symbols must not outgrow the data.
TEST(BioMD, FewAtomsManyFrames) {
    const size_t frames = 5000, atoms = 2;
    std::vector<float> x(frames * atoms * 3);
    std::mt19937 rng(3);
    std::normal_distribution<float> step(0, 0.01f);
    for (size_t i = 0; i < x.size(); i++) x[i] = (i < atoms * 3 ? 1.0f : x[i - atoms * 3]) + step(rng);
    const auto r = round_trip(x, frames, atoms, 5e-4);
    EXPECT_EQ(r.algo, SZ3::ALGO_BIOMD);
    EXPECT_LE(r.max_err, 5e-4);
}

}  // namespace
