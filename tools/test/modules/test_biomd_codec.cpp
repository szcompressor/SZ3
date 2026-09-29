// ALGO_BIOMD on synthetic molecular-dynamics systems: rigid 3- and 4-site water, bonded chains and ions, as GROMACS
// lays them out (molecules in order, water in one block). Through the public API: the bound, the input it refuses,
// and recompression. Through the codec: every SIMD level gives the same bytes, and corrupt or
// truncated streams are refused without reading or writing out of bounds.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

#include "SZ3/api/sz.hpp"
#include "SZ3/compressor/specialized/biomd/BioMDCodec.hpp"
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
                if (eb == 5e-4) EXPECT_GT(double(x.size() * 4) / double(r.bytes.size()), 3.0);
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
        if (atoms > 1)
            EXPECT_EQ(r.algo, SZ3::ALGO_BIOMD) << c.first;  // {1, 3} collapses to 1D, which goes to ALGO_INTERP_LORENZO
        EXPECT_LE(r.max_err, 5e-4) << c.first;
    }
}

TEST(BioMD, RefusesWhatItCannotCode) {
    auto compress = [](const std::vector<float> &x, const std::vector<size_t> &dims) {
        SZ3::Config conf;
        conf.setDims(dims.begin(), dims.end());
        conf.cmprAlgo = SZ3::ALGO_BIOMD;
        conf.errorBoundMode = SZ3::EB_ABS;
        conf.absErrorBound = 5e-4;
        size_t size = 0;
        delete[] SZ_compress(conf, x.data(), size);
    };
    SystemSpec s;
    size_t n;
    const auto x = make_system(s, &n);
    EXPECT_THROW(compress(x, {x.size()}), std::invalid_argument);
    EXPECT_NO_THROW(compress(x, {n, 3}));
    EXPECT_THROW(compress(x, {3, n}), std::invalid_argument);
    {
        auto y = x;
        for (auto &v : y) v += 1e6f;  // 2^28 lattice steps of 2 eb are not enough
        EXPECT_THROW(compress(y, {n, 3}), std::runtime_error);
    }
    for (float bad : {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::infinity()}) {
        auto y = x;
        y[100] = bad;
        EXPECT_THROW(compress(y, {n, 3}), std::runtime_error);
    }
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

TEST(BioMD, Avx2GivesTheSameBytesAsScalar) {
    for (bool four : {false, true})
        for (size_t frames : {size_t(1), size_t(3)}) {
            SystemSpec s;
            s.four_site = four;
            s.frames = frames;
            s.waters = 1001;  // an odd count exercises the vector tails
            size_t n;
            auto x = make_system(s, &n);
            for (size_t i = s.chains * s.chain_len * 2; i + 2 < n; i += 70) x[3 * (i + 2) + 1] -= 2.5f;  // broken
            std::vector<uint8_t> a(SZ3::biomd::compress_bound(frames, n)), b(a.size());
            const size_t na = SZ3::biomd::compress(x.data(), frames, n, 5e-4, a.data(), false);
            const size_t nb = SZ3::biomd::compress(x.data(), frames, n, 5e-4, b.data(), true);
            ASSERT_EQ(na, nb);
            EXPECT_EQ(0, memcmp(a.data(), b.data(), na));
        }
}

TEST(BioMD, CorruptStreamsAreRefusedOrDecodeInBounds) {
    SystemSpec s;
    s.frames = 3;
    s.four_site = true;
    size_t n;
    const auto x = make_system(s, &n);
    std::vector<uint8_t> buf(SZ3::biomd::compress_bound(s.frames, n));
    const size_t size = SZ3::biomd::compress(x.data(), s.frames, n, 5e-4, buf.data());
    buf.resize(size);
    std::vector<float> out(x.size());
    std::mt19937 rng(9);
    auto attempt = [&](const std::vector<uint8_t> &b) {
        try {
            if (SZ3::biomd::stored_values(b.data(), b.size()) == out.size())
                SZ3::biomd::decompress(b.data(), b.size(), out.data());
        } catch (std::exception &) {
        }
    };
    for (size_t cut = 0; cut < size; cut += std::max<size_t>(1, size / 300)) {  // truncated
        std::vector<uint8_t> b(buf.begin(), buf.begin() + cut);
        attempt(b);
    }
    for (int k = 0; k < 300; k++) {  // bytes changed after the size fields
        auto b = buf;
        for (int m = 0; m < 4; m++) b[12 + rng() % (size - 12)] ^= uint8_t(1 + rng() % 255);
        attempt(b);
    }
    SUCCEED();
}

}  // namespace
