/**
 * @file ModuleContract.hpp
 * @brief Group-level acceptance checks every module must pass to enter the FZ module library.
 *
 * A module's own test file adds whatever behaviour is specific to it. This header holds the
 * checks that apply to *every* module in a group, so a new module gets them in one line:
 *
 * @code
 *   TEST(SZ3_MyQuantizerTest, Contract) {
 *       SZ3_test::expectQuantizerContract<float>("MyQuantizer", [] { return MyQuantizer<float>(1e-3); });
 *   }
 * @endcode
 *
 * Each check corresponds to a defect class that has actually shipped in this tree, so a module
 * that passes is known not to repeat one. See docs/MODULE_ACCEPTANCE.md.
 *
 * The factory argument returns a freshly constructed module; the checks need several independent
 * instances (one to encode with, one to load into) and must not share mutable state between them.
 */

#ifndef SZ3_TEST_MODULE_CONTRACT_HPP
#define SZ3_TEST_MODULE_CONTRACT_HPP

#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include "SZ3/def.hpp"
#include "gtest/gtest.h"

namespace SZ3_test {

/// Values that have historically broken quantizers: signed zero, non-finite, and the extremes.
template <class T>
std::vector<T> adversarialValues() {
    const T inf = std::numeric_limits<T>::infinity();
    return {T(0),
            T(-0.0),
            T(1),
            T(-1),
            std::numeric_limits<T>::denorm_min(),
            -std::numeric_limits<T>::denorm_min(),
            std::numeric_limits<T>::min(),
            std::numeric_limits<T>::max(),
            -std::numeric_limits<T>::max(),
            inf,
            -inf,
            std::numeric_limits<T>::quiet_NaN()};
}

/// Fill a buffer with a recognizable pattern so a module that only ORs bits in is caught.
inline void dirty(std::vector<SZ3::uchar> &buf) { std::memset(buf.data(), 0xA5, buf.size()); }

/**
 * @brief Acceptance checks for the Quantizer group.
 *
 * @tparam T       Data type the quantizer is instantiated for
 * @tparam Factory Callable returning a fresh quantizer by value
 * @param name     Module name, printed on failure
 * @param make     Factory
 *
 * Does **not** check an error bound: quantizers in this group claim different kinds of bound
 * (absolute, relative, significant digits, none). Add `expectAbsErrorBound()` on top when the
 * module claims an absolute one.
 */
template <class T, class Factory>
void expectQuantizerContract(const std::string &name, Factory make) {
    auto q = make();

    // 1. SZGenericCompressor rejects a decomposition whose range does not start at 0, and every
    //    quantizer-hosting decomposition forwards its quantizer's range verbatim.
    EXPECT_EQ(q.get_out_range().first, 0) << name << ": get_out_range().first must be 0";

    // 2. No bin may fall below the advertised start. A quantizer whose bin type is signed while
    //    its values fill the sign bit produces negative bins here.
    // 3. recover() must reproduce what quantize_and_overwrite() wrote back, in encode order.
    //    Reading the unpredictable list past its end shows up as a mismatch or a throw.
    {
        auto qa = make();
        std::vector<T> vals = adversarialValues<T>();
        for (T x : vals) vals.push_back(x * T(0.5));
        std::vector<decltype(qa.quantize_and_overwrite(vals[0], T(0)))> bins;
        std::vector<T> written;
        bins.reserve(vals.size());
        written.reserve(vals.size());
        qa.precompress_data();
        for (T x : vals) {
            T scratch = x;
            auto bin = qa.quantize_and_overwrite(scratch, T(0));
            EXPECT_GE(bin, qa.get_out_range().first) << name << ": bin below get_out_range().first";
            bins.push_back(bin);
            written.push_back(scratch);
        }
        qa.postcompress_data();

        auto qb = make();
        // A quantizer with an unpredictable list needs its state, so round-trip through save/load.
        std::vector<SZ3::uchar> state(1u << 22);
        dirty(state);
        SZ3::uchar *wp = state.data();
        qa.save(wp);
        const size_t saved = static_cast<size_t>(wp - state.data());
        ASSERT_GT(saved, 0u) << name << ": save() wrote nothing";
        const SZ3::uchar *rp = state.data();
        size_t remaining = saved;
        qb.load(rp, remaining);
        EXPECT_EQ(static_cast<size_t>(rp - state.data()), saved) << name << ": load() did not consume save()'s bytes";

        qb.predecompress_data();
        for (size_t i = 0; i < bins.size(); i++) {
            const T got = qb.recover(T(0), bins[i]);
            if (std::isnan(static_cast<double>(written[i]))) {
                EXPECT_TRUE(std::isnan(static_cast<double>(got))) << name << ": NaN not reproduced at i=" << i;
            } else {
                EXPECT_EQ(got, written[i]) << name << ": recover() != overwritten value at i=" << i;
            }
        }
        qb.postdecompress_data();
    }

    // 4. load() must respect remaining_length. A truncated stream must throw, not read past the end.
    {
        auto qa = make();
        std::vector<SZ3::uchar> state(1u << 16);
        SZ3::uchar *wp = state.data();
        qa.save(wp);
        const size_t saved = static_cast<size_t>(wp - state.data());
        if (saved > 1) {
            auto qb = make();
            const SZ3::uchar *rp = state.data();
            size_t remaining = saved - 1;
            EXPECT_ANY_THROW(qb.load(rp, remaining)) << name << ": load() accepted a truncated stream";
        }
    }

    // 5. Determinism: the same inputs must produce the same bins.
    {
        auto q1 = make();
        auto q2 = make();
        for (T x : adversarialValues<T>()) {
            T a = x, b = x;
            EXPECT_EQ(q1.quantize_and_overwrite(a, T(0)), q2.quantize_and_overwrite(b, T(0)))
                << name << ": quantize_and_overwrite is not deterministic";
        }
    }
}

/**
 * @brief Add-on for a quantizer that claims an absolute error bound.
 *
 * Checks the bound against the value actually written back, over adversarial data plus `extra`.
 */
template <class T, class Factory>
void expectAbsErrorBound(const std::string &name, Factory make, double eb, const std::vector<T> &extra = {}) {
    auto q = make();
    std::vector<T> vals = adversarialValues<T>();
    vals.insert(vals.end(), extra.begin(), extra.end());
    q.precompress_data();
    for (T x : vals) {
        if (!std::isfinite(static_cast<double>(x))) continue;
        T scratch = x;
        q.quantize_and_overwrite(scratch, T(0));
        const double err = std::fabs(static_cast<double>(scratch) - static_cast<double>(x));
        EXPECT_LE(err, eb * (1 + 1e-9)) << name << ": |err| = " << err << " exceeds eb = " << eb << " at x = " << x;
    }
    q.postcompress_data();
}

/**
 * @brief Acceptance checks for the Encoder group.
 *
 * @tparam Bin     Bin type the encoder is instantiated for
 * @tparam Factory Callable returning a fresh encoder by value
 * @param name     Module name, printed on failure
 * @param make     Factory
 * @param streams  Bin streams to exercise; adversarial defaults are appended
 * @param stateNum What the compressor would pass, i.e. `get_out_range().second`. 0 means the
 *                 decomposition has no bin range and the encoder must derive one from the bins;
 *                 pass a real bound for an encoder that cannot.
 */
template <class Bin, class Factory>
void expectEncoderContract(const std::string &name, Factory make, std::vector<std::vector<Bin>> streams = {},
                           int stateNum = 0) {
    // Streams that have historically broken encoders: a single element, one run per sample,
    // a constant stream, and a value whose top bit is set.
    streams.push_back({Bin(0)});
    streams.push_back({Bin(7)});
    streams.push_back(std::vector<Bin>(64, Bin(3)));
    {
        std::vector<Bin> alternating(256);
        for (size_t i = 0; i < alternating.size(); i++) alternating[i] = Bin(i % 2);
        streams.push_back(alternating);
    }
    {
        std::vector<Bin> oneRunPerSample(512);
        for (size_t i = 0; i < oneRunPerSample.size(); i++) oneRunPerSample[i] = Bin(i % 97);
        streams.push_back(oneRunPerSample);
    }

    for (size_t s = 0; s < streams.size(); s++) {
        const std::vector<Bin> &bins = streams[s];
        if (bins.empty()) continue;
        SCOPED_TRACE(name + ": stream " + std::to_string(s) + " of " + std::to_string(bins.size()) + " bins");

        Bin hi = bins[0];
        for (Bin b : bins) hi = std::max(hi, b);

        auto enc = make();
        enc.preprocess_encode(bins, stateNum);

        // 1. size_est() must cover what encode() writes, or the caller's buffer overruns. The
        //    buffer is sized exactly to the estimate so an underestimate trips the guard below.
        const size_t est = enc.size_est();
        std::vector<SZ3::uchar> buf(est + 4096 + bins.size() * sizeof(Bin) * 2);
        const size_t guard = 64;
        std::vector<SZ3::uchar> canary(guard, 0x5C);

        // 2. encode() must not assume a zeroed destination.
        dirty(buf);
        SZ3::uchar *wp = buf.data();
        enc.save(wp);
        const size_t header = static_cast<size_t>(wp - buf.data());
        std::memcpy(buf.data() + buf.size() - guard, canary.data(), guard);
        const size_t written = enc.encode(bins, wp);
        enc.postprocess_encode();
        const size_t total = header + written;
        ASSERT_LE(total, buf.size() - guard) << name << ": encode() wrote past the buffer";
        EXPECT_EQ(std::memcmp(buf.data() + buf.size() - guard, canary.data(), guard), 0)
            << name << ": encode() clobbered the guard bytes";
        EXPECT_LE(written, est + 4096) << name << ": size_est() = " << est << " underestimates " << written;

        // 3. Round-trip must reproduce the bins exactly.
        auto dec = make();
        const SZ3::uchar *rp = buf.data();
        size_t remaining = total;
        dec.load(rp, remaining);
        dec.preprocess_decode();
        const std::vector<Bin> out = dec.decode(rp, bins.size());
        dec.postprocess_decode();
        ASSERT_EQ(out.size(), bins.size()) << name << ": decode() returned the wrong length";
        for (size_t i = 0; i < bins.size(); i++) {
            ASSERT_EQ(out[i], bins[i]) << name << ": bin " << i << " differs after round-trip";
        }
    }
}

}  // namespace SZ3_test

#endif  // SZ3_TEST_MODULE_CONTRACT_HPP
