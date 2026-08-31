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
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include "SZ3/def.hpp"
#include "SZ3/utils/Config.hpp"
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

/// Fields that have historically broken decompositions: flat, ramp, single spike, extremes.
template <class T>
std::vector<std::vector<T>> adversarialFields(size_t num) {
    std::vector<std::vector<T>> out;
    out.push_back(std::vector<T>(num, T(3.5)));  // constant
    std::vector<T> ramp(num);
    for (size_t i = 0; i < num; i++) ramp[i] = static_cast<T>(i) * T(1e-3);
    out.push_back(ramp);
    std::vector<T> spike(num, T(0));
    spike[num / 2] = T(1e6);
    out.push_back(spike);
    std::vector<T> alternating(num);
    for (size_t i = 0; i < num; i++) alternating[i] = (i % 2) ? T(-1e4) : T(1e4);
    out.push_back(alternating);
    std::vector<T> smooth(num);
    for (size_t i = 0; i < num; i++)
        smooth[i] = std::sin(static_cast<double>(i) * 0.013) * 3.0 + std::cos(static_cast<double>(i) * 0.007);
    out.push_back(smooth);
    return out;
}

/**
 * @brief Acceptance checks for the Decomposition group.
 *
 * @tparam T       Data type
 * @tparam Factory Callable returning a fresh decomposition by value
 * @param name     Module name, printed on failure
 * @param conf     Configuration the factory's decomposition was built for
 * @param make     Factory
 * @param eb       Absolute error bound the module claims, or 0 to skip the bound check
 * @param fields   Input fields to exercise; adversarial defaults are used when empty
 */
template <class T, class Factory>
void expectDecompositionContract(const std::string &name, const SZ3::Config &conf, Factory make, double eb,
                                 std::vector<std::vector<T>> fields = {}) {
    if (fields.empty()) fields = adversarialFields<T>(conf.num);

    // 1. SZGenericCompressor requires the range to start at 0 and to fit in the int that
    //    `preprocess_encode` takes; 0 means "no bin range" and is allowed.
    {
        auto d = make();
        const auto range = d.get_out_range();
        EXPECT_EQ(range.first, 0) << name << ": get_out_range().first must be 0";
        EXPECT_TRUE(range.second == 0 || static_cast<int64_t>(range.second) <= std::numeric_limits<int>::max())
            << name << ": get_out_range().second must fit in int, or be 0 when there is no range";
    }

    for (size_t f = 0; f < fields.size(); f++) {
        SCOPED_TRACE(name + ": field " + std::to_string(f));
        const std::vector<T> &original = fields[f];
        ASSERT_EQ(original.size(), conf.num) << name << ": field size does not match conf.num";

        auto enc = make();
        std::vector<T> scratch = original;
        auto bins = enc.compress(conf, scratch.data());

        // 2. Every bin must sit inside the advertised range.
        const auto range = enc.get_out_range();
        if (range.second != 0) {
            for (size_t i = 0; i < bins.size(); i++) {
                ASSERT_GE(bins[i], range.first) << name << ": bin " << i << " below range";
                ASSERT_LE(bins[i], range.second) << name << ": bin " << i << " above range";
            }
        }

        // 3. size_est() must cover what save() writes, checked with guard bytes.
        const size_t est = enc.size_est();
        const size_t guard = 64;
        std::vector<SZ3::uchar> state(est + 4096 + guard);
        dirty(state);
        std::vector<SZ3::uchar> canary(state.end() - guard, state.end());
        SZ3::uchar *wp = state.data();
        enc.save(wp);
        const size_t saved = static_cast<size_t>(wp - state.data());
        ASSERT_LE(saved, state.size() - guard) << name << ": save() wrote past the buffer";
        EXPECT_EQ(std::memcmp(state.data() + state.size() - guard, canary.data(), guard), 0)
            << name << ": save() clobbered the guard bytes";
        EXPECT_LE(saved, est + 4096) << name << ": size_est() = " << est << " underestimates " << saved;

        // 4. A fresh instance loaded from that state must reconstruct within the bound. Going
        //    through save/load is what the compressor does, and it is where a decomposition that
        //    silently depends on its own compress-side state fails.
        auto dec = make();
        const SZ3::uchar *rp = state.data();
        size_t remaining = saved;
        dec.load(rp, remaining);
        std::vector<T> recovered(conf.num, T(0));
        dec.decompress(conf, bins, recovered.data());

        if (eb > 0) {
            double worst = 0;
            size_t worst_at = 0;
            for (size_t i = 0; i < conf.num; i++) {
                const double err = std::fabs(static_cast<double>(recovered[i]) - static_cast<double>(original[i]));
                if (err > worst) {
                    worst = err;
                    worst_at = i;
                }
            }
            EXPECT_LE(worst, eb * (1 + 1e-6))
                << name << ": |err| = " << worst << " exceeds eb = " << eb << " at index " << worst_at;
        }

        // 5. load() must reject a truncated stream rather than read past remaining_length.
        if (saved > 1) {
            auto trunc = make();
            const SZ3::uchar *tp = state.data();
            size_t short_len = saved - 1;
            EXPECT_ANY_THROW(trunc.load(tp, short_len)) << name << ": load() accepted a truncated stream";
        }

        // 6. Determinism.
        auto again = make();
        std::vector<T> scratch2 = original;
        auto bins2 = again.compress(conf, scratch2.data());
        ASSERT_EQ(bins2.size(), bins.size()) << name << ": compress() is not deterministic in length";
        for (size_t i = 0; i < bins.size(); i++) {
            ASSERT_EQ(bins2[i], bins[i]) << name << ": compress() is not deterministic at bin " << i;
        }
    }
}

}  // namespace SZ3_test

#endif  // SZ3_TEST_MODULE_CONTRACT_HPP
