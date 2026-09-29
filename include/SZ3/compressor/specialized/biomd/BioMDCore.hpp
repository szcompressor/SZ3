#ifndef SZ3_BIOMD_CORE_HPP
#define SZ3_BIOMD_CORE_HPP

// ALGO_BIOMD, part 1: lattice, bit I/O, Huffman coding, sphere / circle geometry and the detection of rigid water and
// bonds.
//
// Every coordinate is put on the lattice q = round(x / step), |x - q step| <= eb. All prediction runs on that lattice
// in integer arithmetic plus correctly rounded IEEE operations, and every rounding the decoder does goes through
// rnd_pred, so decoding gives the same integers on every IEEE machine and compiler, with or without FMA.
//
//  * rigid 3-site water (O, H, H): H1 on the sphere |H1 - O| = r (cube face, two kept coordinates, radial residual);
//    H2 on the circle fixed by |H2 - O| = r and the H-O-H angle (one kept coordinate, side bit, two residuals); the
//    virtual site of 4-site models from M = O + a (H1 + H2 - 2 O).
//  * other atoms: on the sphere of a bond-length class around one of the previous MAXOFF atoms, else as a delta.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>
#if defined(__x86_64__) || defined(_M_X64) || defined(__SSE2__) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2)
#include <emmintrin.h>
#define SZ3_BIOMD_SSE2_ROUND 1
#elif defined(__aarch64__) || defined(_M_ARM64)
#include <arm_neon.h>
#define SZ3_BIOMD_NEON_ROUND 1
#endif
#if defined(_MSC_VER) && !defined(__clang__)
#include <intrin.h>
#endif

#if defined(__GNUC__) || defined(__clang__)
#define SZ3_BIOMD_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define SZ3_BIOMD_INLINE __forceinline
#else
#define SZ3_BIOMD_INLINE inline
#endif
// GCC 13 at -O3 miscompiled Huffman table construction through ipa-modref when it was inlined (correct at -O2, with
// -fno-ipa-modref, or under sanitizers); the functions it hit carry this.
#if defined(__GNUC__) && !defined(__clang__)
#define SZ3_BIOMD_NOIPA __attribute__((noipa))
#else
#define SZ3_BIOMD_NOIPA
#endif
// The frame-wide loops and the kernels of BioMDSimd.hpp are also compiled for AVX2 and picked at run time on x86. Not
// with MinGW: its GCC does not align the stack for 32-byte AVX spills.
#if (defined(__x86_64__) || defined(_M_X64)) && (defined(__GNUC__) || defined(__clang__)) && !defined(__MINGW32__)
#define SZ3_BIOMD_X86_DISPATCH 1
#endif

namespace SZ3 {
namespace biomd {

// ------------------------------------------------------------------------------------------------ bit I/O
struct BitWriter {
    uint8_t *p;
    uint64_t acc = 0;
    int n = 0;
    explicit BitWriter(uint8_t *out) : p(out) {}
    inline void put(uint64_t v, int bits) {  // bits <= 32
        acc |= v << n;
        n += bits;
        if (n >= 32) {
            uint32_t w = uint32_t(acc);
            memcpy(p, &w, 4);
            p += 4;
            acc >>= 32;
            n -= 32;
        }
    }
    inline void put64(uint64_t v, int bits) {
        if (bits > 32) {
            put(v & 0xffffffffu, 32);
            put(v >> 32, bits - 32);
        } else {
            put(v, bits);
        }
    }
    uint8_t *finish() {
        while (n > 0) {
            *p++ = uint8_t(acc);
            acc >>= 8;
            n -= 8;
        }
        n = 0;
        acc = 0;
        return p;
    }
};

// Reads zeros past the end, so a truncated stream decodes to garbage but never reads out of bounds.
struct BitReader {
    const uint8_t *p, *end;
    uint64_t acc = 0;
    int n = 0;
    BitReader(const uint8_t *b, const uint8_t *e) : p(b), end(e) {}
    inline void refill() {
        if (end - p >= 4) {
            uint32_t w;
            memcpy(&w, p, 4);
            p += 4;
            acc |= uint64_t(w) << n;
            n += 32;
        } else {
            while (n <= 56) {
                uint64_t b = (p < end) ? *p : 0;
                p++;
                acc |= b << n;
                n += 8;
            }
        }
    }
    inline uint32_t peek(int bits) {
        if (n < bits) refill();
        return uint32_t(acc & ((1ull << bits) - 1));
    }
    inline void skip(int bits) {
        acc >>= bits;
        n -= bits;
    }
    inline uint32_t get(int bits) {
        if (!bits) return 0;
        uint32_t v = peek(bits);
        skip(bits);
        return v;
    }
    inline uint64_t get64(int bits) {
        if (bits > 32) {
            uint64_t lo = get(32);
            return lo | (uint64_t(get(bits - 32)) << 32);
        }
        return get(bits);
    }
    // back to the byte after the last bit read; past the end means the stream was truncated
    const uint8_t *align() {
        skip(n & 7);
        p -= n / 8;
        n = 0;
        acc = 0;
        if (p > end) throw std::runtime_error("SZ3 BioMD: truncated stream");
        return p;
    }
};

static inline uint32_t zz(int64_t v) { return uint32_t((uint64_t(v) << 1) ^ uint64_t(v >> 63)); }
static inline int64_t unzz(uint32_t u) { return int64_t(u >> 1) ^ -int64_t(u & 1); }
static inline int64_t rnd(double y) { return int64_t(y + std::copysign(0.5, y)); }
// The rounding of every value the decoder predicts. A product is converted to the nearest integer (ties to even) by
// one instruction, with no add a compiler could fuse with the multiply into an FMA: builds with and without FMA, and
// x86 and ARM, predict the same integers. Arguments stay well inside int32.
static inline int64_t rnd_pred(double y) {
#if defined(SZ3_BIOMD_SSE2_ROUND)
    return _mm_cvtsd_si32(_mm_set_sd(y));
#elif defined(SZ3_BIOMD_NEON_ROUND)
    return int64_t(vcvtnd_s64_f64(y));
#else
    volatile double v = y;  // the store rounds the product before nearbyint sees it
    return int64_t(std::nearbyint(v));
#endif
}

static inline int bit_width(uint32_t v) {  // v > 0
#if defined(_MSC_VER) && !defined(__clang__)
    unsigned long i;
    _BitScanReverse(&i, v);
    return int(i) + 1;
#else
    return 32 - __builtin_clz(v);
#endif
}

template <class T>
static inline void put_raw(uint8_t *&p, const T &v) {
    memcpy(p, &v, sizeof(T));
    p += sizeof(T);
}
template <class T>
static inline void get_raw(const uint8_t *&p, const uint8_t *end, T &v) {
    if (size_t(end - p) < sizeof(T)) throw std::runtime_error("SZ3 BioMD: truncated stream");
    memcpy(&v, p, sizeof(T));
    p += sizeof(T);
}
static inline void put_varint(uint8_t *&p, uint64_t v) {
    while (v >= 0x80) {
        *p++ = uint8_t(v | 0x80);
        v >>= 7;
    }
    *p++ = uint8_t(v);
}
static inline uint64_t get_varint(const uint8_t *&p, const uint8_t *end) {
    uint64_t v = 0;
    for (int s = 0;; s += 7) {
        if (p >= end || s > 63) throw std::runtime_error("SZ3 BioMD: truncated stream");
        const uint8_t b = *p++;
        v |= uint64_t(b & 0x7f) << s;
        if (!(b & 0x80)) return v;
    }
}

// ------------------------------------------------------------------------------------------------ Huffman
// Values below DIRECT are symbols; a larger value is a (bit length, next two bits) symbol followed by its low bits.
constexpr int HMAXLEN = 12;
constexpr uint32_t DIRECT = 1024;
constexpr uint32_t ALPHA = DIRECT + 22 * 4;
static inline uint32_t symof(uint32_t v) {
    if (v < DIRECT) return v;
    int nb = bit_width(v);
    return DIRECT + uint32_t(nb - 11) * 4 + ((v >> (nb - 3)) & 3);
}
static inline int rawbits(uint32_t s) { return s < DIRECT ? 0 : int((s - DIRECT) >> 2) + 8; }

// Length-limited Huffman code lengths for m >= 2 weights: sort packed (weight, index) keys, then the in-place
// Moffat-Katajainen algorithm; if the longest code exceeds HMAXLEN, flatten the weights and retry.
SZ3_BIOMD_NOIPA static void huff_lengths(uint64_t *f, uint32_t m, uint8_t *out) {
    std::vector<uint64_t> key(m);
    std::vector<uint32_t> A(m);
    for (;;) {
        for (uint32_t i = 0; i < m; i++) key[i] = (f[i] << 11) | i;
        std::sort(key.begin(), key.end());
        for (uint32_t i = 0; i < m; i++) A[i] = uint32_t(key[i] >> 11);
        uint32_t root = 0, leaf = 0;
        for (uint32_t next = 0; next + 1 < m; next++) {
            if (leaf >= m || (root < next && A[root] < A[leaf])) {
                A[next] = A[root];
                A[root++] = next;
            } else {
                A[next] = A[leaf++];
            }
            if (leaf >= m || (root < next && A[root] < A[leaf])) {
                A[next] += A[root];
                A[root++] = next;
            } else {
                A[next] += A[leaf++];
            }
        }
        A[m - 2] = 0;
        for (int32_t k = int32_t(m) - 3; k >= 0; k--) A[k] = A[A[k]] + 1;
        int32_t avail = 1, used = 0, depth = 0, r = int32_t(m) - 2, nx = int32_t(m) - 1;
        while (avail > 0) {
            while (r >= 0 && int32_t(A[r]) == depth) {
                used++;
                r--;
            }
            while (avail > used) {
                A[nx--] = uint32_t(depth);
                avail--;
            }
            avail = 2 * used;
            depth++;
            used = 0;
        }
        // A[i] is now the length of the i-th lightest symbol
        if (A[0] <= uint32_t(HMAXLEN)) {
            for (uint32_t i = 0; i < m; i++) out[key[i] & 2047] = uint8_t(A[i]);
            return;
        }
        for (uint32_t i = 0; i < m; i++) f[i] = (f[i] >> 1) | 1;
    }
}

struct Huff {
    uint32_t hist[ALPHA];
    uint32_t enc[ALPHA];  // code | len << 16
    uint8_t len[ALPHA];
    uint32_t maxsym = 0;
    bool used = false;
    int32_t single = -1;           // the only symbol of the stream (coded with zero bits), or -1
    std::vector<uint16_t> dtab;    // sym << 4 | len, indexed by the next HMAXLEN bits (or one period, see read_table)
    const uint16_t *dt = nullptr;  // dtab, or a table of zeros for a stream without symbols

    Huff() { reset(); }
    void reset() { memset(hist, 0, sizeof(hist)); }
    inline void count(uint32_t v) { hist[symof(v)]++; }

    void build() {
        uint32_t syms[ALPHA];
        uint32_t m = 0;
        maxsym = 0;
        for (uint32_t s = 0; s < ALPHA; s++)
            if (hist[s]) {
                syms[m++] = s;
                maxsym = s;
            }
        for (uint32_t s = 0; s <= maxsym; s++) len[s] = 0;
        used = m > 0;
        single = m == 1 ? int32_t(syms[0]) : -1;
        if (m == 0) return;
        if (m == 1) {
            enc[syms[0]] = 0;
            return;
        }
        uint64_t f[ALPHA];
        uint8_t l[ALPHA];
        for (uint32_t i = 0; i < m; i++) f[i] = hist[syms[i]];
        huff_lengths(f, m, l);
        for (uint32_t i = 0; i < m; i++) len[syms[i]] = l[i];
        make_codes();
    }
    SZ3_BIOMD_NOIPA void make_codes() {
        uint32_t blcount[HMAXLEN + 2] = {0}, next[HMAXLEN + 2] = {0};
        for (uint32_t s = 0; s <= maxsym; s++) blcount[len[s]]++;
        blcount[0] = 0;
        uint32_t c = 0;
        for (int b = 1; b <= HMAXLEN; b++) {
            c = (c + blcount[b - 1]) << 1;
            next[b] = c;
        }
        for (uint32_t s = 0; s <= maxsym; s++) {
            int l = len[s];
            if (!l) {
                enc[s] = 0;
                continue;
            }
            uint32_t v = next[l]++, r = 0;  // canonical code, bit-reversed for the LSB-first writer
            for (int i = 0; i < l; i++) {
                r = (r << 1) | (v & 1);
                v >>= 1;
            }
            enc[s] = r | (uint32_t(l) << 16);
        }
    }

    // Table format: used, single, then the code lengths up to maxsym. A length is coded against the previous one --
    // '0' same, '10' + sign one apart, '110' + 4 bits any other -- and a run of absent symbols as '111' + Elias
    // gamma(run).
    template <class W>
    void write_table(W &bw) const {
        bw.put(used, 1);
        if (!used) return;
        bw.put(single >= 0, 1);
        if (single >= 0) {
            bw.put(uint32_t(single), 11);
            return;
        }
        bw.put(maxsym, 11);
        int prev = 6;
        for (uint32_t s = 0; s <= maxsym;) {
            if (len[s] == 0) {
                uint32_t r = 0;
                while (s + r <= maxsym && len[s + r] == 0) r++;
                bw.put(7, 3);
                const int nb = bit_width(r);
                bw.put(uint32_t(1) << (nb - 1), nb);  // nb - 1 zeros, then a one (LSB first)
                if (nb > 1) bw.put(r & ((1u << (nb - 1)) - 1), nb - 1);
                s += r;
            } else {
                const int l = len[s];
                if (l == prev) {
                    bw.put(0, 1);
                } else if (l == prev + 1 || l == prev - 1) {
                    bw.put(1, 2);
                    bw.put(uint32_t(l > prev), 1);
                } else {
                    bw.put(3, 3);
                    bw.put(uint32_t(l), 4);
                }
                prev = l;
                s++;
            }
        }
    }
    // full = false: dtab holds one period of 2^(longest code) entries only, for a caller that masks the bits itself
    void read_table(BitReader &br, bool full = true) {
        static const uint16_t zeros[1u << HMAXLEN] = {0};
        dt = zeros;
        used = br.get(1);
        single = -1;
        if (!used) return;
        if (br.get(1)) {
            single = int32_t(br.get(11));
            if (single >= int32_t(ALPHA)) throw std::runtime_error("SZ3 BioMD: corrupt table");
            dtab.assign(1u << HMAXLEN, uint16_t(single << 4));
            dt = dtab.data();
            return;
        }
        maxsym = br.get(11);
        if (maxsym >= ALPHA) throw std::runtime_error("SZ3 BioMD: corrupt table");
        memset(len, 0, maxsym + 1);
        int prev = 6;
        for (uint32_t s = 0; s <= maxsym;) {
            if (br.get(1) == 0) {  // '0'
                len[s++] = uint8_t(prev);
                continue;
            }
            if (br.get(1) == 0) {  // '10' + sign
                prev += br.get(1) ? 1 : -1;
            } else if (br.get(1) == 0) {  // '110' + length
                prev = int(br.get(4));
            } else {  // '111' + gamma(run)
                int z = 0;
                while (br.get(1) == 0)
                    if (++z > 12) throw std::runtime_error("SZ3 BioMD: corrupt table");
                s += (1u << z) | (z ? br.get(z) : 0);
                continue;
            }
            if (prev < 1 || prev > HMAXLEN) throw std::runtime_error("SZ3 BioMD: corrupt table");
            len[s++] = uint8_t(prev);
        }
        uint64_t kraft = 0;  // a decodable prefix code has sum 2^-len <= 1
        for (uint32_t k = 0; k <= maxsym; k++)
            if (len[k]) kraft += uint64_t(1) << (HMAXLEN - len[k]);
        if (kraft > (uint64_t(1) << HMAXLEN)) throw std::runtime_error("SZ3 BioMD: corrupt table");
        make_codes();
        // entries repeat with period 2^(longest code): fill one period, then copy it
        int dbits = 1;
        for (uint32_t k = 0; k <= maxsym; k++) dbits = std::max(dbits, int(len[k]));
        const size_t per = size_t(1) << dbits;
        dtab.resize(full ? size_t(1) << HMAXLEN : per);
        std::fill(dtab.begin(), dtab.begin() + per, uint16_t(0));  // 0 (sym 0, len 0): bits no code starts with
        for (uint32_t s = 0; s <= maxsym; s++) {
            int l = len[s];
            if (!l) continue;
            for (uint32_t k = enc[s] & 0xffff; k < per; k += (1u << l)) dtab[k] = uint16_t((s << 4) | uint32_t(l));
        }
        if (full)
            for (size_t n = per; n < dtab.size(); n *= 2) memcpy(&dtab[n], &dtab[0], n * sizeof(uint16_t));
        dt = dtab.data();
    }
    inline void put(BitWriter &bw, uint32_t v) const {
        uint32_t s = symof(v), e = enc[s];
        bw.put(e & 0xffff, int(e >> 16));
        if (s >= DIRECT) {
            int nb = rawbits(s);
            bw.put(v & ((1u << nb) - 1), nb);
        }
    }
    inline uint32_t get(BitReader &br) const {
        uint32_t e = dt[br.peek(HMAXLEN)];
        br.skip(int(e & 15));
        uint32_t s = e >> 4;
        if (s >= DIRECT) {
            int nb = rawbits(s);
            s = ((4u | ((s - DIRECT) & 3)) << nb) | br.get(nb);
        }
        return s;
    }
};

// ------------------------------------------------------------------------------------------------ geometry
static inline int64_t isqrt_round(int64_t n) { return n <= 0 ? 0 : int64_t(std::sqrt(double(n)) + 0.5); }

// point on a sphere of squared radius R2 around the origin: cube face + 2 kept coords + radial residual
struct SphereCode {
    uint32_t face;
    int64_t a, b, e;
};
static inline SphereCode sphere_encode(const int64_t d[3], int64_t R2) {
    int f = 0;
    int64_t m = std::llabs(d[0]), m1 = std::llabs(d[1]), m2 = std::llabs(d[2]);
    if (m1 > m) {
        f = 1;
        m = m1;
    }
    if (m2 > m) {
        f = 2;
        m = m2;
    }
    int i = f == 2 ? 0 : f + 1, j = f == 0 ? 2 : (f == 1 ? 0 : 1);
    SphereCode c;
    c.face = uint32_t(f * 2 + (d[f] < 0));
    c.a = d[i];
    c.b = d[j];
    c.e = m - isqrt_round(R2 - c.a * c.a - c.b * c.b);
    return c;
}
static inline void sphere_decode(uint32_t face, int64_t a, int64_t b, int64_t e, int64_t R2, int64_t d[3]) {
    int f = int(face >> 1) % 3;
    int i = f == 2 ? 0 : f + 1, j = f == 0 ? 2 : (f == 1 ? 0 : 1);
    d[i] = a;
    d[j] = b;
    int64_t m = isqrt_round(R2 - a * a - b * b) + e;
    d[f] = (face & 1) ? -m : m;
}

// H2 of rigid water, relative to O, given u = H1 - O. Keep coordinate k = argmin|u_k| (relative to the circle
// centre's k coordinate); the other two follow from h.u = P and |h|^2 = R2 (two roots, side bit).
struct Circle {
    int k, i, j;
    int64_t ck, P, al, be;
    double invA;  // 1 / (al^2 + be^2), 0 for a molecule that is not intact
};
// |u| is within ~1 lattice unit of R, so |u| ~ (|u|^2 + R^2) / (2R) to ~1e-2 units: P = Rc |u| needs no sqrt.
// Rc2R = Rc / (2R), invR2 = 1 / R^2.
static inline void circle_setup(const int64_t u[3], double Rc2R, int64_t R2, double invR2, Circle &c) {
    int64_t uu = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    if (uu > 4 * R2) {  // broken molecule: no prediction (circle_solve sees invA == 0)
        c.k = 0;
        c.i = 1;
        c.j = 2;
        c.P = 0;
        c.ck = 0;
        c.al = u[1];
        c.be = u[2];
        c.invA = 0;
        return;
    }
    c.P = rnd_pred(Rc2R * double(uu + R2));
    int k = 0;
    int64_t a0 = std::llabs(u[0]), a1 = std::llabs(u[1]), a2 = std::llabs(u[2]);
    if (a1 < a0) {
        k = 1;
        a0 = a1;
    }
    if (a2 < a0) k = 2;
    c.k = k;
    c.i = k == 2 ? 0 : k + 1;
    c.j = k == 0 ? 2 : (k == 1 ? 0 : 1);
    c.ck = rnd_pred(double(c.P * u[k]) * invR2);  // |u|^2 ~ R2; the <0.5 unit error only shifts the kept symbol
    c.al = u[c.i];
    c.be = u[c.j];
    const int64_t A = c.al * c.al + c.be * c.be;
    c.invA = A > 0 ? 1.0 / double(A) : 0.0;
}
static inline void circle_solve(const int64_t u[3], int64_t R2, const Circle &c, int64_t hk, int64_t x[2],
                                int64_t y[2]) {
    // a molecule that is not intact (|H1-O| or |h_k| beyond 2R, e.g. split over a periodic boundary): no prediction
    if (c.invA == 0 || hk * hk > 4 * R2 || (c.al * c.al + c.be * c.be) > 4 * R2) {
        x[0] = x[1] = y[0] = y[1] = 0;
        return;
    }
    int64_t L = c.P - u[c.k] * hk;
    int64_t S = R2 - hk * hk;
    int64_t sq = isqrt_round(S * (c.al * c.al + c.be * c.be) - L * L);
    int64_t la = L * c.al, lb = L * c.be, bs = c.be * sq, as = c.al * sq;
    x[0] = rnd_pred(double(la + bs) * c.invA);
    y[0] = rnd_pred(double(lb - as) * c.invA);
    x[1] = rnd_pred(double(la - bs) * c.invA);
    y[1] = rnd_pred(double(lb + as) * c.invA);
}

// ------------------------------------------------------------------------------------------------ frame-wide loops
// Written once, compiled for the baseline and (on x86) for AVX2, picked at run time. All give the same results.
static inline bool cpu_avx2() {
#if defined(SZ3_BIOMD_X86_DISPATCH)
    static const bool ok = __builtin_cpu_supports("avx2");
    return ok;
#else
    return false;
#endif
}

template <class T>
using bits_of = typename std::conditional<sizeof(T) == 4, uint32_t, uint64_t>::type;

// the largest |x| as bits (NaN and inf above every finite value)
template <class T>
SZ3_BIOMD_INLINE bits_of<T> max_abs_bits_body(const T *x, size_t n) {
    using U = bits_of<T>;
    const U absmask = U(~U(0)) >> 1;
    U m = 0;
    for (size_t i = 0; i < n; i++) {
        U b;
        memcpy(&b, &x[i], sizeof(T));
        b &= absmask;
        m = b > m ? b : m;
    }
    return m;
}
// q = round(x / step)
template <class T>
SZ3_BIOMD_INLINE void quantize_body(const T *x, size_t n, double inv, int32_t *q) {
    for (size_t i = 0; i < n; i++) {
        const double y = double(x[i]) * inv;
        q[i] = int32_t(y + std::copysign(0.5, y));
    }
}
// the same, and true if a value is more than slack from the lattice
template <class T>
SZ3_BIOMD_INLINE bool quantize_checked_body(const T *x, size_t n, double inv, double slack, int32_t *q) {
    int bad = 0;
    for (size_t i = 0; i < n; i++) {
        const double y = double(x[i]) * inv;
        const int32_t qi = int32_t(y + std::copysign(0.5, y));
        q[i] = qi;
        bad |= int(std::fabs(y - double(qi)) > slack);
    }
    return bad != 0;
}
constexpr int MAXOFF = 4;   // a bond partner is one of the previous MAXOFF atoms
constexpr int MAXCLS = 15;  // bond-length classes
#if defined(SZ3_BIOMD_X86_DISPATCH)
template <class T>
__attribute__((target("avx2"))) static bits_of<T> max_abs_bits_avx2(const T *x, size_t n) {
    return max_abs_bits_body(x, n);
}
template <class T>
__attribute__((target("avx2"))) static void quantize_avx2(const T *x, size_t n, double inv, int32_t *q) {
    quantize_body(x, n, inv, q);
}
template <class T>
__attribute__((target("avx2"))) static bool quantize_checked_avx2(const T *x, size_t n, double inv, double slack,
                                                                  int32_t *q) {
    return quantize_checked_body(x, n, inv, slack, q);
}
#endif
template <class T>
static bits_of<T> max_abs_bits(const T *x, size_t n) {
#if defined(SZ3_BIOMD_X86_DISPATCH)
    if (cpu_avx2()) return max_abs_bits_avx2(x, n);
#endif
    return max_abs_bits_body(x, n);
}
template <class T>
static void quantize(const T *x, size_t n, double inv, int32_t *q) {
#if defined(SZ3_BIOMD_X86_DISPATCH)
    if (cpu_avx2()) return quantize_avx2(x, n, inv, q);
#endif
    quantize_body(x, n, inv, q);
}
template <class T>
static bool quantize_checked(const T *x, size_t n, double inv, double slack, int32_t *q) {
#if defined(SZ3_BIOMD_X86_DISPATCH)
    if (cpu_avx2()) return quantize_checked_avx2(x, n, inv, slack, q);
#endif
    return quantize_checked_body(x, n, inv, slack, q);
}

// ------------------------------------------------------------------------------------------------ layout
struct Layout {
    // per atom: 0 = water O (atoms i+1 .. i+nsite-1 are its H, H and, for 4-site models, M), 1 = the rest of a water
    // (coded with its O), 2 = other
    std::vector<uint8_t> kind;
    int nsite = 3;  // 3: O,H,H   4: O,H,H,M with M = O + a (H1 - O + H2 - O)  (TIP4P-style virtual site)
    double vs_a = 0;
    // other atoms: 0 = no bond, else off * 16 + cls + 1: on the sphere of class cls around atom i - off
    std::vector<uint16_t> ref;
    double r = 0, rhh = 0;     // water O-H and H-H (nm)
    std::vector<double> blen;  // bond-length classes (nm)
};

// Rigid water on frame x: kind, nsite, vs_a, r, rhh. The tolerances grow with the lattice step, so rounded input (xtc
// files, or data this codec decompressed) still fits.
template <class T>
static void detect_water(const T *x, size_t N, Layout &L, double step) {
    L.kind.assign(N, 2);
    L.r = L.rhh = 0;
    L.nsite = 3;
    L.vs_a = 0;
    auto d2 = [x](size_t a, size_t b) {
        float dx = float(x[3 * a] - x[3 * b]), dy = float(x[3 * a + 1] - x[3 * b + 1]),
              dz = float(x[3 * a + 2] - x[3 * b + 2]);
        return dx * dx + dy * dy + dz * dz;
    };
    // Candidates: O-H in [0.08, 0.125] nm, H-H in [0.13, 0.2] nm. The O-H / H-H modes are estimated from ~2000
    // candidates probed across the frame; rigid water gives a sharp H-H peak, flexible CH2/NH2 groups a broad one.
    auto cand = [&](size_t i) {
        float a = d2(i, i + 1), b = d2(i, i + 2);
        if (!(a > 0.0064f && a < 0.015625f && b > 0.0064f && b < 0.015625f)) return false;
        float c = d2(i + 1, i + 2);
        return c > 0.0169f && c < 0.04f;
    };
    std::vector<float> cand_oh, cand_hh;
    if (N >= 3) {
        const size_t probe = std::max<size_t>(1, N / 300);
        for (size_t i = 0; i + 2 < N; i += probe)
            for (size_t o = 0; o < 3 && i + o + 2 < N; o++)
                if (cand(i + o)) {
                    cand_hh.push_back(std::sqrt(d2(i + o + 1, i + o + 2)));
                    cand_oh.push_back(std::sqrt(d2(i + o, i + o + 1)));
                    cand_oh.push_back(std::sqrt(d2(i + o, i + o + 2)));
                    break;
                }
    }
    if (cand_hh.size() < 16) return;
    // tolerances stay well below the O-H / C-H difference (> 0.01 nm) that separates water from CH2/NH2
    const double tight = std::max(0.001, step), tol = std::max(0.002, 2.0 * step), loose = std::max(0.01, 3.0 * tight);
    auto wmode = [](const std::vector<float> &v, float lo, size_t nb) {  // centre of the heaviest 0.004 nm window
        std::vector<uint32_t> h(nb, 0);
        for (float y : v) {
            long k = long((y - lo) * 1e4f);
            if (k >= 0 && k < long(nb)) h[k]++;
        }
        size_t m = 0;
        uint64_t b = 0, s = 0;
        for (size_t k = 0; k < nb; k++) {
            s += h[k];
            if (k >= 41) s -= h[k - 41];
            if (s > b) {
                b = s;
                m = k;
            }
        }
        return double(lo) + (double(m >= 20 ? m - 20 : 0) + 0.5) * 1e-4;
    };
    auto refine = [](const std::vector<float> &v, double c) {  // median of the values within 0.003 of c
        std::vector<float> w;
        for (float y : v)
            if (std::fabs(y - c) < 0.003) w.push_back(y);
        if (w.empty()) return c;
        std::nth_element(w.begin(), w.begin() + w.size() / 2, w.end());
        return double(w[w.size() / 2]);
    };
    const double h0 = refine(cand_hh, wmode(cand_hh, 0.13f, 710));
    std::vector<float> ohs;
    for (size_t c = 0; c < cand_hh.size(); c++)
        if (std::fabs(cand_hh[c] - h0) < tol) {
            ohs.push_back(cand_oh[2 * c]);
            ohs.push_back(cand_oh[2 * c + 1]);
        }
    const double r0 = refine(ohs, wmode(ohs, 0.08f, 460));
    // rigidity test: rigid water puts nearly all nearby candidates within 0.001 nm of both modes
    size_t ntight = 0, nloose = 0;
    for (size_t c = 0; c < cand_hh.size(); c++) {
        double dh = std::fabs(cand_hh[c] - h0), da = std::fabs(cand_oh[2 * c] - r0),
               db = std::fabs(cand_oh[2 * c + 1] - r0);
        if (dh < loose && da < loose && db < loose) {
            nloose++;
            ntight += dh < tight && da < tight && db < tight;
        }
    }
    if (nloose < 16 || ntight * 10 < nloose * 6) return;
    const float alo = float((r0 - tol) * (r0 - tol)), ahi = float((r0 + tol) * (r0 + tol));
    const float clo = float((h0 - tol) * (h0 - tol)), chi = float((h0 + tol) * (h0 + tol));
    size_t nw = 0;
    for (size_t i = 0; i + 2 < N;) {
        float a = d2(i, i + 1), b, c;
        if (a > alo && a < ahi && (b = d2(i, i + 2)) > alo && b < ahi && (c = d2(i + 1, i + 2)) > clo && c < chi) {
            L.kind[i] = 0;
            L.kind[i + 1] = L.kind[i + 2] = 1;
            nw++;
            i += 3;
        } else {
            i++;
        }
    }
    if (nw < 16) {
        std::fill(L.kind.begin(), L.kind.end(), 2);
        return;
    }
    // 4-site water? The atom after each O,H,H is then M = O + a (H1 + H2 - 2 O), a fitted on a sample.
    double sxy = 0, sxx = 0;
    std::vector<size_t> samp;
    for (size_t i = 0; i + 3 < N && samp.size() < 512 && (i < 30000 || !samp.empty()); i++)
        if (L.kind[i] == 0) {
            if (L.kind[i + 3] == 2) samp.push_back(i);
            i += 2;
        }
    for (size_t i : samp)
        for (int c = 0; c < 3; c++) {
            double sv = double(x[3 * (i + 1) + c]) + x[3 * (i + 2) + c] - 2.0 * x[3 * i + c],
                   m = double(x[3 * (i + 3) + c]) - x[3 * i + c];
            sxy += sv * m;
            sxx += sv * sv;
        }
    const double a = sxx > 0 ? sxy / sxx : 0;
    auto is_vs = [&](size_t i) {
        for (int c = 0; c < 3; c++) {
            double pr = x[3 * i + c] + a * (double(x[3 * (i + 1) + c]) + x[3 * (i + 2) + c] - 2.0 * x[3 * i + c]);
            if (std::fabs(pr - x[3 * (i + 3) + c]) > 0.002) return false;
        }
        return true;
    };
    size_t good = 0;
    for (size_t i : samp) good += is_vs(i);
    if (samp.size() >= 16 && a > 0.01 && a < 0.5 && good * 10 >= samp.size() * 9) {
        L.nsite = 4;
        L.vs_a = a;
        for (size_t i = 0; i < N; i++)
            if (L.kind[i] == 0 && !(i + 3 < N && L.kind[i + 3] == 2 && is_vs(i)))
                L.kind[i] = L.kind[i + 1] = L.kind[i + 2] = 2;  // not a complete 4-site water: other atoms
        for (size_t i = 0; i < N; i++)
            if (L.kind[i] == 0) L.kind[i + 3] = 1;
    }
    // geometry from the tight candidates (means are exact for rigid water, ~unbiased for rounded input)
    double sr = 0, sh = 0;
    size_t ns = 0;
    for (size_t c = 0; c < cand_hh.size(); c++)
        if (std::fabs(cand_hh[c] - h0) < tol && std::fabs(cand_oh[2 * c] - r0) < tol &&
            std::fabs(cand_oh[2 * c + 1] - r0) < tol) {
            sr += double(cand_oh[2 * c]) + cand_oh[2 * c + 1];
            sh += cand_hh[c];
            ns++;
        }
    L.r = sr / (2.0 * ns);
    L.rhh = sh / ns;
}

// Bonds of the other atoms (kind 2) on frame x: each takes the nearest of its previous MAXOFF atoms if that is within
// bond range, and the bond-length classes are the modes of those distances.
template <class T>
static void detect_bonds(const T *x, size_t N, Layout &L) {
    L.ref.assign(N, 0);
    L.blen.clear();
    // Each atom of a run of other atoms takes the nearest of its previous MAXOFF atoms (in or before the run):
    // offset dof, and the 1e-4 nm bin of the distance over 0.05 .. 0.25 nm if it is in bond range, else NB.
    const size_t NB = 2000;
    thread_local std::vector<uint8_t> dof;
    thread_local std::vector<uint16_t> bin;
    dof.assign(N, 0);
    bin.assign(N, uint16_t(NB));
    for (size_t i = 1; i < N;) {
        if (L.kind[i] != 2) {
            i++;
            continue;
        }
        for (; i < N && L.kind[i] == 2; i++) {
            float best = 1e30f;
            uint8_t bo = 0;
            for (size_t o = 1; o <= size_t(MAXOFF) && o <= i; o++) {  // branchless: the minimum is unpredictable
                const float dx = float(x[3 * i] - x[3 * (i - o)]), dy = float(x[3 * i + 1] - x[3 * (i - o) + 1]),
                            dz = float(x[3 * i + 2] - x[3 * (i - o) + 2]);
                const float v = dx * dx + dy * dy + dz * dz;
                const bool lt = v < best;
                best = lt ? v : best;
                bo = lt ? uint8_t(o) : bo;
            }
            dof[i] = bo;
            if (best > 0.0025f && best < 0.0625f)
                bin[i] = uint16_t(std::min(NB - 1, size_t((std::sqrt(best) - 0.05f) * 1e4f)));
        }
    }
    std::vector<uint32_t> hist(NB, 0);
    for (size_t i = 0; i < N; i++)
        if (bin[i] < NB) hist[bin[i]]++;
    // peaks: take the heaviest +-10-bin window, suppress +-30 bins, up to MAXCLS classes; the window sums are
    // computed once and updated as bins are suppressed
    std::vector<uint32_t> hs = hist;
    std::vector<int64_t> win(NB, 0);
    int64_t remaining = 0;
    {
        int64_t run = 0;
        for (size_t k = 0; k < NB; k++) remaining += hs[k];
        for (size_t k = 0; k < std::min<size_t>(NB, 11); k++) run += hs[k];
        for (size_t k = 0; k < NB; k++) {
            win[k] = run;
            if (k + 11 < NB) run += hs[k + 11];
            if (k >= 10) run -= hs[k - 10];
        }
    }
    for (int c = 0; c < MAXCLS; c++) {
        size_t best = 0;
        int64_t bv = 0;
        for (size_t k = 0; k < NB; k++)
            if (win[k] > bv) {
                bv = win[k];
                best = k;
            }
        if (bv < 8 || bv * 200 < remaining + bv) break;  // ignore classes under 0.5% of the bonded atoms
        double sw = 0, sm = 0;
        for (size_t w = (best >= 10 ? best - 10 : 0); w < std::min(NB, best + 11); w++) {
            sw += hist[w];
            sm += hist[w] * (0.05 + (w + 0.5) * 1e-4);
        }
        L.blen.push_back(sm / sw);
        for (size_t w = (best >= 30 ? best - 30 : 0); w < std::min(NB, best + 31); w++) {
            if (!hs[w]) continue;
            for (size_t k = (w >= 10 ? w - 10 : 0); k < std::min(NB, w + 11); k++) win[k] -= hs[w];
            remaining -= hs[w];
            hs[w] = 0;
        }
    }
    // bin -> nearest class within 0.004 nm (40 bins); bin NB: none
    std::vector<int8_t> bincls(NB + 1, -1);
    std::vector<uint8_t> bdist(NB, 255);
    for (size_t c = 0; c < L.blen.size(); c++) {
        long centre = long((L.blen[c] - 0.05) * 1e4);
        for (long k = std::max(0L, centre - 40); k <= std::min(long(NB) - 1, centre + 40); k++) {
            uint8_t dd = uint8_t(std::labs(k - centre));
            if (dd < bdist[k]) {
                bdist[k] = dd;
                bincls[k] = int8_t(c);
            }
        }
    }
    for (size_t i = 0; i < N; i++) {
        const int bc = bincls[bin[i]];
        L.ref[i] = bc >= 0 ? uint16_t(dof[i] * 16 + bc + 1) : 0;
    }
}

}  // namespace biomd
}  // namespace SZ3
#endif
