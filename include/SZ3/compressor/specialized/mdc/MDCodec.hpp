#ifndef SZ3_MDC_CODEC_HPP
#define SZ3_MDC_CODEC_HPP

// ALGO_MDC, part 3: the stream. One codec for one-frame and multi-frame chunks of molecular-dynamics coordinates
// {frames, atoms, 3}.
//
// The atoms fall into four groups, and every frame picks, per group, the predictor that is cheapest on a sample of
// that group (frame 0 of a chunk, and so every one-frame chunk, is all intra):
//
//   water O      0 intra: box lattice, mixed radix    1 x[t-1]        2 2x[t-1] - x[t-2]
//   water H/M    0 intra: sphere + circle geometry    1 O->H vectors of t-1, then the same sphere/circle constraints
//                                                     2 O->H vectors extrapolated from t-1, t-2, then the constraints
//   bonded       0 intra: bond sphere around parent   1 parent->atom vector of t-1 + bond sphere   2 extrapolated
//                3 x[t-1] (for 'bonds' that are not: molecules split over the periodic boundary)
//   unbonded     0 x[i-1] of this frame   1 x[t-1]   2 2x[t-1] - x[t-2]   3 x[t-1] + (x[i-1] - x[i-1] at t-1)
//
// The predicted vector of a temporal mode also fixes the dropped axis and its sign, so those spend no face symbol.
//
// Stream: header (atoms, frames, trailing fill frames, lattice step, water geometry, bond-length classes), the water
// runs, the bond references of the other atoms (Huffman coded in the context of the previous atom's), then per frame
// the group modes, the Huffman tables (or "keep the previous frame's") and the symbols.

#include "SZ3/compressor/specialized/mdc/MDCCore.hpp"
#include "SZ3/compressor/specialized/mdc/MDCSimd.hpp"

namespace SZ3 {
namespace mdc {

enum { G_O, G_WH, G_NB, G_NU, NGROUP };
// S_FACE / S_BF: cube face and radial residual in one symbol; S_RS: the two circle residuals and the side bit;
// S_MJ: the three virtual-site residuals. Values too large for a joint symbol escape to S_E / S_BE, S_R, S_M.
enum { S_O, S_FACE, S_KEPT, S_E, S_V, S_R, S_M, S_BF, S_BK, S_BE, S_U, S_RS, S_MJ, NS };
constexpr int NMODE[NGROUP] = {3, 3, 4, 4};
constexpr uint32_t MAGIC = 0x3143444d;  // "MDC1"
constexpr int64_t PCLAMP = int64_t(1) << 29;
constexpr int64_t RMAXB2 = int64_t(1) << 28;  // bonds of 16384 lattice units or more are not coded as spheres
constexpr size_t MAXCTX = 63;                 // contexts of the reference coding
static_assert(16 * (MAXOFF + 1) <= 256, "reference values must fit a byte");

struct Options {
    int simd = SIMD_AVX512;  // widest instruction set used: SIMD_SCALAR .. SIMD_AVX512 (all give the same bytes)
};

static inline int64_t clampp(int64_t v) { return v > PCLAMP ? PCLAMP : (v < -PCLAMP ? -PCLAMP : v); }

// sphere point d predicted by p: p fixes the dropped axis and its sign; the two kept coordinates are residuals
// against p, the dropped one is sphere-reconstructed and corrected by e.
static inline void sph_axes(const int64_t p[3], int &f, int &i, int &j, int64_t &s) {
    f = 0;
    int64_t m = std::llabs(p[0]);
    if (std::llabs(p[1]) > m) {
        f = 1;
        m = std::llabs(p[1]);
    }
    if (std::llabs(p[2]) > m) f = 2;
    i = f == 2 ? 0 : f + 1;
    j = f == 0 ? 2 : (f == 1 ? 0 : 1);
    s = p[f] < 0 ? -1 : 1;
}
static inline void sph_pred_enc(const int64_t d[3], const int64_t p[3], int64_t R2, int64_t &ra, int64_t &rb,
                                int64_t &e) {
    int f, i, j;
    int64_t s;
    sph_axes(p, f, i, j, s);
    ra = d[i] - p[i];
    rb = d[j] - p[j];
    e = s * d[f] - isqrt_round(R2 - d[i] * d[i] - d[j] * d[j]);
}
static inline void sph_pred_dec(const int64_t p[3], int64_t R2, int64_t ra, int64_t rb, int64_t e, int64_t d[3]) {
    int f, i, j;
    int64_t s;
    sph_axes(p, f, i, j, s);
    d[i] = p[i] + ra;
    d[j] = p[j] + rb;
    d[f] = s * (isqrt_round(R2 - d[i] * d[i] - d[j] * d[j]) + e);
}

// Huffman symbol k of stream s, or (s == 255) nb raw bits
struct Tok {
    uint32_t v;
    uint16_t k;
    uint8_t s;
    uint8_t nb;
};
// writes tokens into a buffer sized for the worst case (12 per atom) and counts the symbols
struct StoreSink {
    Tok *t;
    Huff *H;
    inline void sym(int s, uint32_t v) {
        const uint32_t k = symof(v);
        *t++ = {v, uint16_t(k), uint8_t(s), 0};
        H[s].hist[k]++;
    }
    inline void raw(uint64_t v, int nb) {  // up to 48 bits per token (the high 16 in k)
        while (nb > 48) {
            *t++ = {uint32_t(v), 0, 255, 32};
            v >>= 32;
            nb -= 32;
        }
        *t++ = {uint32_t(v), uint16_t(v >> 32), 255, uint8_t(nb)};
    }
};
// cost estimate on a sample: empirical entropy of each stream + raw bits
struct EstSink {
    std::vector<uint32_t> *hist;     // NS histograms of ALPHA
    std::vector<uint32_t> *touched;  // NS lists of nonzero bins
    double rawb = 0;
    inline void sym(int s, uint32_t v) {
        uint32_t k = symof(v);
        if (hist[s][k]++ == 0) touched[s].push_back(k);
        rawb += rawbits(k);
    }
    inline void raw(uint64_t, int nb) { rawb += nb; }
    // Estimated bits for the whole group (scale = units / sampled units); resets the histograms. The sample entropy
    // gets the Miller-Madow correction, (K - 1) / (2 ln 2) bits per stream for K distinct symbols: without it a small
    // sample favours wide distributions. The table costs ~4 bits per symbol once per frame, not per unit.
    double take(double scale) {
        double b = rawb * scale;
        rawb = 0;
        for (int s = 0; s < NS; s++) {
            if (touched[s].empty()) continue;
            double n = 0, h = 0;
            for (uint32_t k : touched[s]) n += hist[s][k];
            for (uint32_t k : touched[s]) {
                double c = hist[s][k];
                h -= c * std::log2(c / n);
                hist[s][k] = 0;
            }
            const double K = double(touched[s].size());
            b += scale * (h + (K - 1) / (2 * 0.6931471805599453)) + 4.0 * K;
            touched[s].clear();
        }
        return b;
    }
};
template <class Sink>
static inline void put_fe(Sink &k, int s, int se, uint32_t face, uint32_t ze) {
    if (ze < 15) {
        k.sym(s, face * 16 + ze);
    } else {
        k.sym(s, face * 16 + 15);
        k.sym(se, ze);
    }
}
template <class Sink>
static inline void put_rs(Sink &k, uint32_t r1, uint32_t r2, uint32_t side) {
    if (r1 < 7 && r2 < 7) {
        k.sym(S_RS, ((r1 * 7 + r2) << 1) | side);
    } else {
        k.sym(S_RS, 98 + side);
        k.sym(S_R, r1);
        k.sym(S_R, r2);
    }
}
template <class Sink>
static inline void put_m3(Sink &k, const uint32_t z[3]) {
    if (z[0] < 5 && z[1] < 5 && z[2] < 5) {
        k.sym(S_MJ, (z[0] * 5 + z[1]) * 5 + z[2]);
    } else {
        k.sym(S_MJ, 125);
        for (int c = 0; c < 3; c++) k.sym(S_M, z[c]);
    }
}

// Bits for a mixed-radix index below Rx Ry Rz, and whether it fits 62 bits (else per-axis widths).
static inline void box_bits(uint64_t Rx, uint64_t Ry, uint64_t Rz, int &obits, int ob[3], bool &mixed) {
    auto bl = [](uint64_t r) {
        int b = 0;
        while (b < 64 && (uint64_t(1) << b) < r) b++;
        return b;
    };
    ob[0] = bl(Rx);
    ob[1] = bl(Ry);
    ob[2] = bl(Rz);
    obits = 0;
    mixed = ob[0] + ob[1] + ob[2] <= 64;
    if (mixed) {  // then Rx Ry fits 64 bits unless it is 2^64, and Rx Ry Rz is below 2^64 or its bit length is > 62
        const uint64_t xy = Rx * Ry;
        const bool xy_ovf = Ry != 0 && xy / Ry != Rx;
        const uint64_t xyz = xy * Rz;
        if (xy_ovf || (Rz != 0 && xyz / Rz != xy)) {
            mixed = false;
        } else {
            while (obits < 64 && (uint64_t(1) << obits) < xyz) obits++;
            mixed = obits <= 62;
        }
    }
    if (mixed) ob[0] = ob[1] = ob[2] = 0;
}

struct FrameCtx {
    const int32_t *q, *qp, *qpp;
    const Layout *L;
    int64_t R2;
    double Rc2R, invR2;  // see circle_setup
    const std::vector<int64_t> *BR2;
    // water O box
    int32_t omin[3];
    uint64_t Rx, Ry, Rz;
    int obits, ob[3];
    bool mixed;
};

static inline void rel(const int32_t *a, size_t i, size_t j, int64_t o[3]) {
    for (int c = 0; c < 3; c++) o[c] = int64_t(a[3 * i + c]) - a[3 * j + c];
}

// ------------------------------------------------------------------------------------------------ encoders
// Sink = StoreSink (the frame) or EstSink (the predictor choice)
template <class Sink>
static inline void enc_O(const FrameCtx &C, int m, size_t i, Sink &k) {
    const int32_t *O = &C.q[3 * i];
    if (m == 0) {
        if (C.mixed)
            k.raw(uint64_t(O[0] - C.omin[0]) + C.Rx * (uint64_t(O[1] - C.omin[1]) + C.Ry * uint64_t(O[2] - C.omin[2])),
                  C.obits);
        else
            for (int c = 0; c < 3; c++) k.raw(uint64_t(int64_t(O[c]) - C.omin[c]), C.ob[c]);
    } else {
        for (int c = 0; c < 3; c++) {
            int64_t pr = m == 1 ? C.qp[3 * i + c] : clampp(2 * int64_t(C.qp[3 * i + c]) - C.qpp[3 * i + c]);
            k.sym(S_O, zz(int64_t(O[c]) - pr));
        }
    }
}

template <class Sink>
static inline void enc_WH(const FrameCtx &C, int m, size_t i, Sink &k) {
    const int32_t *O = &C.q[3 * i];
    int64_t d[3], hv[3];
    rel(C.q, i + 1, i, d);
    rel(C.q, i + 2, i, hv);
    Circle cc;
    int64_t xs[2], ys[2];
    if (m == 0) {
        SphereCode sc = sphere_encode(d, C.R2);
        put_fe(k, S_FACE, S_E, sc.face, zz(sc.e));
        k.sym(S_KEPT, zz(sc.a));
        k.sym(S_KEPT, zz(sc.b));
        circle_setup(d, C.Rc2R, C.R2, C.invR2, cc);
        int64_t hk = hv[cc.k];
        circle_solve(d, C.R2, cc, hk, xs, ys);
        int sg = std::llabs(hv[cc.i] - xs[1]) + std::llabs(hv[cc.j] - ys[1]) <
                 std::llabs(hv[cc.i] - xs[0]) + std::llabs(hv[cc.j] - ys[0]);
        k.sym(S_V, zz(hk - cc.ck));
        put_rs(k, zz(hv[cc.i] - xs[sg]), zz(hv[cc.j] - ys[sg]), uint32_t(sg));
    } else {
        int64_t p[3], p2[3];
        rel(C.qp, i + 1, i, p);
        rel(C.qp, i + 2, i, p2);
        if (m == 2) {
            int64_t pp[3], pp2[3];
            rel(C.qpp, i + 1, i, pp);
            rel(C.qpp, i + 2, i, pp2);
            for (int c = 0; c < 3; c++) {
                p[c] = clampp(2 * p[c] - pp[c]);
                p2[c] = clampp(2 * p2[c] - pp2[c]);
            }
        }
        int64_t ra, rb, e;
        sph_pred_enc(d, p, C.R2, ra, rb, e);
        k.sym(S_KEPT, zz(ra));
        k.sym(S_KEPT, zz(rb));
        k.sym(S_E, zz(e));
        circle_setup(d, C.Rc2R, C.R2, C.invR2, cc);
        int64_t hk = hv[cc.k];
        circle_solve(d, C.R2, cc, hk, xs, ys);
        int sg = std::llabs(hv[cc.i] - xs[1]) + std::llabs(hv[cc.j] - ys[1]) <
                 std::llabs(hv[cc.i] - xs[0]) + std::llabs(hv[cc.j] - ys[0]);
        int sp = std::llabs(p2[cc.i] - xs[1]) + std::llabs(p2[cc.j] - ys[1]) <
                 std::llabs(p2[cc.i] - xs[0]) + std::llabs(p2[cc.j] - ys[0]);
        k.sym(S_V, zz(hk - p2[cc.k]));
        put_rs(k, zz(hv[cc.i] - xs[sg]), zz(hv[cc.j] - ys[sg]), uint32_t(sg ^ sp));
    }
    if (C.L->nsite == 4) {
        uint32_t z[3];
        for (int c = 0; c < 3; c++) z[c] = zz(int64_t(O[9 + c]) - O[c] - rnd_pred(C.L->vs_a * double(d[c] + hv[c])));
        put_m3(k, z);
    }
}

template <class Sink>
static inline void enc_NB(const FrameCtx &C, int m, size_t i, Sink &k) {
    const uint16_t rf = C.L->ref[i];
    const size_t j = i - (rf >> 4);
    const int64_t R2c = (*C.BR2)[(rf & 15) - 1];
    if (m == 3) {
        for (int c = 0; c < 3; c++) k.sym(S_BK, zz(int64_t(C.q[3 * i + c]) - C.qp[3 * i + c]));
        return;
    }
    int64_t d[3];
    rel(C.q, i, j, d);
    if (m == 0) {
        SphereCode sc = sphere_encode(d, R2c);
        put_fe(k, S_BF, S_BE, sc.face, zz(sc.e));
        k.sym(S_BK, zz(sc.a));
        k.sym(S_BK, zz(sc.b));
    } else {
        int64_t p[3];
        rel(C.qp, i, j, p);
        if (m == 2) {
            int64_t pp[3];
            rel(C.qpp, i, j, pp);
            for (int c = 0; c < 3; c++) p[c] = clampp(2 * p[c] - pp[c]);
        }
        int64_t ra, rb, e;
        sph_pred_enc(d, p, R2c, ra, rb, e);
        k.sym(S_BK, zz(ra));
        k.sym(S_BK, zz(rb));
        k.sym(S_BE, zz(e));
    }
}

static inline int64_t pred_NU(const FrameCtx &C, const int32_t *q, int m, size_t i, int c) {
    switch (m) {
        case 0:
            return i ? q[3 * (i - 1) + c] : 0;
        case 1:
            return C.qp[3 * i + c];
        case 2:
            return clampp(2 * int64_t(C.qp[3 * i + c]) - C.qpp[3 * i + c]);
        default:
            return i ? clampp(int64_t(C.qp[3 * i + c]) + q[3 * (i - 1) + c] - C.qp[3 * (i - 1) + c]) : C.qp[3 * i + c];
    }
}
template <class Sink>
static inline void enc_NU(const FrameCtx &C, int m, size_t i, Sink &k) {
    for (int c = 0; c < 3; c++) k.sym(S_U, zz(int64_t(C.q[3 * i + c]) - pred_NU(C, C.q, m, i, c)));
}

// cached water molecules (their O) that do not fit frame x: they are coded as other atoms in this chunk
template <class T>
static void water_misfits(const T *x, const std::vector<uint32_t> &wo, double r, double rhh, double step,
                          std::vector<uint32_t> &broken) {
    broken.clear();
    auto d2f = [x](size_t a, size_t b) {
        float dx = float(x[3 * a] - x[3 * b]), dy = float(x[3 * a + 1] - x[3 * b + 1]),
              dz = float(x[3 * a + 2] - x[3 * b + 2]);
        return dx * dx + dy * dy + dz * dz;
    };
    const double tol = std::max(0.002, 2.0 * step);
    const float alo = float((r - tol) * (r - tol)), ahi = float((r + tol) * (r + tol));
    const float clo = float((rhh - tol) * (rhh - tol)), chi = float((rhh + tol) * (rhh + tol));
    for (uint32_t o : wo) {
        const float a = d2f(o, o + 1), b = d2f(o, o + 2), c = d2f(o + 1, o + 2);
        if (!(a > alo && a < ahi && b > alo && b < ahi && c > clo && c < chi)) broken.push_back(o);
    }
}

// the water runs (first O, count) of the waters in wo, as the stream stores them
static inline void encode_runs(const std::vector<uint32_t> &wo, int nsite, std::vector<uint8_t> &out) {
    out.resize(16 + 20 * (wo.size() + 1));
    uint8_t *tp = out.data();
    std::vector<uint32_t> runs;
    for (size_t k = 0; k < wo.size();) {
        const size_t k0 = k;
        while (k + 1 < wo.size() && wo[k + 1] == wo[k] + uint32_t(nsite)) k++;
        runs.push_back(wo[k0]);
        runs.push_back(uint32_t(k - k0 + 1));
        k++;
    }
    put_varint(tp, runs.size() / 2);
    uint32_t last = 0;
    for (size_t k = 0; k < runs.size(); k += 2) {
        put_varint(tp, runs[k] - last);
        put_varint(tp, runs[k + 1]);
        last = runs[k] + uint32_t(nsite) * runs[k + 1];
    }
    out.resize(size_t(tp - out.data()));
}

// ------------------------------------------------------------------------------------------------ stream
inline size_t compress_bound(size_t F, size_t N) { return 512 + N * 8 + F * (N * 3 * 9 + 12 * 4096); }

template <class T>
size_t compress_impl(const T *src, size_t F, size_t N, double eb, uint8_t *out, const Options &opt, double forced_step);

// Compresses F frames of N atoms (x, y, z per atom, nm) with absolute bound eb into out, which must hold
// compress_bound(F, N) bytes. Throws std::runtime_error for input it cannot represent (non-finite values, a bound too
// small for the coordinate range, too many atoms); the caller codes those some other way.
template <class T>
size_t compress(const T *src, size_t F, size_t N, double eb, uint8_t *out, const Options &opt = Options()) {
    static_assert(std::is_floating_point<T>::value, "mdc: float or double coordinates");
    if (N > (size_t(1) << 30) || F > (size_t(1) << 30) || !(eb > 0))
        throw std::runtime_error("mdc: unsupported size or bound");
    return compress_impl(src, F, N, eb, out, opt, 0);
}

template <class T>
size_t compress_impl(const T *src, size_t F, size_t N, double eb, uint8_t *out, const Options &opt,
                     double forced_step) {
    uint8_t *p = out;
    // trailing frames that are all one value (the unwritten rest of a chunk) are stored as that value
    size_t nfill = 0;
    T fill = 0;
    if (F > 1 && N > 0) {
        fill = src[(F - 1) * N * 3];
        auto isfill = [src, N, fill](size_t t) {
            const T *X = src + t * N * 3;
            for (size_t i = 0; i < N * 3; i++)
                if (memcmp(&X[i], &fill, sizeof(T)) != 0) return false;
            return true;
        };
        while (nfill + 1 < F && isfill(F - 1 - nfill)) nfill++;
    }
    const size_t Ftot = F;
    F -= nfill;
    // --- lattice step. Candidates: 2 eb, and 2 (eb - ulp) for the binade of max|x| and its neighbours; the first on
    // whose lattice the data already sits is taken (so decompressed data recompresses to the same values), else the
    // candidate of max|x|'s own binade (or eb when the float spacing reaches eb).
    double step, grid_slack = 0.5, proven_step = 0;
    bool step_proven = true;
    {
        const size_t n = F * N * 3;
        const bits_of<T> bmax = max_abs_bits(src, n);
        T amax;
        memcpy(&amax, &bmax, sizeof(T));
        if (!(amax < T(1e30))) throw std::runtime_error("mdc: non-finite input");
        auto step_for = [eb](int e2) {  // e2: binade exponent
            const double ulp = std::ldexp(1.0, e2 - (std::numeric_limits<T>::digits - 1));
            return ulp < 0.5 * eb ? 2.0 * (eb - ulp) : eb;
        };
        int e0;
        std::frexp(double(amax), &e0);
        e0 -= 1;  // amax in [2^e0, 2^(e0+1))
        // step_for(e0) and step_for(e0 + 1) keep every value within eb by construction; 2 eb and step_for(e0 - 1) only
        // when the data sits on their lattice, which the quantization checks
        const double cands[4] = {2.0 * eb, step_for(e0), step_for(e0 - 1), step_for(e0 + 1)};
        const bool proven[4] = {false, true, false, true};
        step = amax > 0 ? cands[1] : 2.0 * eb;
        for (int ci = 0; ci < 4; ci++) {
            const double c = cands[ci];
            if (n == 0 || !(c > 0)) break;
            const double inv0 = 1.0 / c;
            bool on = true;
            for (size_t i = 0, sst = std::max<size_t>(1, n / 256); i < n && on; i += sst) {
                double y = double(src[i]) * inv0;
                on = std::fabs(y - double(rnd(y))) < 0.01;
            }
            if (on) {
                step = c;
                step_proven = proven[ci];
                break;
            }
        }
        if (forced_step > 0) {
            step = forced_step;
            step_proven = true;
        }
        // a lattice that is not proven is checked while quantizing: |x / step - q| must stay below slack, which keeps
        // |q step - x| + ulp(x) / 2 <= eb; one value outside restarts with the proven step
        grid_slack = 0.5 * (1.0 - (double(std::nextafter(amax, T(INFINITY))) - double(amax)) / eb);
        proven_step = cands[1];
        if (!(step > 0) || double(amax) / step > double(1 << 28))
            throw std::runtime_error("mdc: error bound too small for the coordinate range");
    }
    const double inv = 1.0 / step;
    // --- layout. The rigid-water part (which atoms, geometry, virtual site) is kept per thread for the last few
    // (atom count, step) keys and reused while the cached molecules still fit the first frame, which spares one-frame
    // chunks the water detection. The bonds of the other atoms are found afresh for every chunk: their best partners
    // change from frame to frame. The layout is stored in the stream, so decoding does not depend on this; only the
    // compressed bytes can depend on chunks compressed before in the same thread.
    struct WaterEntry {
        size_t N = 0;
        double step = 0;
        bool valid = false;
        Layout L;  // kind, nsite, vs_a, r, rhh
        int64_t R2 = 0;
        double Rc = 0;
        std::vector<uint32_t> wo;   // water O's
        std::vector<uint8_t> runs;  // encoded water runs
        size_t nbroken0 = 0;        // waters that did not fit the frame the layout was found in
    };
    thread_local WaterEntry cache[4];
    thread_local unsigned next_slot = 0;
    WaterEntry *E = nullptr;
    thread_local std::vector<uint32_t> broken;  // cached waters that do not fit this chunk
    broken.clear();
    for (auto &c : cache)
        if (c.valid && c.N == N && c.step == step) {
            E = &c;
            break;
        }
    if (E) {
        water_misfits(src, E->wo, E->L.r, E->L.rhh, step, broken);
        // molecules split over the periodic boundary break a steady share of the waters in atom-wrapped data: detect
        // again only when clearly more are broken than when the layout was found
        if (broken.size() > 2 * E->nbroken0 + E->wo.size() / 500) E = nullptr;
    }
    if (!E) {
        E = &cache[next_slot++ % 4];
        E->valid = false;
        E->N = N;
        E->step = step;
        Layout &L0 = E->L;
        detect_water(src, N, L0, step);
        E->R2 = 0;
        E->Rc = 0;
        constexpr double RMAX = 16384;
        if (L0.r > 0 && L0.r / step < RMAX) {
            double Ru = L0.r / step, cth = (2 * L0.r * L0.r - L0.rhh * L0.rhh) / (2 * L0.r * L0.r);
            E->R2 = rnd(Ru * Ru);
            E->Rc = Ru * cth;
        } else if (L0.r > 0) {
            std::fill(L0.kind.begin(), L0.kind.end(), 2);
            L0.r = 0;
            L0.nsite = 3;
        }
        // a gap of one or two molecules between waters is water split over the boundary in this frame
        if (L0.r > 0) {
            const size_t ns = size_t(L0.nsite);
            size_t prev_end = 0;  // one past the last water seen, 0: none yet
            for (size_t i = 0; i < N;) {
                if (L0.kind[i] != 0) {
                    i++;
                    continue;
                }
                const size_t gap = i - prev_end;
                if (prev_end && gap > 0 && gap % ns == 0 && gap <= 2 * ns) {
                    bool free_ = true;
                    for (size_t a = prev_end; a < i; a++) free_ &= L0.kind[a] == 2;
                    if (free_)
                        for (size_t o = prev_end; o < i; o += ns) {
                            L0.kind[o] = 0;
                            for (size_t a = 1; a < ns; a++) L0.kind[o + a] = 1;
                        }
                }
                i += ns;
                prev_end = i;
            }
        }
        E->wo.clear();
        for (size_t i = 0; i < N; i++)
            if (L0.kind[i] == 0) E->wo.push_back(uint32_t(i));
        encode_runs(E->wo, L0.nsite, E->runs);
        water_misfits(src, E->wo, E->L.r, E->L.rhh, step, broken);
        E->nbroken0 = broken.size();
        E->valid = true;
    }
    // this chunk's water: the cached one without the molecules that do not fit
    thread_local Layout LB;
    LB.kind = E->L.kind;
    LB.nsite = E->L.nsite;
    LB.vs_a = E->L.vs_a;
    LB.r = E->L.r;
    LB.rhh = E->L.rhh;
    const std::vector<uint32_t> *wo = &E->wo;
    const std::vector<uint8_t> *wruns = &E->runs;
    thread_local std::vector<uint32_t> wo_here;
    thread_local std::vector<uint8_t> runs_here;
    if (!broken.empty()) {
        for (uint32_t o : broken)
            for (int a = 0; a < LB.nsite; a++) LB.kind[o + a] = 2;
        wo_here.clear();
        for (uint32_t o : E->wo)
            if (LB.kind[o] == 0) wo_here.push_back(o);
        encode_runs(wo_here, LB.nsite, runs_here);
        wo = &wo_here;
        wruns = &runs_here;
    }
    // bonds of the other atoms, for this chunk
    detect_bonds(src, N, LB);
    thread_local std::vector<int64_t> BR2;
    BR2.assign(LB.blen.size(), 0);
    for (size_t c = 0; c < LB.blen.size(); c++) BR2[c] = rnd((LB.blen[c] / step) * (LB.blen[c] / step));
    for (size_t i = 0; i < N; i++)
        if (LB.ref[i] && BR2[(LB.ref[i] & 15) - 1] >= RMAXB2) LB.ref[i] = 0;
    thread_local std::vector<uint32_t> ubuf[2];  // bonded, unbonded
    for (auto &u : ubuf) u.clear();
    for (size_t i = 0; i < N; i++)
        if (LB.kind[i] == 2) ubuf[LB.ref[i] ? 0 : 1].push_back(uint32_t(i));
    const Layout &L = LB;
    const int64_t R2 = E->R2;
    const double Rc = E->Rc;
    const double Rc2R = R2 > 0 ? Rc / (2.0 * std::sqrt(double(R2))) : 0.0, invR2 = R2 > 0 ? 1.0 / double(R2) : 0.0;

    // --- header
    put_raw(p, MAGIC);
    put_raw(p, uint32_t(N));
    put_raw(p, uint32_t(Ftot));
    put_raw(p, uint32_t(nfill));
    put_raw(p, fill);
    put_raw(p, step);
    put_raw(p, R2);
    put_raw(p, Rc);
    *p++ = uint8_t(L.nsite);
    put_raw(p, L.vs_a);
    *p++ = uint8_t(L.blen.size());
    for (auto b : BR2) put_raw(p, b);
    memcpy(p, wruns->data(), wruns->size());
    p += wruns->size();
    {  // references of the other atoms, coded in the context of the previous one's: the previous values seen at
       // least 16 times (up to MAXCTX, most frequent first) have their own tables, the rest share one
        constexpr size_t REFV = 16 * (MAXOFF + 1);
        uint32_t pc[REFV] = {0};
        uint16_t prev = 0;
        for (size_t i = 0; i < N; i++)
            if (L.kind[i] == 2) {
                pc[prev]++;
                prev = L.ref[i];
            }
        uint16_t cval[MAXCTX];
        size_t nc = 0;
        {
            // most frequent first, lowest value on ties (a sort with a by-reference comparator lambda here was
            // miscompiled by GCC 13 -O3, ipa-modref)
            uint32_t left[REFV];
            memcpy(left, pc, sizeof(pc));
            while (nc < MAXCTX) {
                uint32_t b = 0;
                for (uint32_t v = 1; v < REFV; v++)
                    if (left[v] > left[b]) b = v;
                if (left[b] < 16) break;
                cval[nc++] = uint16_t(b);
                left[b] = 0;
            }
        }
        uint8_t cid[REFV] = {0};
        for (size_t c = 0; c < nc; c++) cid[cval[c]] = uint8_t(c + 1);
        thread_local std::vector<Huff> hc;
        if (hc.size() < MAXCTX + 1) hc.resize(MAXCTX + 1);
        for (size_t c = 0; c <= nc; c++) hc[c].reset();
        prev = 0;
        for (size_t i = 0; i < N; i++)
            if (L.kind[i] == 2) {
                hc[cid[prev]].count(L.ref[i]);
                prev = L.ref[i];
            }
        *p++ = uint8_t(nc);
        for (size_t c = 0; c < nc; c++) *p++ = uint8_t(cval[c]);
        BitWriter bw(p);
        for (size_t c = 0; c <= nc; c++) {
            hc[c].build();
            hc[c].write_table(bw);
        }
        prev = 0;
        for (size_t i = 0; i < N; i++)
            if (L.kind[i] == 2) {
                hc[cid[prev]].put(bw, L.ref[i]);
                prev = L.ref[i];
            }
        p = bw.finish();
    }
    const std::vector<uint32_t> *units[NGROUP] = {wo, wo, &ubuf[0], &ubuf[1]};

    struct Work {
        std::vector<int32_t> b[3];
        std::vector<Tok> seq;  // sized for the largest frame seen
        Huff H[NS];            // this frame's tables
        Huff hprev[NS];        // the previous frame's
        WaterBatch wb;
        SphereBatch sb;
        std::vector<uint32_t> hist[NS], touched[NS];
    };
    thread_local Work W;
    for (int k = 0; k < 3; k++)
        if (W.b[k].size() < N * 3) W.b[k].resize(N * 3);
    for (int s = 0; s < NS; s++)
        if (W.hist[s].size() < ALPHA) W.hist[s].assign(ALPHA, 0);
    if (W.seq.size() < N * 12 + 16) W.seq.resize(N * 12 + 16);
    int32_t *q = W.b[0].data(), *qp = W.b[1].data(), *qpp = W.b[2].data();
    EstSink est{W.hist, W.touched};

    const Kernels K = simd_kernels(opt.simd);
    bool nb_small = true;  // the batch needs exact products: squared bond lengths below 2^22
    for (int64_t b : BR2) nb_small &= b < (int64_t(1) << 22);
    int mode[NGROUP] = {0, 0, 0, 0};
    for (size_t t = 0; t < F; t++) {
        const T *X = src + t * N * 3;
        if (step_proven)
            quantize(X, N * 3, inv, q);
        else if (quantize_checked(X, N * 3, inv, grid_slack, q))
            return compress_impl(src, Ftot, N, eb, out, opt, proven_step);
        FrameCtx C{q, qp, qpp, &L, R2, Rc2R, invR2, &BR2, {0, 0, 0}, 1, 1, 1, 0, {0, 0, 0}, true};
        if (!units[G_O]->empty()) {
            int32_t mn[3] = {INT32_MAX, INT32_MAX, INT32_MAX}, mx[3] = {INT32_MIN, INT32_MIN, INT32_MIN};
            for (uint32_t i : *units[G_O])
                for (int c = 0; c < 3; c++) {
                    mn[c] = std::min(mn[c], q[3 * i + c]);
                    mx[c] = std::max(mx[c], q[3 * i + c]);
                }
            for (int c = 0; c < 3; c++) C.omin[c] = mn[c];
            C.Rx = uint64_t(int64_t(mx[0]) - mn[0] + 1);
            C.Ry = uint64_t(int64_t(mx[1]) - mn[1] + 1);
            C.Rz = uint64_t(int64_t(mx[2]) - mn[2] + 1);
            box_bits(C.Rx, C.Ry, C.Rz, C.obits, C.ob, C.mixed);
        }
        // --- per-group predictor choice on a sample: at frames 1 and 2 (the second-order modes need two past
        // frames), then every 8th frame; frame 0 is intra
        if (t == 0)
            for (int g = 0; g < NGROUP; g++) mode[g] = 0;
        if (t > 0 && (t <= 2 || t % 8 == 0)) {
            for (int g = 0; g < NGROUP; g++) {
                const auto &u = *units[g];
                if (u.empty()) continue;
                const size_t stride = std::max<size_t>(1, u.size() / 256), ns = (u.size() + stride - 1) / stride;
                double best = 1e300;
                for (int m = 0; m < NMODE[g]; m++) {
                    if (m == 2 && t < 2) continue;
                    double c;
                    const double scale = double(u.size()) / double(ns);
                    if (g == G_O && m == 0) {
                        c = double(C.mixed ? C.obits : C.ob[0] + C.ob[1] + C.ob[2]) * double(u.size());
                    } else {
                        for (size_t k = 0; k < u.size(); k += stride) {
                            size_t i = u[k];
                            if (g == G_O)
                                enc_O(C, m, i, est);
                            else if (g == G_WH)
                                enc_WH(C, m, i, est);
                            else if (g == G_NB)
                                enc_NB(C, m, i, est);
                            else
                                enc_NU(C, m, i, est);
                        }
                        c = est.take(scale);
                    }
                    if (c < best) {
                        best = c;
                        mode[g] = m;
                    }
                }
            }
        }
        // --- intra water geometry of the whole frame in one SIMD batch (bit-identical to enc_WH, mode 0)
        const bool batch = mode[G_WH] == 0 && K.water && R2 > 0 && R2 < (int64_t(1) << 22) && !units[G_O]->empty();
        if (batch) {
            const auto &u = *units[G_O];
            WaterBatch &B = W.wb;
            B.resize(u.size());
            for (size_t w = 0; w < u.size(); w++) {
                const int32_t *O = &q[3 * u[w]];
                B.dx[w] = double(O[3] - O[0]);
                B.dy[w] = double(O[4] - O[1]);
                B.dz[w] = double(O[5] - O[2]);
                B.hx[w] = double(O[6] - O[0]);
                B.hy[w] = double(O[7] - O[1]);
                B.hz[w] = double(O[8] - O[2]);
                if (L.nsite == 4) {
                    B.mx[w] = double(O[9] - O[0]);
                    B.my[w] = double(O[10] - O[1]);
                    B.mz[w] = double(O[11] - O[2]);
                }
            }
            for (size_t w = u.size(); w < u.size() + SIMD_PAD; w++)
                B.dx[w] = B.dy[w] = B.dz[w] = B.hx[w] = B.hy[w] = B.hz[w] = B.mx[w] = B.my[w] = B.mz[w] = 0;
            K.water(B, u.size(), R2, Rc2R, invR2, L.nsite == 4, L.vs_a);
        }
        // --- intra sphere code of the bonded atoms in one batch (bit-identical to enc_NB, mode 0)
        const bool nbatch = mode[G_NB] == 0 && K.sphere && nb_small && !units[G_NB]->empty();
        if (nbatch) {
            const auto &u = *units[G_NB];
            SphereBatch &B = W.sb;
            B.resize(u.size());
            for (size_t w = 0; w < u.size(); w++) {
                const size_t i = u[w], j = i - (L.ref[i] >> 4);
                B.dx[w] = double(q[3 * i] - q[3 * j]);
                B.dy[w] = double(q[3 * i + 1] - q[3 * j + 1]);
                B.dz[w] = double(q[3 * i + 2] - q[3 * j + 2]);
                B.r2[w] = double(BR2[(L.ref[i] & 15) - 1]);
            }
            for (size_t w = u.size(); w < u.size() + SIMD_PAD; w++) B.dx[w] = B.dy[w] = B.dz[w] = B.r2[w] = 0;
            K.sphere(B, u.size());
        }
        // --- symbols in atom order
        Huff *H = W.H;
        for (int s = 0; s < NS; s++) H[s].reset();
        StoreSink sink{W.seq.data(), H};
        size_t iw = 0, ib = 0;
        for (size_t i = 0; i < N; i++) {
            const int k = L.kind[i];
            if (k == 0 && batch) {
                enc_O(C, mode[G_O], i, sink);
                const WaterBatch &B = W.wb;
                put_fe(sink, S_FACE, S_E, uint32_t(B.face[iw]), zz(B.e[iw]));
                sink.sym(S_KEPT, zz(B.a[iw]));
                sink.sym(S_KEPT, zz(B.b[iw]));
                sink.sym(S_V, zz(B.v[iw]));
                put_rs(sink, zz(B.r1[iw]), zz(B.r2[iw]), uint32_t(B.side[iw]));
                if (L.nsite == 4) {
                    const uint32_t z[3] = {zz(B.m0[iw]), zz(B.m1[iw]), zz(B.m2[iw])};
                    put_m3(sink, z);
                }
                iw++;
            } else if (k == 0) {
                enc_O(C, mode[G_O], i, sink);
                enc_WH(C, mode[G_WH], i, sink);
            } else if (k == 2) {
                if (!L.ref[i]) {
                    enc_NU(C, mode[G_NU], i, sink);
                } else if (nbatch) {
                    const SphereBatch &B = W.sb;
                    put_fe(sink, S_BF, S_BE, uint32_t(B.face[ib]), zz(B.e[ib]));
                    sink.sym(S_BK, zz(B.a[ib]));
                    sink.sym(S_BK, zz(B.b[ib]));
                    ib++;
                } else {
                    enc_NB(C, mode[G_NB], i, sink);
                }
            }
        }
        for (int s = 0; s < NS; s++) H[s].build();
        // frames after the first may keep a stream's table from the previous frame when that codes this frame no worse
        bool keep[NS] = {false};
        if (t > 0)
            for (int s = 0; s < NS; s++) {
                if (!W.hprev[s].used || !H[s].used) continue;
                const uint64_t old_bits = H[s].coded_bits_with(W.hprev[s]);
                if (old_bits == ~uint64_t(0)) continue;
                const uint64_t new_bits = H[s].coded_bits_with(H[s]) + H[s].table_bits();
                if (old_bits <= new_bits) {
                    keep[s] = true;
                    H[s].use_code_of(W.hprev[s]);
                }
            }
        BitWriter bw(p);
        for (int g = 0; g < NGROUP; g++) bw.put(uint32_t(mode[g]), 2);
        if (!units[G_O]->empty() && mode[G_O] == 0) {
            for (int c = 0; c < 3; c++) bw.put(uint32_t(C.omin[c]), 32);
            bw.put(uint32_t(C.Rx - 1), 32);
            bw.put(uint32_t(C.Ry - 1), 32);
            bw.put(uint32_t(C.Rz - 1), 32);
        }
        for (int s = 0; s < NS; s++) {
            if (t > 0 && W.hprev[s].used) bw.put(keep[s], 1);
            if (!keep[s]) H[s].write_table(bw);
        }
        for (int s = 0; s < NS; s++) W.hprev[s].use_code_of(H[s]);
        for (const Tok *tk = W.seq.data(); tk != sink.t; tk++) {
            if (tk->s == 255) {
                bw.put64(uint64_t(tk->v) | (uint64_t(tk->k) << 32), tk->nb);
                continue;
            }
            const uint32_t e = H[tk->s].enc[tk->k];
            const int len = int(e >> 16);
            if (tk->k < DIRECT) {
                bw.put(e & 0xffff, len);
            } else {  // the code and the raw low bits in one write
                const int nb = rawbits(tk->k);
                bw.put64(uint64_t(e & 0xffff) | (uint64_t(tk->v & ((1u << nb) - 1)) << len), len + nb);
            }
        }
        p = bw.finish();
        int32_t *tmp = qpp;
        qpp = qp;
        qp = q;
        q = tmp;
    }
    return size_t(p - out);
}

// ------------------------------------------------------------------------------------------------ decoder
struct Src {
    BitReader *br;
    const Huff *H;
    inline uint32_t sym(int s) { return H[s].get(*br); }
    inline uint64_t raw(int nb) { return br->get64(nb); }
    inline void fe(int s, int se, uint32_t &face, uint32_t &ze) {
        uint32_t v = sym(s);
        face = v >> 4;
        ze = v & 15;
        if (ze == 15) ze = sym(se);
    }
    inline void rs(uint32_t &r1, uint32_t &r2, uint32_t &side) {
        uint32_t v = sym(S_RS);
        if (v < 98) {
            side = v & 1;
            v >>= 1;
            r1 = v / 7;
            r2 = v % 7;
        } else {
            side = (v - 98) & 1;
            r1 = sym(S_R);
            r2 = sym(S_R);
        }
    }
    inline void m3(uint32_t z[3]) {
        uint32_t v = sym(S_MJ);
        if (v < 125) {
            z[2] = v % 5;
            v /= 5;
            z[1] = v % 5;
            z[0] = v / 5;
        } else {
            for (int c = 0; c < 3; c++) z[c] = sym(S_M);
        }
    }
};

// The number of values (frames * atoms * 3) a stream holds, 0 if it is not an ALGO_MDC stream.
inline size_t stored_values(const uint8_t *in, size_t size) {
    if (size < 12) return 0;
    uint32_t magic, N, F;
    memcpy(&magic, in, 4);
    memcpy(&N, in + 4, 4);
    memcpy(&F, in + 8, 4);
    return magic == MAGIC ? size_t(N) * F * 3 : 0;
}

// Decompresses into dst, which holds stored_values(in, size) values. Throws std::runtime_error on a corrupt stream.
template <class T>
void decompress(const uint8_t *in, size_t size, T *dst) {
    if (stored_values(in, size) == 0) throw std::runtime_error("mdc: not an ALGO_MDC stream");
    const uint8_t *p = in, *end = in + size;
    uint32_t magic, N32, F32, nfill32;
    T fill;
    double step, Rc, vs_a;
    int64_t R2;
    get_raw(p, end, magic);
    get_raw(p, end, N32);
    get_raw(p, end, F32);
    get_raw(p, end, nfill32);
    get_raw(p, end, fill);
    get_raw(p, end, step);
    get_raw(p, end, R2);
    get_raw(p, end, Rc);
    if (nfill32 >= F32 && F32 > 0) throw std::runtime_error("mdc: corrupt fill count");
    for (size_t i = size_t(F32 - nfill32) * N32 * 3; i < size_t(F32) * N32 * 3; i++) dst[i] = fill;
    const size_t N = N32, F = F32 - nfill32;
    Layout L;
    uint8_t nsite, ncls;
    get_raw(p, end, nsite);
    get_raw(p, end, vs_a);
    get_raw(p, end, ncls);
    L.nsite = nsite;
    L.vs_a = vs_a;
    if (L.nsite != 3 && L.nsite != 4) throw std::runtime_error("mdc: bad water model");
    if (R2 < 0 || R2 >= (int64_t(1) << 30)) throw std::runtime_error("mdc: bad water geometry");
    std::vector<int64_t> BR2(ncls);
    for (auto &b : BR2) {
        get_raw(p, end, b);  // classes of RMAXB2 or more are stored but no atom refers to them
        if (b < 0) throw std::runtime_error("mdc: bad bond length");
    }
    const double Rc2R = R2 > 0 ? Rc / (2.0 * std::sqrt(double(R2))) : 0.0, invR2 = R2 > 0 ? 1.0 / double(R2) : 0.0;
    L.kind.assign(N, 2);
    {
        const uint64_t nr = get_varint(p, end);
        uint64_t last = 0;
        for (uint64_t k = 0; k < nr; k++) {
            const uint64_t s = get_varint(p, end) + last, c = get_varint(p, end);
            if (s > N || c > N || s + uint64_t(L.nsite) * c > N) throw std::runtime_error("mdc: corrupt water runs");
            for (uint64_t m = 0; m < c; m++) {
                L.kind[s + L.nsite * m] = 0;
                for (int a = 1; a < L.nsite; a++) L.kind[s + L.nsite * m + a] = 1;
            }
            last = s + L.nsite * c;
        }
    }
    L.ref.assign(N, 0);
    size_t nwat = 0;
    {
        uint8_t nc8;
        get_raw(p, end, nc8);
        const size_t nc = nc8;
        if (nc > MAXCTX || size_t(end - p) < nc) throw std::runtime_error("mdc: corrupt references");
        uint8_t cid[256] = {0};  // context of each previous value (they are below 16 * (MAXOFF + 1) <= 256)
        for (size_t c = 0; c < nc; c++) cid[*p++] = uint8_t(c + 1);
        BitReader br(p, end);
        thread_local std::vector<Huff> hc;
        if (hc.size() < nc + 1) hc.resize(nc + 1);
        for (size_t c = 0; c <= nc; c++) hc[c].read_table(br, false);
        // one flat table: each context's entries over one period of its longest code, every entry holding the value,
        // its code length and where the next value's context starts (a single load per value)
        struct Ctx {
            uint32_t base, mask;
        };
        Ctx cx[MAXCTX + 1];
        thread_local std::vector<uint64_t> flat;
        size_t tot = 0;
        for (size_t c = 0; c <= nc; c++) {
            uint32_t lmax = 0;
            if (hc[c].used && hc[c].single < 0)
                for (uint32_t k = 0; k <= hc[c].maxsym; k++) lmax = std::max(lmax, uint32_t(hc[c].len[k]));
            cx[c] = {uint32_t(tot), (1u << lmax) - 1};
            tot += size_t(1) << lmax;
        }
        flat.resize(tot);
        for (size_t c = 0; c <= nc; c++)
            for (uint32_t k = 0; k <= cx[c].mask; k++) {
                const uint32_t e = hc[c].dt[k], v = e >> 4;
                if (v >= 256) throw std::runtime_error("mdc: corrupt references");
                const Ctx &nx = cx[cid[v]];
                flat[cx[c].base + k] =
                    uint64_t(v) | uint64_t(e & 15) << 16 | uint64_t(nx.mask) << 20 | uint64_t(nx.base) << 32;
            }
        uint32_t base = cx[cid[0]].base, mask = cx[cid[0]].mask;
        for (size_t i = 0; i < N; i++) {
            nwat += L.kind[i] == 0;
            if (L.kind[i] == 2) {
                const uint64_t e = flat[base + (br.peek(HMAXLEN) & mask)];
                br.skip(int((e >> 16) & 15));
                const uint16_t r = uint16_t(e & 0xffff);
                base = uint32_t(e >> 32);
                mask = uint32_t(e >> 20) & 0xfff;
                if (r && (size_t(r >> 4) > i || (r & 15) == 0 || size_t(r & 15) > ncls))
                    throw std::runtime_error("mdc: corrupt references");
                L.ref[i] = r;
            }
        }
        p = br.align();
    }
    std::vector<int32_t> B[3];
    for (auto &b : B) b.assign(N * 3, 0);
    std::vector<Huff> HD(NS);  // tables persist across the frames of the chunk
    int32_t *q = B[0].data(), *qp = B[1].data(), *qpp = B[2].data();
    for (size_t t = 0; t < F; t++) {
        BitReader br(p, end);
        int mode[NGROUP];
        for (int g = 0; g < NGROUP; g++) mode[g] = int(br.get(2));
        FrameCtx C{q, qp, qpp, &L, R2, Rc2R, invR2, &BR2, {0, 0, 0}, 1, 1, 1, 0, {0, 0, 0}, true};
        if (nwat && mode[G_O] == 0) {
            for (int c = 0; c < 3; c++) C.omin[c] = int32_t(br.get(32));
            C.Rx = uint64_t(br.get(32)) + 1;
            C.Ry = uint64_t(br.get(32)) + 1;
            C.Rz = uint64_t(br.get(32)) + 1;
            box_bits(C.Rx, C.Ry, C.Rz, C.obits, C.ob, C.mixed);
        }
        Huff *H = HD.data();
        for (int s = 0; s < NS; s++) {
            const bool kept = t > 0 && H[s].used && br.get(1);
            if (!kept) H[s].read_table(br);
        }
        Src S{&br, H};
        for (size_t i = 0; i < N; i++) {
            const int k = L.kind[i];
            if (k == 0) {
                int32_t *O = &q[3 * i];
                const int mo = mode[G_O];
                if (mo == 0) {
                    if (C.mixed) {
                        uint64_t idx = S.raw(C.obits);
                        O[0] = int32_t(int64_t(idx % C.Rx) + C.omin[0]);
                        idx /= C.Rx;
                        O[1] = int32_t(int64_t(idx % C.Ry) + C.omin[1]);
                        idx /= C.Ry;
                        O[2] = int32_t(int64_t(idx) + C.omin[2]);
                    } else {
                        for (int c = 0; c < 3; c++) O[c] = int32_t(int64_t(S.raw(C.ob[c])) + C.omin[c]);
                    }
                } else {
                    for (int c = 0; c < 3; c++) {
                        int64_t pr = mo == 1 ? qp[3 * i + c] : clampp(2 * int64_t(qp[3 * i + c]) - qpp[3 * i + c]);
                        O[c] = int32_t(pr + unzz(S.sym(S_O)));
                    }
                }
                const int mw = mode[G_WH];
                int64_t d[3], hv[3];
                Circle cc;
                int64_t xs[2], ys[2];
                if (mw == 0) {
                    uint32_t face, ze, r1, r2, side;
                    S.fe(S_FACE, S_E, face, ze);
                    int64_t a = unzz(S.sym(S_KEPT)), b = unzz(S.sym(S_KEPT));
                    sphere_decode(face, a, b, unzz(ze), R2, d);
                    circle_setup(d, Rc2R, R2, invR2, cc);
                    int64_t hk = cc.ck + unzz(S.sym(S_V));
                    S.rs(r1, r2, side);
                    const int sg = int(side);
                    circle_solve(d, R2, cc, hk, xs, ys);
                    hv[cc.k] = hk;
                    hv[cc.i] = xs[sg] + unzz(r1);
                    hv[cc.j] = ys[sg] + unzz(r2);
                } else {
                    int64_t pv[3], p2[3];
                    rel(qp, i + 1, i, pv);
                    rel(qp, i + 2, i, p2);
                    if (mw == 2) {
                        int64_t pp[3], pp2[3];
                        rel(qpp, i + 1, i, pp);
                        rel(qpp, i + 2, i, pp2);
                        for (int c = 0; c < 3; c++) {
                            pv[c] = clampp(2 * pv[c] - pp[c]);
                            p2[c] = clampp(2 * p2[c] - pp2[c]);
                        }
                    }
                    int64_t ra = unzz(S.sym(S_KEPT)), rb = unzz(S.sym(S_KEPT)), e = unzz(S.sym(S_E));
                    sph_pred_dec(pv, R2, ra, rb, e, d);
                    circle_setup(d, Rc2R, R2, invR2, cc);
                    int64_t hk = p2[cc.k] + unzz(S.sym(S_V));
                    uint32_t r1, r2, side;
                    S.rs(r1, r2, side);
                    circle_solve(d, R2, cc, hk, xs, ys);
                    int sp = std::llabs(p2[cc.i] - xs[1]) + std::llabs(p2[cc.j] - ys[1]) <
                             std::llabs(p2[cc.i] - xs[0]) + std::llabs(p2[cc.j] - ys[0]);
                    int sg = int(side) ^ sp;
                    hv[cc.k] = hk;
                    hv[cc.i] = xs[sg] + unzz(r1);
                    hv[cc.j] = ys[sg] + unzz(r2);
                }
                for (int c = 0; c < 3; c++) {
                    O[3 + c] = int32_t(O[c] + d[c]);
                    O[6 + c] = int32_t(O[c] + hv[c]);
                }
                if (L.nsite == 4) {
                    uint32_t z[3];
                    S.m3(z);
                    for (int c = 0; c < 3; c++)
                        O[9 + c] = int32_t(O[c] + rnd_pred(vs_a * double(d[c] + hv[c])) + unzz(z[c]));
                }
            } else if (k == 2) {
                int32_t *A = &q[3 * i];
                if (L.ref[i]) {
                    const size_t j = i - (L.ref[i] >> 4);
                    const int64_t R2c = BR2[(L.ref[i] & 15) - 1];
                    int64_t d[3];
                    if (mode[G_NB] == 3) {
                        for (int c = 0; c < 3; c++) A[c] = int32_t(qp[3 * i + c] + unzz(S.sym(S_BK)));
                        continue;
                    }
                    if (mode[G_NB] == 0) {
                        uint32_t face, ze;
                        S.fe(S_BF, S_BE, face, ze);
                        int64_t a = unzz(S.sym(S_BK)), b = unzz(S.sym(S_BK));
                        sphere_decode(face, a, b, unzz(ze), R2c, d);
                    } else {
                        int64_t pv[3];
                        rel(qp, i, j, pv);
                        if (mode[G_NB] == 2) {
                            int64_t pp[3];
                            rel(qpp, i, j, pp);
                            for (int c = 0; c < 3; c++) pv[c] = clampp(2 * pv[c] - pp[c]);
                        }
                        int64_t ra = unzz(S.sym(S_BK)), rb = unzz(S.sym(S_BK)), e = unzz(S.sym(S_BE));
                        sph_pred_dec(pv, R2c, ra, rb, e, d);
                    }
                    for (int c = 0; c < 3; c++) A[c] = int32_t(q[3 * j + c] + d[c]);
                } else {
                    for (int c = 0; c < 3; c++) A[c] = int32_t(pred_NU(C, q, mode[G_NU], i, c) + unzz(S.sym(S_U)));
                }
            }
        }
        T *D = dst + t * N * 3;
        for (size_t i = 0; i < N * 3; i++) D[i] = T(double(q[i]) * step);
        p = br.align();
        int32_t *tmp = qpp;
        qpp = qp;
        qp = q;
        q = tmp;
    }
}

}  // namespace mdc
}  // namespace SZ3
#endif
