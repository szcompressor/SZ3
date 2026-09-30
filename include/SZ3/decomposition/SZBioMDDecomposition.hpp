#ifndef SZ3_BIOMD_DECOMPOSITION_HPP
#define SZ3_BIOMD_DECOMPOSITION_HPP

// ALGO_BIOMD: molecular-dynamics coordinates {frames, atoms, 3} (nm, absolute bound) as integer symbols for the
// Huffman encoder. Every coordinate goes on the lattice q = round(x / step), |x - q step| <= eb, and all prediction is
// integer arithmetic on that lattice, which the decoder repeats exactly.
//  * rigid water (O, H, H): each H on the sphere |H - O| = r (the virtual site of 4-site models is a bonded atom)
//  * other atoms: on the sphere of a bond-length class around one of the previous MAXOFF atoms, else a delta
//  * per frame and atom group, intra or from the previous frame, whichever is cheaper on a sample (frame 0 is intra)
// The layout (which atoms are water, the bond of each other atom) is found on a chunk's first frame and stored with it.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "SZ3/decomposition/Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {
namespace biomd {

constexpr int MAXOFF = 4;                     // a bond partner is one of the previous MAXOFF atoms
constexpr int MAXCLS = 15;                    // bond-length classes
constexpr int64_t RMAXB2 = int64_t(1) << 28;  // bonds of 16384 lattice units or more are not coded as spheres
enum { G_O, G_WH, G_NB, G_NU, NGROUP };
enum { S_W, S_REF, S_O, S_FE, S_KEPT, S_E, S_CA, S_CE, S_BFE, S_BK, S_BE, S_U, NS };
// each stream is split by the mode (0 intra, 1 previous frame) of its group, so each (stream, mode) gets its own code
constexpr int SGROUP[NS] = {-1, -1, G_O, G_WH, G_WH, G_WH, G_WH, G_WH, G_NB, G_NB, G_NB, G_NU};
constexpr int PER_UNIT[NS] = {1, 1, 3, 2, 4, 2, 1, 2, 1, 3, 1, 3};  // most symbols per unit of the group, per frame

inline uint32_t zz(int64_t v) { return uint32_t((uint64_t(v) << 1) ^ uint64_t(v >> 63)); }
inline int64_t unzz(uint32_t u) { return int64_t(u >> 1) ^ -int64_t(u & 1); }
inline int64_t rnd(double y) { return int64_t(y + std::copysign(0.5, y)); }
inline int64_t isqrt_round(int64_t n) { return n <= 0 ? 0 : int64_t(std::sqrt(double(n)) + 0.5); }
inline void rel(const int32_t *a, size_t i, size_t j, int64_t o[3]) {
    for (int c = 0; c < 3; c++) o[c] = int64_t(a[3 * i + c]) - a[3 * j + c];
}

// The largest axis f of p, and the other two.
ALWAYS_INLINE void axes(const int64_t p[3], int &f, int &i, int &j) {
    f = std::llabs(p[1]) > std::llabs(p[0]) ? 1 : 0;
    if (std::llabs(p[2]) > std::llabs(p[f])) f = 2;
    i = f == 2 ? 0 : f + 1;
    j = f == 0 ? 2 : (f == 1 ? 0 : 1);
}

// ------------------------------------------------------------------------------------------------ layout
struct Layout {
    std::vector<uint8_t> kind;  // per atom: 0 = water O (then its H at i + 1, i + 2), 1 = water H, 2 = other
    std::vector<uint16_t>
        ref;               // other atoms: 0 = no bond, else off * 16 + cls + 1 (the sphere of class cls around i - off)
    double r = 0, hh = 0;  // water O-H and H-H (nm)
    std::vector<double> blen;  // bond-length classes (nm)
};

// Rigid water on frame x: kind and r. The tolerances grow with the lattice step, so rounded input (xtc files,
// or data this codec decompressed) still fits.
template <class T>
void detect_water(const T *x, size_t N, Layout &L, double step) {
    L.kind.assign(N, 2);
    L.r = 0;
    auto d2 = [x](size_t a, size_t b) {
        float dx = float(x[3 * a] - x[3 * b]), dy = float(x[3 * a + 1] - x[3 * b + 1]),
              dz = float(x[3 * a + 2] - x[3 * b + 2]);
        return dx * dx + dy * dy + dz * dz;
    };
    // Candidates: O-H in [0.08, 0.125] nm, H-H in [0.13, 0.2] nm, probed across the frame. Rigid water gives sharp
    // O-H / H-H peaks, flexible CH2/NH2 groups broad ones.
    std::vector<float> cand_oh, cand_hh;
    const size_t probe = std::max<size_t>(1, N / 300);
    for (size_t i = 0; i + 2 < N; i += probe)
        for (size_t o = 0; o < 3 && i + o + 2 < N; o++) {
            const size_t k = i + o;
            const float a = d2(k, k + 1), b = d2(k, k + 2), c = d2(k + 1, k + 2);
            if (a > 0.0064f && a < 0.015625f && b > 0.0064f && b < 0.015625f && c > 0.0169f && c < 0.04f) {
                cand_hh.push_back(std::sqrt(c));
                cand_oh.push_back(std::sqrt(a));
                cand_oh.push_back(std::sqrt(b));
                break;
            }
        }
    if (cand_hh.size() < 16) return;
    // tolerances stay well below the O-H / C-H difference (> 0.01 nm) that separates water from CH2/NH2
    const double tight = std::max(0.001, step), tol = std::max(0.002, 2.0 * step), loose = std::max(0.01, 3.0 * tight);
    auto mode = [](std::vector<float> v) {  // median of the values within 0.003 nm of the heaviest 0.004 nm window
        std::sort(v.begin(), v.end());
        size_t best = 0, n = 0;
        for (size_t lo = 0, hi = 0; hi < v.size(); hi++) {
            while (v[hi] - v[lo] > 0.004f) lo++;
            if (hi - lo + 1 > n) n = hi - lo + 1, best = (lo + hi) / 2;
        }
        std::vector<float> w;
        for (float y : v)
            if (std::fabs(y - v[best]) < 0.003f) w.push_back(y);
        return double(w[w.size() / 2]);
    };
    const double h0 = mode(cand_hh);
    std::vector<float> ohs;
    for (size_t c = 0; c < cand_hh.size(); c++)
        if (std::fabs(cand_hh[c] - h0) < tol) ohs.push_back(cand_oh[2 * c]), ohs.push_back(cand_oh[2 * c + 1]);
    if (ohs.empty()) return;
    const double r0 = mode(ohs);
    // rigidity test: rigid water puts nearly all nearby candidates within `tight` of both modes
    size_t ntight = 0, nloose = 0, ns = 0;
    double sr = 0, sh = 0;
    for (size_t c = 0; c < cand_hh.size(); c++) {
        const double dh = std::fabs(cand_hh[c] - h0), da = std::fabs(cand_oh[2 * c] - r0),
                     db = std::fabs(cand_oh[2 * c + 1] - r0);
        if (dh < loose && da < loose && db < loose) {
            nloose++;
            ntight += dh < tight && da < tight && db < tight;
        }
        if (dh < tol && da < tol && db < tol) sr += double(cand_oh[2 * c]) + cand_oh[2 * c + 1], sh += cand_hh[c], ns++;
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
    L.r = sr / (2.0 * ns);
    L.hh = sh / double(ns);
}

// Bonds of the other atoms on frame x: each takes the nearest of its previous MAXOFF atoms if that is within bond
// range, and the bond-length classes are the peaks of those distances (1e-4 nm bins over 0.01 .. 0.25 nm; the virtual
// site of 4-site water is a 0.015 nm bond to its O).
template <class T>
void detect_bonds(const T *x, size_t N, Layout &L) {
    constexpr size_t NB = 2400;
    L.ref.assign(N, 0);
    L.blen.clear();
    std::vector<uint8_t> dof(N, 0);
    std::vector<uint16_t> bin(N, uint16_t(NB));
    std::vector<int64_t> hist(NB, 0);
    for (size_t i = 1; i < N; i++) {
        if (L.kind[i] != 2) continue;
        float best = 1e30f;
        unsigned off = 0;  // locals and no branch: a byte store may alias x, and which one is nearest is random
        for (unsigned o = 1; o <= unsigned(MAXOFF) && o <= i; o++) {
            const float dx = float(x[3 * i] - x[3 * (i - o)]), dy = float(x[3 * i + 1] - x[3 * (i - o) + 1]),
                        dz = float(x[3 * i + 2] - x[3 * (i - o) + 2]), v = dx * dx + dy * dy + dz * dz;
            off = v < best ? o : off;
            best = std::min(best, v);
        }
        dof[i] = uint8_t(off);
        if (best > 0.0001f && best < 0.0625f) {
            bin[i] = uint16_t(std::min(NB - 1, size_t((std::sqrt(best) - 0.01f) * 1e4f)));
            hist[bin[i]]++;
        }
    }
    // peaks: the heaviest +-10-bin window, then suppress +-30 bins; classes under 0.5% of the bonded atoms are ignored.
    // Bins within 0.004 nm (40 bins) of a class go to the nearest one.
    std::vector<int8_t> bincls(NB + 1, -1);
    std::vector<uint8_t> bdist(NB, 255);
    std::vector<int64_t> hs = hist, pre(NB + 1);
    int64_t total = 0;
    for (int64_t h : hist) total += h;
    for (int c = 0; c < MAXCLS; c++) {
        for (size_t k = 0; k < NB; k++) pre[k + 1] = pre[k] + hs[k];
        size_t best = 0;
        int64_t bv = 0;
        for (size_t k = 0; k < NB; k++) {
            const int64_t w = pre[std::min(NB, k + 11)] - pre[k >= 10 ? k - 10 : 0];
            if (w > bv) bv = w, best = k;
        }
        if (bv < 8 || bv * 200 < total) break;
        double sw = 0, sm = 0;
        for (size_t w = best >= 10 ? best - 10 : 0; w < std::min(NB, best + 11); w++)
            sw += hist[w], sm += hist[w] * (0.01 + (w + 0.5) * 1e-4);
        L.blen.push_back(sm / sw);
        for (size_t w = best >= 30 ? best - 30 : 0; w < std::min(NB, best + 31); w++) hs[w] = 0;
        const long centre = long((L.blen.back() - 0.01) * 1e4);
        for (long k = std::max(0L, centre - 40); k <= std::min(long(NB) - 1, centre + 40); k++)
            if (uint8_t(std::labs(k - centre)) < bdist[k])
                bdist[k] = uint8_t(std::labs(k - centre)), bincls[k] = int8_t(c);
    }
    for (size_t i = 0; i < N; i++)
        if (bincls[bin[i]] >= 0) L.ref[i] = uint16_t(dof[i] * 16 + bincls[bin[i]] + 1);
}

// ------------------------------------------------------------------------------------------------ symbols
// a buffer that grows by doubling and is never cleared (every value is written before it is read)
struct Grow {
    std::unique_ptr<int[]> p;
    size_t cap = 0;
    int *fit(size_t len, size_t need) {  // room for need values, keeping the first len
        if (cap < need) {
            std::unique_ptr<int[]> q(new int[2 * need]);
            if (len) std::copy(p.get(), p.get() + len, q.get());
            p = std::move(q);
            cap = 2 * need;
        }
        return p.get();
    }
};
struct Out {  // the frame's buffer of each stream, and raw bits (LSB first)
    int *p[NS];
    std::vector<uchar> *r;
    uint64_t acc = 0;
    int n = 0;
    ALWAYS_INLINE void sym(int s, uint32_t x) { *p[s]++ = int(x); }
    ALWAYS_INLINE void raw(uint64_t v, int b) {  // b <= 56
        acc |= uint64_t(v) << n;
        for (n += b; n >= 8; n -= 8, acc >>= 8) r->push_back(uchar(acc));
    }
};
struct Est {  // bits on a sample, for the predictor choice
    double bits = 0;
    inline void sym(int, uint32_t v) { bits += std::log2(double(v) + 1.0) + 1; }
    inline void raw(uint64_t, int b) { bits += b; }
};
struct In {  // a cursor per stream, and one on the raw bits
    const int *p[NS], *end[NS];
    const uchar *r, *rend;
    uint64_t acc = 0;
    int n = 0;
    ALWAYS_INLINE uint32_t sym(int s) {
        if (p[s] == end[s]) throw std::runtime_error("SZ3 BioMD: corrupt stream");
        return uint32_t(*p[s]++);
    }
    ALWAYS_INLINE uint64_t raw(int b) {  // b <= 56
        for (; n < b; n += 8) {
            if (r == rend) throw std::runtime_error("SZ3 BioMD: corrupt stream");
            acc |= uint64_t(*r++) << n;
        }
        const uint64_t v = acc & ((uint64_t(1) << b) - 1);
        acc >>= b;
        n -= b;
        return v;
    }
};

// v as one symbol below 4096, else as 4096 + v / 256 and its low byte as raw bits
template <class S>
ALWAYS_INLINE void sym_big(S &k, int s, uint32_t v) {
    if (v < 4096) return k.sym(s, v);
    k.sym(s, 4096 + (v >> 8));
    k.raw(v & 255, 8);
}
ALWAYS_INLINE uint32_t get_big(In &k, int s) {
    const uint32_t v = k.sym(s);
    return v < 4096 ? v : (v - 4096) << 8 | uint32_t(k.raw(8));
}

struct Frame {
    int32_t *q, *qp;  // this frame and the one before, on the lattice
    const Layout *L;
    int64_t R2, D2;  // squared O-H and H-H distances of water (lattice units)
    const int64_t *BR2;
    int32_t omin[3];  // the box of the water O's: corner, sides, and the bits of a point in it (x + Rx (y + Ry z))
    uint64_t R[3];
    int ob;  // or, if that takes more than 56, -1: then each coordinate in its own bits
};

inline int bits_of(uint64_t v) {  // v < 2^63
    int b = 0;
    while (v >> b) b++;
    return b;
}
inline void box(Frame &C, const uint32_t *side) {  // sides of at most 2^30
    for (int c = 0; c < 3; c++) C.R[c] = side[c];
    C.ob = C.R[0] * C.R[1] < (uint64_t(1) << 56) / C.R[2] ? bits_of(C.R[0] * C.R[1] * C.R[2] - 1) : -1;
}

// a water O of an intra frame: its point in the box
template <class S>
inline void enc_box(const Frame &C, const int32_t *p, S &k) {
    uint64_t u[3];
    for (int c = 0; c < 3; c++) u[c] = uint64_t(int64_t(p[c]) - C.omin[c]);
    if (C.ob >= 0) return k.raw(u[0] + C.R[0] * (u[1] + C.R[1] * u[2]), C.ob);
    for (int c = 0; c < 3; c++) k.raw(u[c], bits_of(C.R[c] - 1));
}
inline void dec_box(const Frame &C, int32_t *p, In &k) {
    uint64_t v = C.ob >= 0 ? k.raw(C.ob) : 0;
    for (int c = 0; c < 3; c++) {
        const uint64_t u = C.ob >= 0 ? v % C.R[c] : k.raw(bits_of(C.R[c] - 1));
        v /= C.R[c];
        p[c] = int32_t(C.omin[c] + int64_t(u));
    }
}

ALWAYS_INLINE int64_t pred_NU(const Frame &C, int m, size_t i, int c) {
    return m ? C.qp[3 * i + c] : (i ? C.q[3 * (i - 1) + c] : 0);
}

// Atom i on the sphere |d|^2 = R2 around atom j: the largest axis f of d is dropped and rebuilt from the other two,
// with a residual e. Intra (m == 0) the face (2 f + sign) is sent with e and the two kept coordinates as they are;
// else the vector p predicted by mode m fixes the axis and the sign, and the kept coordinates are residuals against p.
template <class S>
ALWAYS_INLINE void enc_sph(const Frame &C, int m, size_t i, size_t j, int64_t R2, int sfe, int sk, int se, S &k) {
    int64_t d[3], p[3] = {0, 0, 0};
    int f, x, y;
    rel(C.q, i, j, d);
    if (m) rel(C.qp, i, j, p);  // the vector of the previous frame
    axes(m ? p : d, f, x, y);
    const bool neg = (m ? p[f] : d[f]) < 0;
    const uint32_t ze = zz((neg ? -d[f] : d[f]) - isqrt_round(R2 - d[x] * d[x] - d[y] * d[y]));
    if (m == 0) k.sym(sfe, uint32_t(2 * f + neg) + 6 * std::min(ze, 15u));
    if (m != 0 || ze >= 15) k.sym(se, m == 0 ? ze - 15 : ze);
    k.sym(sk, zz(d[x] - p[x]));
    k.sym(sk, zz(d[y] - p[y]));
}
ALWAYS_INLINE void dec_sph(const Frame &C, int m, size_t i, size_t j, int64_t R2, int sfe, int sk, int se, In &k) {
    int64_t d[3], p[3] = {0, 0, 0};
    int f, x, y;
    bool neg;
    uint32_t ze;
    if (m == 0) {
        const uint32_t v = k.sym(sfe);
        f = int(v % 6 / 2);
        x = f == 2 ? 0 : f + 1;
        y = f == 0 ? 2 : (f == 1 ? 0 : 1);
        neg = v & 1;
        ze = v / 6;
        if (ze == 15) ze += k.sym(se);
    } else {
        rel(C.qp, i, j, p);  // the vector of the previous frame
        axes(p, f, x, y);
        neg = p[f] < 0;
        ze = k.sym(se);
    }
    d[x] = p[x] + unzz(k.sym(sk));
    d[y] = p[y] + unzz(k.sym(sk));
    const int64_t r = isqrt_round(R2 - d[x] * d[x] - d[y] * d[y]) + unzz(ze);
    d[f] = neg ? -r : r;
    for (int c = 0; c < 3; c++) C.q[3 * i + c] = int32_t(C.q[3 * j + c] + d[c]);
}

// H2 - O = d on the circle |d|^2 = R2, 2 d.d1 = R2 + |d1|^2 - D2, d1 = H1 - O. The coordinate a of d on the axis k
// where d1 is smallest leaves two points of the circle, one on each side of the plane through d1 and axis k; the side
// of d (intra: sent, else as a change from the side of the prediction p) picks one, and two residuals correct it.
// Integer products and a square root and a division of exact doubles, so the decoder finds the same point. Beyond the
// lattice sizes of a water (a molecule split by the boundary) the point is {a on k, 0, 0}.
ALWAYS_INLINE int circ_axis(const int64_t d1[3]) {
    int k = std::llabs(d1[1]) < std::llabs(d1[0]) ? 1 : 0;
    return std::llabs(d1[2]) < std::llabs(d1[k]) ? 2 : k;
}
// the side of v: 1 when (-d1[j], d1[i]) . (v[i], v[j]) < 0
ALWAYS_INLINE int circ_side(const int64_t d1[3], int k, const int64_t v[3]) {
    const int i = k == 2 ? 0 : k + 1, j = k == 0 ? 2 : (k == 1 ? 0 : 1);
    return d1[i] * v[j] - d1[j] * v[i] < 0;
}
ALWAYS_INLINE void circ_pt(const int64_t d1[3], int64_t R2, int64_t D2, int k, int64_t a, int r, int64_t P[3]) {
    const int i = k == 2 ? 0 : k + 1, j = k == 0 ? 2 : (k == 1 ? 0 : 1);
    const int64_t u = d1[i], w = d1[j], g2 = u * u + w * w, lim = 16384;
    const bool ok =
        g2 > 0 && std::llabs(u) < lim && std::llabs(w) < lim && std::llabs(d1[k]) < lim && std::llabs(a) < lim;
    const int64_t M = ok ? R2 + g2 + d1[k] * d1[k] - D2 - 2 * d1[k] * a : 0;  // 2 (u P_i + w P_j)
    const int64_t s = ok ? (r ? -1 : 1) * isqrt_round(4 * g2 * (R2 - a * a) - M * M) : 0;
    P[k] = a;
    P[i] = ok ? rnd(double(M * u - s * w) / double(2 * g2)) : 0;
    P[j] = ok ? rnd(double(M * w + s * u) / double(2 * g2)) : 0;
}
template <class S>
ALWAYS_INLINE void enc_circ(const Frame &C, int m, size_t o, S &k) {
    int64_t d1[3], d[3], p[3], P[3];
    rel(C.q, o + 1, o, d1);
    rel(C.q, o + 2, o, d);
    const int ax = circ_axis(d1), r = circ_side(d1, ax, d);
    if (m) {
        rel(C.qp, o + 2, o, p);
        k.sym(S_CA, zz(d[ax] - p[ax]) * 2 + uint32_t(r != circ_side(d1, ax, p)));
    } else {
        k.sym(S_CA, zz(d[ax]) * 2 + uint32_t(r));
    }
    circ_pt(d1, C.R2, C.D2, ax, d[ax], r, P);
    for (int c = 0; c < 3; c++)
        if (c != ax) k.sym(S_CE, zz(d[c] - P[c]));
}
ALWAYS_INLINE void dec_circ(const Frame &C, int m, size_t o, In &k) {
    int64_t d1[3], p[3] = {0, 0, 0}, P[3];
    rel(C.q, o + 1, o, d1);
    const int ax = circ_axis(d1);
    if (m) rel(C.qp, o + 2, o, p);
    const uint32_t v = k.sym(S_CA);
    circ_pt(d1, C.R2, C.D2, ax, p[ax] + unzz(v >> 1), int(v & 1) ^ (m ? circ_side(d1, ax, p) : 0), P);
    for (int c = 0; c < 3; c++)
        C.q[3 * (o + 2) + c] = int32_t(C.q[3 * o + c] + (c == ax ? P[c] : P[c] + unzz(k.sym(S_CE))));
}

// One atom: a water (O and the rest of its molecule), or another atom. kind 1 atoms come with their O.
template <class S>
ALWAYS_INLINE void enc_atom(const Frame &C, const int *mode, size_t i, S &k) {
    const Layout &L = *C.L;
    const int32_t *q = C.q;
    if (L.kind[i] == 0) {
        if (mode[G_O])
            for (int c = 0; c < 3; c++) k.sym(S_O, zz(q[3 * i + c] - C.qp[3 * i + c]));
        else
            enc_box(C, q + 3 * i, k);
        enc_sph(C, mode[G_WH], i + 1, i, C.R2, S_FE, S_KEPT, S_E, k);
        enc_circ(C, mode[G_WH], i, k);
    } else if (L.kind[i] == 2) {
        if (L.ref[i]) {
            enc_sph(C, mode[G_NB], i, i - (L.ref[i] >> 4), C.BR2[(L.ref[i] & 15) - 1], S_BFE, S_BK, S_BE, k);
        } else {
            for (int c = 0; c < 3; c++) sym_big(k, S_U, zz(q[3 * i + c] - pred_NU(C, mode[G_NU], i, c)));
        }
    }
}
ALWAYS_INLINE void dec_atom(const Frame &C, const int *mode, size_t i, In &k) {
    const Layout &L = *C.L;
    int32_t *q = C.q;
    if (L.kind[i] == 0) {
        if (mode[G_O])
            for (int c = 0; c < 3; c++) q[3 * i + c] = int32_t(C.qp[3 * i + c] + unzz(k.sym(S_O)));
        else
            dec_box(C, q + 3 * i, k);
        dec_sph(C, mode[G_WH], i + 1, i, C.R2, S_FE, S_KEPT, S_E, k);
        dec_circ(C, mode[G_WH], i, k);
    } else if (L.kind[i] == 2) {
        if (L.ref[i]) {
            dec_sph(C, mode[G_NB], i, i - (L.ref[i] >> 4), C.BR2[(L.ref[i] & 15) - 1], S_BFE, S_BK, S_BE, k);
        } else {
            for (int c = 0; c < 3; c++) q[3 * i + c] = int32_t(pred_NU(C, mode[G_NU], i, c) + unzz(get_big(k, S_U)));
        }
    }
}

}  // namespace biomd

template <class T, uint N>
class SZBioMDDecomposition : public concepts::DecompositionInterface<T, int, N> {
   public:
    explicit SZBioMDDecomposition(const Config &conf)
        : F_(N == 3 ? conf.dims[0] : 1), A_(N >= 2 ? conf.dims[N - 2] : 1), eb_(conf.absErrorBound) {
        if (N > 3 || conf.dims[N - 1] != 3) throw std::invalid_argument("SZ3 BioMD: data must be {frames, atoms, 3}");
        if (!std::is_floating_point<T>::value) throw std::invalid_argument("SZ3 BioMD: data must be float or double");
    }

    std::vector<int> compress(const Config & /*conf*/, T *data) override {
        using namespace biomd;
        // trailing frames all of one value (the unwritten rest of a chunk) are stored as that value
        fill_ = data[(F_ - 1) * A_ * 3];
        auto is_fill = [this, data](size_t t) {
            for (size_t i = t * A_ * 3; i < (t + 1) * A_ * 3; i++)
                if (memcmp(&data[i], &fill_, sizeof(T)) != 0) return false;
            return true;
        };
        for (Fc_ = F_; Fc_ > 1 && is_fill(Fc_ - 1);) Fc_--;
        // lattice step: 2 (eb - ulp) for the binade of max|x|, so that q step rounded to T stays within eb
        using U = typename std::conditional<sizeof(T) == 4, uint32_t, uint64_t>::type;
        U bmax = 0;  // the largest |x| as bits: NaN and inf are above every finite value
        for (size_t i = 0; i < Fc_ * A_ * 3; i++) {
            U b;
            memcpy(&b, &data[i], sizeof(T));
            b &= U(~U(0)) >> 1;
            bmax = b > bmax ? b : bmax;
        }
        T amax;
        memcpy(&amax, &bmax, sizeof(T));
        if (!(amax <= std::numeric_limits<T>::max())) throw std::runtime_error("SZ3 BioMD: non-finite input");
        int e2;
        std::frexp(double(amax), &e2);
        const double ulp = std::ldexp(1.0, e2 - std::numeric_limits<T>::digits);
        step_ = ulp < 0.5 * eb_ ? 2.0 * (eb_ - ulp) : eb_;
        if (!(step_ > 0) || amax / step_ > double(1 << 28))
            throw std::runtime_error("SZ3 BioMD: error bound too small for the coordinate range");
        const double inv = 1.0 / step_;

        // the layout, from the first frame
        Layout L;
        detect_water(data, A_, L, step_);
        R2_ = D2_ = 0;
        if (L.r > 0 && L.r / step_ < 16384) {
            R2_ = rnd((L.r / step_) * (L.r / step_));
            D2_ = rnd((L.hh / step_) * (L.hh / step_));
            D2_ = D2_ < (int64_t(1) << 30) ? D2_ : 0;
        } else {  // no water, or too many lattice steps across one for exact products
            std::fill(L.kind.begin(), L.kind.end(), 2);
        }
        detect_bonds(data, A_, L);
        BR2_.clear();
        for (double b : L.blen) BR2_.push_back(rnd((b / step_) * (b / step_)));
        std::vector<uint32_t> wo, ub[2];  // water O's; bonded, unbonded other atoms
        for (size_t i = 0; i < A_; i++) {
            if (L.ref[i] && BR2_[(L.ref[i] & 15) - 1] >= RMAXB2) L.ref[i] = 0;
            if (L.kind[i] == 0) wo.push_back(uint32_t(i));
            if (L.kind[i] == 2) ub[L.ref[i] ? 0 : 1].push_back(uint32_t(i));
        }
        nwat_ = wo.size();

        // per (stream, mode) a buffer, kept long enough for the next frame's symbols
        biomd::Grow buf_s[NS * 2];
        size_t len[NS * 2] = {0};
        const size_t nunit[NGROUP + 1] = {wo.size(), wo.size(), ub[0].size(), ub[1].size() + ub[0].size(), A_};
        Out o;
        o.r = &raw_;
        raw_.clear();
        auto open = [&](const int *mode) {
            for (int s = 0; s < NS; s++) {
                const size_t k = s * 2 + (SGROUP[s] < 0 ? 0 : mode[SGROUP[s]]);
                const size_t need = len[k] + PER_UNIT[s] * nunit[SGROUP[s] < 0 ? NGROUP : SGROUP[s]] + 16;
                o.p[s] = buf_s[k].fit(len[k], need) + len[k];
            }
        };
        auto close = [&](const int *mode) {
            for (int s = 0; s < NS; s++) {
                const size_t k = s * 2 + (SGROUP[s] < 0 ? 0 : mode[SGROUP[s]]);
                len[k] = size_t(o.p[s] - buf_s[k].p.get());
            }
        };
        const int mode0[NGROUP] = {0, 0, 0, 0};
        open(mode0);
        for (size_t k = 0; k < wo.size(); k++) o.sym(S_W, wo[k] - (k ? wo[k - 1] : 0));  // the waters as gaps
        for (size_t i = 0; i < A_; i++)
            if (L.kind[i] == 2) o.sym(S_REF, L.ref[i]);
        close(mode0);
        std::unique_ptr<int32_t[]> buf(new int32_t[6 * A_]);  // written before read
        Frame C{buf.get(), buf.get() + 3 * A_, &L, R2_, D2_, BR2_.data(), {0, 0, 0}, {1, 1, 1}, 0};
        const std::vector<uint32_t> *units[NGROUP] = {&wo, &wo, &ub[0], &ub[1]};
        int mode[NGROUP] = {0, 0, 0, 0};
        modes_.clear();
        omin_.clear();
        oside_.clear();
        for (size_t t = 0; t < Fc_; t++) {
            const T *X = data + t * A_ * 3;
            for (size_t i = 0; i < A_ * 3; i++) {
                const double y = double(X[i]) * inv;
                C.q[i] = int32_t(y + std::copysign(0.5, y));
            }
            int32_t mx[3] = {0, 0, 0};
            for (int c = 0; c < 3; c++) C.omin[c] = mx[c] = wo.empty() ? 0 : C.q[3 * wo[0] + c];
            for (uint32_t i : wo)
                for (int c = 0; c < 3; c++) {
                    C.omin[c] = std::min(C.omin[c], C.q[3 * i + c]);
                    mx[c] = std::max(mx[c], C.q[3 * i + c]);
                }
            // per group, the predictor that is cheapest on a sample: at frames 1 and 2, then every 8th frame
            if (t > 0 && (t <= 2 || t % 8 == 0))
                for (int g = 0; g < NGROUP; g++) {
                    const auto &u = *units[g];
                    double best = 1e300;
                    for (int m = 0; m < 2 && !u.empty(); m++) {
                        Est est;
                        int md[NGROUP] = {mode[0], mode[1], mode[2], mode[3]};
                        md[g] = m;
                        const size_t st = std::max<size_t>(1, u.size() / 256);
                        for (size_t k = 0; k < u.size(); k += st)
                            if (g != G_O)
                                enc_atom(C, md, u[k], est);
                            else if (m == 0)  // the box index
                                for (int c = 0; c < 3; c++) est.bits += std::log2(double(mx[c] - C.omin[c]) + 1.0);
                            else
                                for (int c = 0; c < 3; c++) est.sym(S_O, zz(C.q[3 * u[k] + c] - C.qp[3 * u[k] + c]));
                        if (est.bits < best) best = est.bits, mode[g] = m;
                    }
                }
            if (t == 0) std::fill(mode, mode + NGROUP, 0);
            for (int g = 0; g < NGROUP; g++) modes_.push_back(uint8_t(mode[g]));
            for (int c = 0; c < 3; c++) {
                omin_.push_back(C.omin[c]);
                oside_.push_back(uint32_t(mx[c] - C.omin[c]) + 1);
            }
            box(C, &oside_[t * 3]);
            open(mode);
            for (size_t i = 0; i < A_; i++) enc_atom(C, mode, i, o);
            close(mode);
            std::swap(C.q, C.qp);
        }
        if (o.n > 0) raw_.push_back(uchar(o.acc));
        // the streams one after the other, each as [count, symbols], for SegmentedEncoder
        std::vector<int> seg;
        size_t n = NS * 2;
        for (size_t k = 0; k < NS * 2; k++) n += len[k];
        seg.reserve(n);
        for (size_t k = 0; k < NS * 2; k++) {
            seg.push_back(int(len[k]));
            seg.insert(seg.end(), buf_s[k].p.get(), buf_s[k].p.get() + len[k]);
        }
        return seg;
    }

    T *decompress(const Config & /*conf*/, std::vector<int> &quant_inds, T *dec_data) override {
        using namespace biomd;
        Layout L;
        const int *cur[NS * 2], *end[NS * 2];
        for (size_t k = 0, s = 0; s < NS * 2; s++) {
            const size_t n = k < quant_inds.size() ? size_t(quant_inds[k]) : 0;
            if (k >= quant_inds.size() || n > quant_inds.size() - k - 1)
                throw std::runtime_error("SZ3 BioMD: corrupt stream");
            cur[s] = quant_inds.data() + k + 1;
            end[s] = cur[s] + n;
            k += n + 1;
        }
        In in;
        auto open = [&](const int *mode) {
            for (int s = 0; s < NS; s++) {
                const size_t k = s * 2 + (SGROUP[s] < 0 ? 0 : mode[SGROUP[s]]);
                in.p[s] = cur[k];
                in.end[s] = end[k];
            }
        };
        auto close = [&](const int *mode) {
            for (int s = 0; s < NS; s++) cur[s * 2 + (SGROUP[s] < 0 ? 0 : mode[SGROUP[s]])] = in.p[s];
        };
        const int mode0[NGROUP] = {0, 0, 0, 0};
        open(mode0);
        L.kind.assign(A_, 2);
        for (size_t k = 0, o = 0; k < nwat_; k++) {
            const uint32_t g = in.sym(S_W);
            if ((k && g < 3) || g > A_ || o + g + 2 >= A_) throw std::runtime_error("SZ3 BioMD: corrupt stream");
            o += g;
            L.kind[o] = 0;
            L.kind[o + 1] = L.kind[o + 2] = 1;
        }
        L.ref.assign(A_, 0);
        for (size_t i = 0; i < A_; i++)
            if (L.kind[i] == 2) {
                const uint32_t r = in.sym(S_REF);
                if (r && ((r >> 4) > i || (r & 15) == 0 || (r & 15) > BR2_.size()))
                    throw std::runtime_error("SZ3 BioMD: corrupt stream");
                L.ref[i] = uint16_t(r);
            }
        close(mode0);
        std::unique_ptr<int32_t[]> buf(new int32_t[6 * A_]);  // written before read
        Frame C{buf.get(), buf.get() + 3 * A_, &L, R2_, D2_, BR2_.data(), {0, 0, 0}, {1, 1, 1}, 0};
        in.r = raw_.data();
        in.rend = raw_.data() + raw_.size();
        for (size_t t = 0; t < Fc_; t++) {
            int mode[NGROUP];
            for (int g = 0; g < NGROUP; g++) mode[g] = modes_[t * NGROUP + g] & 1;
            for (int c = 0; c < 3; c++) C.omin[c] = omin_[t * 3 + c];
            box(C, &oside_[t * 3]);
            open(mode);
            for (size_t i = 0; i < A_; i++) dec_atom(C, mode, i, in);
            close(mode);
            T *D = dec_data + t * A_ * 3;
            for (size_t i = 0; i < A_ * 3; i++) D[i] = T(double(C.q[i]) * step_);
            std::swap(C.q, C.qp);
        }
        std::fill(dec_data + Fc_ * A_ * 3, dec_data + F_ * A_ * 3, fill_);
        return dec_data;
    }

    void save(uchar *&c) override {
        write(uint64_t(Fc_), c);
        write(fill_, c);
        write(step_, c);
        write(R2_, c);
        write(D2_, c);
        write(uint64_t(nwat_), c);
        write(uint8_t(BR2_.size()), c);
        if (!BR2_.empty()) write(BR2_.data(), BR2_.size(), c);
        write(modes_.data(), modes_.size(), c);
        write(omin_.data(), omin_.size(), c);
        write(oside_.data(), oside_.size(), c);
        write(uint64_t(raw_.size()), c);
        if (!raw_.empty()) write(raw_.data(), raw_.size(), c);
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        uint8_t nb;
        uint64_t nw, fc;
        read(fc, c, remaining_length);
        if (fc < 1 || fc > F_) throw std::runtime_error("SZ3 BioMD: corrupt stream");
        Fc_ = size_t(fc);
        read(fill_, c, remaining_length);
        read(step_, c, remaining_length);
        read(R2_, c, remaining_length);
        read(D2_, c, remaining_length);
        read(nw, c, remaining_length);
        nwat_ = size_t(std::min<uint64_t>(nw, A_));
        read(nb, c, remaining_length);
        BR2_.resize(nb);
        if (nb) read(BR2_.data(), nb, c, remaining_length);
        modes_.resize(Fc_ * biomd::NGROUP);
        read(modes_.data(), modes_.size(), c, remaining_length);
        omin_.resize(Fc_ * 3);
        read(omin_.data(), omin_.size(), c, remaining_length);
        oside_.resize(Fc_ * 3);
        read(oside_.data(), oside_.size(), c, remaining_length);
        uint64_t nr = 0;
        read(nr, c, remaining_length);
        if (nr > remaining_length) throw std::runtime_error("SZ3 BioMD: corrupt stream");
        raw_.resize(size_t(nr));
        if (nr) read(raw_.data(), raw_.size(), c, remaining_length);
        bool ok = nw <= A_ && R2_ >= 0 && R2_ < (int64_t(1) << 28) && D2_ >= 0 && D2_ < (int64_t(1) << 30) &&
                  nb <= biomd::MAXCLS;
        for (int64_t b : BR2_) ok = ok && b >= 0;  // classes of RMAXB2 or more are stored, no atom refers to them
        for (uint32_t r : oside_) ok = ok && r >= 1 && r <= (uint32_t(1) << 30);
        if (!ok) throw std::runtime_error("SZ3 BioMD: corrupt stream");
    }

    size_t size_est() override { return 64 + 8 * BR2_.size() + 24 * F_ + raw_.size(); }

    std::pair<int, int> get_out_range() override { return {0, 0}; }  // SegmentedEncoder takes any int

   private:
    size_t F_, A_, Fc_ = 1;  // frames, atoms, frames before the fill
    T fill_ = 0;
    double eb_, step_ = 0;
    int64_t R2_ = 0, D2_ = 0;
    size_t nwat_ = 0;
    std::vector<int64_t> BR2_;
    std::vector<uint8_t> modes_;   // per frame and group
    std::vector<int32_t> omin_;    // per frame: the corner of the water O box
    std::vector<uint32_t> oside_;  // per frame: its sides
    std::vector<uchar> raw_;       // the water O's of intra frames
};

template <class T, uint N>
SZBioMDDecomposition<T, N> make_decomposition_biomd(const Config &conf) {
    return SZBioMDDecomposition<T, N>(conf);
}

}  // namespace SZ3
#endif
