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
enum { S_W, S_REF, S_O, S_OH, S_OL, S_FE, S_KEPT, S_E, S_CA, S_CE, S_BFE, S_BK, S_BE, S_U, NS };
constexpr int NMODE[NGROUP] = {2, 2, 2, 2};
// each stream is split by the mode of its group, so each (stream, mode) gets its own code
constexpr int SGROUP[NS] = {-1, -1, G_O, G_O, G_O, G_WH, G_WH, G_WH, G_WH, G_WH, G_NB, G_NB, G_NB, G_NU};
constexpr int PER_UNIT[NS] = {1, 1, 3, 3, 3, 2, 4,
                              2, 1, 2, 1, 3, 1, 3};  // most symbols per unit of the group, per frame

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
        for (size_t o = 1; o <= size_t(MAXOFF) && o <= i; o++) {
            const float dx = float(x[3 * i] - x[3 * (i - o)]), dy = float(x[3 * i + 1] - x[3 * (i - o) + 1]),
                        dz = float(x[3 * i + 2] - x[3 * (i - o) + 2]), v = dx * dx + dy * dy + dz * dz;
            if (v < best) best = v, dof[i] = uint8_t(o);
        }
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
struct Out {  // the frame's buffer of each stream
    int *p[NS];
    ALWAYS_INLINE void sym(int s, uint32_t x) { *p[s]++ = int(x); }
};
struct Est {  // bits on a sample, for the predictor choice
    double bits = 0;
    inline void sym(int, uint32_t v) { bits += std::log2(double(v) + 1.0) + 1; }
};
struct In {  // a cursor per stream
    const int *p[NS], *end[NS];
    ALWAYS_INLINE uint32_t sym(int s) {
        if (p[s] == end[s]) throw std::runtime_error("SZ3 BioMD: corrupt stream");
        return uint32_t(*p[s]++);
    }
};

struct Frame {
    int32_t *q, *qp;  // this frame and the one before, on the lattice
    const Layout *L;
    int64_t R2, D2;  // squared O-H and H-H distances of water (lattice units)
    const int64_t *BR2;
    int32_t omin[3];
};

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
// where d1 is smallest fixes two points P of the circle (integer arithmetic and one square root, so the decoder finds
// the same); the one meant (intra: a bit, else the one nearer the prediction p) is corrected by two residuals.
// Beyond the lattice sizes of a water (a molecule split by the boundary) both points are {a on k, 0, 0}.
ALWAYS_INLINE int64_t rdiv(int64_t n, int64_t d) { return (n >= 0 ? n + d / 2 : n - d / 2) / d; }  // d > 0
ALWAYS_INLINE void circ_pts(const int64_t d1[3], int64_t R2, int64_t D2, int k, int64_t a, int64_t P[2][3]) {
    const int i = k == 2 ? 0 : k + 1, j = k == 0 ? 2 : (k == 1 ? 0 : 1);
    const int64_t u = d1[i], w = d1[j], g2 = u * u + w * w, lim = 16384;
    const bool ok =
        g2 > 0 && std::llabs(u) < lim && std::llabs(w) < lim && std::llabs(d1[k]) < lim && std::llabs(a) < lim;
    const int64_t M = ok ? R2 + g2 + d1[k] * d1[k] - D2 - 2 * d1[k] * a : 0;  // 2 (u P_i + w P_j)
    const int64_t s = ok ? isqrt_round(4 * g2 * (R2 - a * a) - M * M) : 0;
    for (int r = 0; r < 2; r++) {
        P[r][k] = a;
        P[r][i] = ok ? rdiv(M * u - (r ? -s : s) * w, 2 * g2) : 0;
        P[r][j] = ok ? rdiv(M * w + (r ? -s : s) * u, 2 * g2) : 0;
    }
}
ALWAYS_INLINE int circ_axis(const int64_t d1[3]) {
    int k = std::llabs(d1[1]) < std::llabs(d1[0]) ? 1 : 0;
    return std::llabs(d1[2]) < std::llabs(d1[k]) ? 2 : k;
}
ALWAYS_INLINE int64_t l1(const int64_t a[3], const int64_t b[3]) {
    return std::llabs(a[0] - b[0]) + std::llabs(a[1] - b[1]) + std::llabs(a[2] - b[2]);
}
template <class S>
ALWAYS_INLINE void enc_circ(const Frame &C, int m, size_t o, S &k) {
    int64_t d1[3], d[3], p[3], P[2][3];
    rel(C.q, o + 1, o, d1);
    rel(C.q, o + 2, o, d);
    const int ax = circ_axis(d1);
    circ_pts(d1, C.R2, C.D2, ax, d[ax], P);
    int r = l1(P[1], d) < l1(P[0], d);
    if (m) {
        rel(C.qp, o + 2, o, p);
        const int rp = l1(P[1], p) < l1(P[0], p);
        k.sym(S_CA, zz(d[ax] - p[ax]) * 2 + uint32_t(r != rp));
    } else {
        k.sym(S_CA, zz(d[ax]) * 2 + uint32_t(r));
    }
    for (int c = 0; c < 3; c++)
        if (c != ax) k.sym(S_CE, zz(d[c] - P[r][c]));
}
ALWAYS_INLINE void dec_circ(const Frame &C, int m, size_t o, In &k) {
    int64_t d1[3], p[3] = {0, 0, 0}, P[2][3];
    rel(C.q, o + 1, o, d1);
    const int ax = circ_axis(d1);
    if (m) rel(C.qp, o + 2, o, p);
    const uint32_t v = k.sym(S_CA);
    circ_pts(d1, C.R2, C.D2, ax, p[ax] + unzz(v >> 1), P);
    int r = int(v & 1);
    if (m) r ^= l1(P[1], p) < l1(P[0], p);
    for (int c = 0; c < 3; c++)
        C.q[3 * (o + 2) + c] = int32_t(C.q[3 * o + c] + (c == ax ? P[r][c] : P[r][c] + unzz(k.sym(S_CE))));
}

// One atom: a water (O and the rest of its molecule), or another atom. kind 1 atoms come with their O.
template <class S>
ALWAYS_INLINE void enc_atom(const Frame &C, const int *mode, size_t i, S &k) {
    const Layout &L = *C.L;
    const int32_t *q = C.q;
    if (L.kind[i] == 0) {
        for (int c = 0; c < 3; c++) {
            if (mode[G_O] == 0) {
                const uint32_t u = uint32_t(q[3 * i + c] - C.omin[c]);
                k.sym(S_OH, u >> 8);
                k.sym(S_OL, u & 255);
            } else {
                k.sym(S_O, zz(q[3 * i + c] - C.qp[3 * i + c]));
            }
        }
        enc_sph(C, mode[G_WH], i + 1, i, C.R2, S_FE, S_KEPT, S_E, k);
        enc_circ(C, mode[G_WH], i, k);
    } else if (L.kind[i] == 2) {
        if (L.ref[i]) {
            enc_sph(C, mode[G_NB], i, i - (L.ref[i] >> 4), C.BR2[(L.ref[i] & 15) - 1], S_BFE, S_BK, S_BE, k);
        } else {
            for (int c = 0; c < 3; c++) k.sym(S_U, zz(q[3 * i + c] - pred_NU(C, mode[G_NU], i, c)));
        }
    }
}
ALWAYS_INLINE void dec_atom(const Frame &C, const int *mode, size_t i, In &k) {
    const Layout &L = *C.L;
    int32_t *q = C.q;
    if (L.kind[i] == 0) {
        for (int c = 0; c < 3; c++) {
            if (mode[G_O] == 0) {
                const uint32_t hi = k.sym(S_OH);
                q[3 * i + c] = int32_t(C.omin[c] + int64_t(hi << 8 | k.sym(S_OL)));
            } else {
                q[3 * i + c] = int32_t(C.qp[3 * i + c] + unzz(k.sym(S_O)));
            }
        }
        dec_sph(C, mode[G_WH], i + 1, i, C.R2, S_FE, S_KEPT, S_E, k);
        dec_circ(C, mode[G_WH], i, k);
    } else if (L.kind[i] == 2) {
        if (L.ref[i]) {
            dec_sph(C, mode[G_NB], i, i - (L.ref[i] >> 4), C.BR2[(L.ref[i] & 15) - 1], S_BFE, S_BK, S_BE, k);
        } else {
            for (int c = 0; c < 3; c++) q[3 * i + c] = int32_t(pred_NU(C, mode[G_NU], i, c) + unzz(k.sym(S_U)));
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

        // layout of the first frame, or the previous chunk's
        Layout L;
        std::vector<uint32_t> wo, ub[2];  // water O's; bonded, unbonded other atoms
        std::vector<int64_t> BR2;
        int64_t R2 = 0, D2 = 0;
        {
            detect_water(data, A_, L, step_);
            R2 = 0;
            if (L.r > 0 && L.r / step_ < 16384) {
                R2 = rnd((L.r / step_) * (L.r / step_));
                D2 = rnd((L.hh / step_) * (L.hh / step_));
            } else {  // no water, or too many lattice steps across one for exact products
                std::fill(L.kind.begin(), L.kind.end(), 2);
            }
            detect_bonds(data, A_, L);
            BR2.clear();
            for (double b : L.blen) BR2.push_back(rnd((b / step_) * (b / step_)));
            wo.clear();
            for (auto &u : ub) u.clear();
            for (size_t i = 0; i < A_; i++) {
                if (L.ref[i] && BR2[(L.ref[i] & 15) - 1] >= RMAXB2) L.ref[i] = 0;
                if (L.kind[i] == 0) wo.push_back(uint32_t(i));
                if (L.kind[i] == 2) ub[L.ref[i] ? 0 : 1].push_back(uint32_t(i));
            }
        }
        R2_ = R2;
        D2_ = D2 < (int64_t(1) << 30) ? D2 : 0;
        BR2_ = BR2;
        nwat_ = wo.size();

        // per (stream, mode) a buffer, kept long enough for the next frame's symbols
        thread_local std::vector<int> buf_s[NS * 4];
        size_t len[NS * 4] = {0};
        const size_t nunit[NGROUP + 1] = {wo.size(), wo.size(), ub[0].size(), ub[1].size() + ub[0].size(), A_};
        Out o;
        auto open = [&](const int *mode) {
            for (int s = 0; s < NS; s++) {
                const size_t k = s * 4 + (SGROUP[s] < 0 ? 0 : mode[SGROUP[s]]);
                const size_t need = len[k] + PER_UNIT[s] * nunit[SGROUP[s] < 0 ? NGROUP : SGROUP[s]] + 16;
                if (buf_s[k].size() < need) buf_s[k].resize(2 * need);
                o.p[s] = buf_s[k].data() + len[k];
            }
        };
        auto close = [&](const int *mode) {
            for (int s = 0; s < NS; s++) {
                const size_t k = s * 4 + (SGROUP[s] < 0 ? 0 : mode[SGROUP[s]]);
                len[k] = size_t(o.p[s] - buf_s[k].data());
            }
        };
        const int mode0[NGROUP] = {0, 0, 0, 0};
        open(mode0);
        for (size_t k = 0; k < wo.size(); k++) o.sym(S_W, wo[k] - (k ? wo[k - 1] : 0));  // the waters as gaps
        for (size_t i = 0; i < A_; i++)
            if (L.kind[i] == 2) o.sym(S_REF, L.ref[i]);
        close(mode0);
        thread_local std::vector<int32_t> buf;
        if (buf.size() < 6 * A_) buf.resize(6 * A_);
        Frame C{buf.data(), buf.data() + 3 * A_, &L, R2, D2_, BR2_.data(), {0, 0, 0}};
        const std::vector<uint32_t> *units[NGROUP] = {&wo, &wo, &ub[0], &ub[1]};
        int mode[NGROUP] = {0, 0, 0, 0};
        modes_.clear();
        omin_.clear();
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
                    for (int m = 0; m < NMODE[g] && !u.empty(); m++) {
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
            for (int c = 0; c < 3; c++) omin_.push_back(C.omin[c]);
            open(mode);
            for (size_t i = 0; i < A_; i++) enc_atom(C, mode, i, o);
            close(mode);
            std::swap(C.q, C.qp);
        }
        // the streams one after the other, each as [count, symbols], for SegmentedEncoder
        std::vector<int> seg;
        size_t n = NS * 4;
        for (size_t k = 0; k < NS * 4; k++) n += len[k];
        seg.reserve(n);
        for (size_t k = 0; k < NS * 4; k++) {
            seg.push_back(int(len[k]));
            seg.insert(seg.end(), buf_s[k].begin(), buf_s[k].begin() + len[k]);
        }
        return seg;
    }

    T *decompress(const Config & /*conf*/, std::vector<int> &quant_inds, T *dec_data) override {
        using namespace biomd;
        Layout L;
        const int *cur[NS * 4], *end[NS * 4];
        for (size_t k = 0, s = 0; s < NS * 4; s++) {
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
                const size_t k = s * 4 + (SGROUP[s] < 0 ? 0 : mode[SGROUP[s]]);
                in.p[s] = cur[k];
                in.end[s] = end[k];
            }
        };
        auto close = [&](const int *mode) {
            for (int s = 0; s < NS; s++) cur[s * 4 + (SGROUP[s] < 0 ? 0 : mode[SGROUP[s]])] = in.p[s];
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
        thread_local std::vector<int32_t> buf;
        if (buf.size() < 6 * A_) buf.resize(6 * A_);
        Frame C{buf.data(), buf.data() + 3 * A_, &L, R2_, D2_, BR2_.data(), {0, 0, 0}};
        for (size_t t = 0; t < Fc_; t++) {
            int mode[NGROUP];
            for (int g = 0; g < NGROUP; g++) mode[g] = modes_[t * NGROUP + g] % NMODE[g];
            for (int c = 0; c < 3; c++) C.omin[c] = omin_[t * 3 + c];
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
        bool ok = nw <= A_ && R2_ >= 0 && R2_ < (int64_t(1) << 28) && D2_ >= 0 && D2_ < (int64_t(1) << 30) &&
                  nb <= biomd::MAXCLS;
        for (int64_t b : BR2_) ok = ok && b >= 0;  // classes of RMAXB2 or more are stored, no atom refers to them
        if (!ok) throw std::runtime_error("SZ3 BioMD: corrupt stream");
    }

    size_t size_est() override { return 64 + 8 * BR2_.size() + 16 * F_; }

    std::pair<int, int> get_out_range() override { return {0, 0}; }  // SegmentedEncoder takes any int

   private:
    size_t F_, A_, Fc_ = 1;  // frames, atoms, frames before the fill
    T fill_ = 0;
    double eb_, step_ = 0;
    int64_t R2_ = 0, D2_ = 0;
    size_t nwat_ = 0;
    std::vector<int64_t> BR2_;
    std::vector<uint8_t> modes_;  // per frame and group
    std::vector<int32_t> omin_;   // per frame: the corner of the water O box
};

template <class T, uint N>
SZBioMDDecomposition<T, N> make_decomposition_biomd(const Config &conf) {
    return SZBioMDDecomposition<T, N>(conf);
}

}  // namespace SZ3
#endif
