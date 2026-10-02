#ifndef SZ3_BIOMD_DECOMPOSITION_HPP
#define SZ3_BIOMD_DECOMPOSITION_HPP

// ALGO_BIOMD: molecular-dynamics coordinates {frames, atoms, 3} (nm, absolute bound) as streams of integer symbols,
// one after the other as [count, symbols]. Every coordinate goes on a lattice,
//   q = round(x / step), |x - q step| <= eb,
// and all prediction is integer arithmetic on that lattice, which the decoder repeats exactly.
//  * rigid water (O, H1, H2): H1 on the sphere |H1 - O| = r, H2 on the circle that r and the H-H distance leave
//  * other atoms: on the sphere of a bond-length class around one of the previous MAX_BOND_OFFSET atoms, else a delta
//  * frame 0 intra; for the frames after it, per atom group intra or from the previous frame, whichever costs less on a
//    sample of frame 1
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
#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {
namespace biomd {

constexpr int MAX_BOND_OFFSET = 4;                 // a bond partner is one of the previous MAX_BOND_OFFSET atoms
constexpr int MAX_BOND_CLASSES = 15;               // at most this many bond-length classes
constexpr int64_t MAX_BOND_R2 = int64_t(1) << 28;  // bonds of 16384 lattice units or more are not coded as spheres
constexpr uint32_t ESCAPE = 4096;                  // unbonded values from here on: an escape symbol and a raw byte
// |q| <= 2^28: displacements stay within 2^29, and every symbol within 32 bits
constexpr double MAX_LATTICE = double(1 << 28);
// the lattice spans |x| < B, the smallest power of two of at least 2^LATTICE_SPAN_BITS eb
#ifndef SZ3_BIOMD_LATTICE_SPAN_BITS
#define SZ3_BIOMD_LATTICE_SPAN_BITS 17
#endif
constexpr int LATTICE_SPAN_BITS = SZ3_BIOMD_LATTICE_SPAN_BITS;

// atom groups, each with its own predictor mode for the frames after the first
enum { G_WATER_O, G_WATER_H, G_BONDED, G_UNBONDED, NUM_GROUPS };
// symbol streams: the layout (gaps between water O's, the bond of each other atom); water O from the previous frame;
// water H1 on the sphere around O and bonded atoms on the sphere around their partner (face, kept coordinates, radial
// residual); water H2 on the circle (coordinate and side, residuals); unbonded atoms as deltas
enum {
    S_WATER_GAP,
    S_BOND_REF,
    S_WATER_O,
    S_WATER_H1_FACE,
    S_WATER_H1_KEPT,
    S_WATER_H1_RADIAL,
    S_WATER_H2_AXIS,
    S_WATER_H2_RESIDUAL,
    S_BOND_FACE,
    S_BOND_KEPT,
    S_BOND_RADIAL,
    S_UNBONDED,
    NUM_STREAMS
};
// the group whose atoms write the stream (-1: the layout streams, one symbol per atom at most), and the most symbols
// one of them (water: one molecule) writes per frame; the two size the streams' buffers
inline constexpr int STREAM_GROUP[NUM_STREAMS] = {-1,        -1,        G_WATER_O, G_WATER_H, G_WATER_H, G_WATER_H,
                                                  G_WATER_H, G_WATER_H, G_BONDED,  G_BONDED,  G_BONDED,  G_UNBONDED};
inline constexpr int MAX_SYMBOLS_PER_UNIT[NUM_STREAMS] = {1, 1, 3, 1, 2, 1, 1, 2, 1, 2, 1, 3};

// Input that BIOMD does not code; SZ_compress gives it to LORENZO_REG with first-order Lorenzo alone, which codes
// coordinates best of its predictors and keeps NaN and Inf exactly.
struct Fallback : std::runtime_error {
    using std::runtime_error::runtime_error;
    static void apply(Config &conf) {
        conf.cmprAlgo = ALGO_LORENZO_REG;
        conf.lorenzo = true, conf.lorenzo2 = false, conf.regression = false;
    }
};

inline int64_t round_half_away(double y) { return int64_t(y + std::copysign(0.5, y)); }
inline int64_t isqrt_round(int64_t n) { return n <= 0 ? 0 : int64_t(std::sqrt(double(n)) + 0.5); }
// d = atom i - atom j
inline void displacement(const int32_t *q, size_t i, size_t j, int64_t d[3]) {
    for (int c = 0; c < 3; c++) d[c] = int64_t(q[3 * i + c]) - q[3 * j + c];
}

// The largest axis of p, dropped by the sphere code, and the other two, kept.
ALWAYS_INLINE void split_axes(const int64_t p[3], int &drop, int &keep0, int &keep1) {
    drop = std::llabs(p[1]) > std::llabs(p[0]) ? 1 : 0;
    if (std::llabs(p[2]) > std::llabs(p[drop])) drop = 2;
    keep0 = (drop + 1) % 3;
    keep1 = (drop + 2) % 3;
}

// ------------------------------------------------------------------------------------------------ layout
enum : uint8_t { K_WATER_O, K_WATER_H, K_OTHER };
struct Layout {
    std::vector<uint8_t> kind;  // per atom: K_WATER_O (then its H at i + 1, i + 2), K_WATER_H, K_OTHER
    std::vector<uint16_t>
        bond;  // other atoms: 0 = none, else offset * 16 + class + 1: the class sphere around i - offset
    double water_oh = 0, water_hh = 0;  // O-H and H-H distances of the water (nm)
    std::vector<double> bond_lengths;   // bond-length classes (nm)
};
inline size_t bond_offset(uint16_t bond) { return bond >> 4; }
inline size_t bond_class(uint16_t bond) { return (bond & 15) - 1; }

// Rigid water on frame x: kind, water_oh, water_hh. The tolerances grow with the lattice step, so rounded input (xtc
// files, or data this codec decompressed) still fits.
template <class T>
void detect_water(const T *x, size_t atoms, Layout &layout, double step) {
    layout.kind.assign(atoms, K_OTHER);
    layout.water_oh = 0;
    auto dist2 = [x](size_t a, size_t b) {
        float dx = float(x[3 * a] - x[3 * b]), dy = float(x[3 * a + 1] - x[3 * b + 1]),
              dz = float(x[3 * a + 2] - x[3 * b + 2]);
        return dx * dx + dy * dy + dz * dz;
    };
    // Candidates: O-H in [0.08, 0.125] nm, H-H in [0.13, 0.2] nm, probed at PROBES places across the frame. Rigid
    // water gives sharp O-H / H-H peaks, flexible CH2/NH2 groups broad ones.
    constexpr float OH2_MIN = 0.0064f, OH2_MAX = 0.015625f, HH2_MIN = 0.0169f, HH2_MAX = 0.04f;  // nm^2
    constexpr size_t PROBES = 300, MIN_COUNT = 16;  // fewer candidates, or waters, than MIN_COUNT: no water
    std::vector<float> cand_oh, cand_hh;
    const size_t probe = std::max<size_t>(1, atoms / PROBES);
    for (size_t i = 0; i + 2 < atoms; i += probe)
        for (size_t o = 0; o < 3 && i + o + 2 < atoms; o++) {
            const size_t k = i + o;
            const float a = dist2(k, k + 1), b = dist2(k, k + 2), c = dist2(k + 1, k + 2);
            if (a > OH2_MIN && a < OH2_MAX && b > OH2_MIN && b < OH2_MAX && c > HH2_MIN && c < HH2_MAX) {
                cand_hh.push_back(std::sqrt(c));
                cand_oh.push_back(std::sqrt(a));
                cand_oh.push_back(std::sqrt(b));
                break;
            }
        }
    if (cand_hh.size() < MIN_COUNT) return;
    // tolerances (nm) stay well below the O-H / C-H difference (> 0.01 nm) that separates water from CH2/NH2
    const double tight = std::max(0.001, step), tol = std::max(0.002, 2.0 * step), loose = std::max(0.01, 3.0 * tight);
    auto peak = [](std::vector<float> v) {  // median of the values within PEAK_SPREAD of the heaviest PEAK_WINDOW
        constexpr float PEAK_WINDOW = 0.004f, PEAK_SPREAD = 0.003f;  // nm
        std::sort(v.begin(), v.end());
        size_t best = 0, n = 0;
        for (size_t lo = 0, hi = 0; hi < v.size(); hi++) {
            while (v[hi] - v[lo] > PEAK_WINDOW) lo++;
            if (hi - lo + 1 > n) n = hi - lo + 1, best = (lo + hi) / 2;
        }
        std::vector<float> w;
        for (float y : v)
            if (std::fabs(y - v[best]) < PEAK_SPREAD) w.push_back(y);
        return double(w[w.size() / 2]);
    };
    const double hh = peak(cand_hh);
    std::vector<float> oh_of_hh;
    for (size_t c = 0; c < cand_hh.size(); c++)
        if (std::fabs(cand_hh[c] - hh) < tol)
            oh_of_hh.push_back(cand_oh[2 * c]), oh_of_hh.push_back(cand_oh[2 * c + 1]);
    if (oh_of_hh.empty()) return;
    const double oh = peak(oh_of_hh);
    // rigidity test: rigid water puts at least six in ten of the nearby candidates within `tight` of both peaks
    size_t n_tight = 0, n_loose = 0, n_rigid = 0;
    double sum_oh = 0, sum_hh = 0;
    for (size_t c = 0; c < cand_hh.size(); c++) {
        const double dh = std::fabs(cand_hh[c] - hh), da = std::fabs(cand_oh[2 * c] - oh),
                     db = std::fabs(cand_oh[2 * c + 1] - oh);
        if (dh < loose && da < loose && db < loose) {
            n_loose++;
            n_tight += dh < tight && da < tight && db < tight;
        }
        if (dh < tol && da < tol && db < tol)
            sum_oh += double(cand_oh[2 * c]) + cand_oh[2 * c + 1], sum_hh += cand_hh[c], n_rigid++;
    }
    if (n_loose < MIN_COUNT || n_tight * 10 < n_loose * 6) return;
    const float oh2_lo = float((oh - tol) * (oh - tol)), oh2_hi = float((oh + tol) * (oh + tol));
    const float hh2_lo = float((hh - tol) * (hh - tol)), hh2_hi = float((hh + tol) * (hh + tol));
    size_t waters = 0;
    for (size_t i = 0; i + 2 < atoms;) {
        float a = dist2(i, i + 1), b, c;
        if (a > oh2_lo && a < oh2_hi && (b = dist2(i, i + 2)) > oh2_lo && b < oh2_hi &&
            (c = dist2(i + 1, i + 2)) > hh2_lo && c < hh2_hi) {
            layout.kind[i] = K_WATER_O;
            layout.kind[i + 1] = layout.kind[i + 2] = K_WATER_H;
            waters++;
            i += 3;
        } else {
            i++;
        }
    }
    if (waters < MIN_COUNT) {
        std::fill(layout.kind.begin(), layout.kind.end(), K_OTHER);
        return;
    }
    layout.water_oh = sum_oh / (2.0 * n_rigid);
    layout.water_hh = sum_hh / double(n_rigid);
}

// Bonds of the other atoms on frame x: each takes the nearest of its previous MAX_BOND_OFFSET atoms if that is within
// bond range, and the bond-length classes are the peaks of those distances (1e-4 nm bins over 0.01 .. 0.25 nm; the
// virtual site of 4-site water is a 0.015 nm bond to its O).
template <class T>
void detect_bonds(const T *x, size_t atoms, Layout &layout) {
    // bond lengths in bins of BIN_WIDTH from BOND_MIN to 0.25 nm; a nearest previous atom outside that range is no bond
    constexpr double BOND_MIN = 0.01, BIN_WIDTH = 1e-4, BINS_PER_NM = 1e4;
    constexpr size_t NUM_BINS = 2400;
    constexpr float BOND2_MIN = 0.0001f, BOND2_MAX = 0.0625f;  // BOND_MIN^2, 0.25^2
    // a class: the heaviest window of +-PEAK_HALF bins, of at least MIN_PEAK atoms and a MIN_SHARE-th of the bonded
    // ones; it then suppresses +-SUPPRESS_HALF bins, and takes the bins within ASSIGN_HALF that are nearest to it
    constexpr size_t PEAK_HALF = 10, SUPPRESS_HALF = 30;
    constexpr long ASSIGN_HALF = 40;
    constexpr int64_t MIN_PEAK = 8, MIN_SHARE = 200;
    layout.bond.assign(atoms, 0);
    layout.bond_lengths.clear();
    std::vector<uint8_t> partner(atoms, 0);  // offset of the nearest previous atom
    std::vector<uint16_t> bin(atoms, uint16_t(NUM_BINS));
    std::vector<int64_t> hist(NUM_BINS, 0);
    for (size_t i = 1; i < atoms; i++) {
        if (layout.kind[i] != K_OTHER) continue;
        float best = 1e30f;
        unsigned offset = 0;  // of the nearest previous atom
        for (unsigned o = 1; o <= unsigned(MAX_BOND_OFFSET) && o <= i; o++) {
            const float dx = float(x[3 * i] - x[3 * (i - o)]), dy = float(x[3 * i + 1] - x[3 * (i - o) + 1]),
                        dz = float(x[3 * i + 2] - x[3 * (i - o) + 2]), v = dx * dx + dy * dy + dz * dz;
            offset = v < best ? o : offset;
            best = std::min(best, v);
        }
        partner[i] = uint8_t(offset);
        if (best > BOND2_MIN && best < BOND2_MAX) {
            bin[i] = uint16_t(std::min(NUM_BINS - 1, size_t((std::sqrt(best) - float(BOND_MIN)) * float(BINS_PER_NM))));
            hist[bin[i]]++;
        }
    }
    std::vector<int8_t> bin_class(NUM_BINS + 1, -1);
    std::vector<uint8_t> bin_distance(NUM_BINS, 255);
    std::vector<int64_t> left = hist, prefix(NUM_BINS + 1);
    int64_t total = 0;
    for (int64_t h : hist) total += h;
    for (int cls = 0; cls < MAX_BOND_CLASSES; cls++) {
        for (size_t k = 0; k < NUM_BINS; k++) prefix[k + 1] = prefix[k] + left[k];
        size_t peak = 0;
        int64_t peak_weight = 0;
        for (size_t k = 0; k < NUM_BINS; k++) {
            const int64_t w =
                prefix[std::min(NUM_BINS, k + PEAK_HALF + 1)] - prefix[k >= PEAK_HALF ? k - PEAK_HALF : 0];
            if (w > peak_weight) peak_weight = w, peak = k;
        }
        if (peak_weight < MIN_PEAK || peak_weight * MIN_SHARE < total) break;
        double weight = 0, length = 0;
        for (size_t w = peak >= PEAK_HALF ? peak - PEAK_HALF : 0; w < std::min(NUM_BINS, peak + PEAK_HALF + 1); w++)
            weight += hist[w], length += hist[w] * (BOND_MIN + (w + 0.5) * BIN_WIDTH);
        layout.bond_lengths.push_back(length / weight);
        for (size_t w = peak >= SUPPRESS_HALF ? peak - SUPPRESS_HALF : 0;
             w < std::min(NUM_BINS, peak + SUPPRESS_HALF + 1); w++)
            left[w] = 0;
        const long centre = long((layout.bond_lengths.back() - BOND_MIN) * BINS_PER_NM);
        for (long k = std::max(0L, centre - ASSIGN_HALF); k <= std::min(long(NUM_BINS) - 1, centre + ASSIGN_HALF); k++)
            if (uint8_t(std::labs(k - centre)) < bin_distance[k])
                bin_distance[k] = uint8_t(std::labs(k - centre)), bin_class[k] = int8_t(cls);
    }
    for (size_t i = 0; i < atoms; i++)
        if (bin_class[bin[i]] >= 0) layout.bond[i] = uint16_t(partner[i] * 16 + bin_class[bin[i]] + 1);
}

// ------------------------------------------------------------------------------------------------ symbols
// the position in each stream, and raw bits (LSB first)
struct SymbolWriter {
    int *cursor[NUM_STREAMS];
    BitAppender *raw;
    ALWAYS_INLINE void put(int s, uint32_t v) { *cursor[s]++ = int(v); }
    ALWAYS_INLINE void put_bits(uint64_t v, int bits) { raw->put(v, bits); }  // bits <= 56
};
// bits of a sample, for the predictor choice
struct CostEstimate {
    double bits = 0;
    inline void put(int, uint32_t v) { bits += std::log2(double(v) + 1.0) + 1; }
    inline void put_bits(uint64_t, int n) { bits += n; }
};
// a cursor per stream, and one on the raw bits
struct SymbolReader {
    const int *cursor[NUM_STREAMS], *end[NUM_STREAMS];
    BitConsumer raw{nullptr, nullptr};
    ALWAYS_INLINE uint32_t get(int s) {
        if (cursor[s] == end[s]) throw std::runtime_error("SZ3 BioMD: corrupt stream");
        return uint32_t(*cursor[s]++);
    }
    ALWAYS_INLINE uint64_t get_bits(int bits) { return raw.get(bits); }  // bits <= 56
};

// v as one symbol below ESCAPE, else as ESCAPE + v / 256 and its low byte as raw bits
template <class Sink>
ALWAYS_INLINE void put_escaped(Sink &out, int s, uint32_t v) {
    if (v < ESCAPE) return out.put(s, v);
    out.put(s, ESCAPE + (v >> 8));
    out.put_bits(v & 255, 8);
}
ALWAYS_INLINE uint32_t get_escaped(SymbolReader &in, int s) {
    const uint32_t v = in.get(s);
    return v < ESCAPE ? v : (v - ESCAPE) << 8 | uint32_t(in.get_bits(8));
}

struct FrameContext {
    int32_t *cur, *prev;  // this frame and the one before, on the lattice
    const Layout *layout;
    int64_t water_oh2, water_hh2;  // squared O-H and H-H distances of the water (lattice units)
    const int64_t *bond_r2;        // squared bond lengths per class (lattice units)
    // the box of the water O's: corner, sides, and the bits of a point in it (x + Rx (y + Ry z)), or -1 if that
    // takes more than 56: then each coordinate in the bits of its side
    int32_t box_min[3];
    uint64_t box_size[3];
    int box_bits, side_bits[3];
};
inline void set_water_box(FrameContext &frame, const uint32_t *size) {  // sides of at most 2^30
    uint64_t *R = frame.box_size;
    for (int c = 0; c < 3; c++) R[c] = size[c];
    // bits of the largest value below n
    auto bits_below = [](uint64_t n) { return int(vector_bit_width(std::vector<uint64_t>{n - 1})); };
    frame.box_bits = R[0] * R[1] < (uint64_t(1) << 56) / R[2] ? bits_below(R[0] * R[1] * R[2]) : -1;
    for (int c = 0; c < 3; c++) frame.side_bits[c] = bits_below(R[c]);
}

// a water O of an intra frame: its point in the box
template <class Sink>
inline void put_box_point(const FrameContext &frame, const int32_t *p, Sink &out) {
    const uint64_t *R = frame.box_size;
    uint64_t u[3];
    for (int c = 0; c < 3; c++) u[c] = uint64_t(int64_t(p[c]) - frame.box_min[c]);
    if (frame.box_bits >= 0) return out.put_bits(u[0] + R[0] * (u[1] + R[1] * u[2]), frame.box_bits);
    for (int c = 0; c < 3; c++) out.put_bits(u[c], frame.side_bits[c]);
}
inline void get_box_point(const FrameContext &frame, int32_t *p, SymbolReader &in) {
    const uint64_t *R = frame.box_size;
    uint64_t v = frame.box_bits >= 0 ? in.get_bits(frame.box_bits) : 0;
    for (int c = 0; c < 3; c++) {
        const uint64_t u = frame.box_bits >= 0 ? v % R[c] : in.get_bits(frame.side_bits[c]);
        v /= R[c];
        p[c] = int32_t(frame.box_min[c] + int64_t(u));
    }
}

// coordinate c of unbonded atom i: from the atom before (intra) or from the previous frame
ALWAYS_INLINE int64_t predict_unbonded(const FrameContext &frame, int mode, size_t i, int c) {
    return mode ? frame.prev[3 * i + c] : (i ? frame.cur[3 * (i - 1) + c] : 0);
}

// Atom i on the sphere |d|^2 = r2 around atom j: the largest axis of d is dropped and rebuilt from the other two, with
// a radial residual. Intra (mode 0) the face (2 axis + sign) is sent with the residual and the two kept coordinates as
// they are; else the vector p of the previous frame fixes the axis and the sign, and the kept coordinates are residuals
// against p.
template <class Sink>
ALWAYS_INLINE void put_sphere(const FrameContext &frame, int mode, size_t i, size_t j, int64_t r2, int s_face,
                              int s_kept, int s_radial, Sink &out) {
    int64_t d[3], p[3] = {0, 0, 0};
    int drop, k0, k1;
    displacement(frame.cur, i, j, d);
    if (mode) displacement(frame.prev, i, j, p);
    split_axes(mode ? p : d, drop, k0, k1);
    const bool neg = (mode ? p[drop] : d[drop]) < 0;
    const uint32_t radial = zigzag((neg ? -d[drop] : d[drop]) - isqrt_round(r2 - d[k0] * d[k0] - d[k1] * d[k1]));
    if (mode == 0) out.put(s_face, uint32_t(2 * drop + neg) + 6 * std::min(radial, 15u));
    if (mode != 0 || radial >= 15) out.put(s_radial, mode == 0 ? radial - 15 : radial);
    out.put(s_kept, zigzag(d[k0] - p[k0]));
    out.put(s_kept, zigzag(d[k1] - p[k1]));
}
ALWAYS_INLINE void get_sphere(const FrameContext &frame, int mode, size_t i, size_t j, int64_t r2, int s_face,
                              int s_kept, int s_radial, SymbolReader &in) {
    int64_t d[3], p[3] = {0, 0, 0};
    int drop, k0, k1;
    bool neg;
    uint32_t radial;
    if (mode == 0) {
        const uint32_t v = in.get(s_face);
        drop = int(v % 6 / 2);
        k0 = (drop + 1) % 3;
        k1 = (drop + 2) % 3;
        neg = v & 1;
        radial = v / 6;
        if (radial == 15) radial += in.get(s_radial);
    } else {
        displacement(frame.prev, i, j, p);
        split_axes(p, drop, k0, k1);
        neg = p[drop] < 0;
        radial = in.get(s_radial);
    }
    d[k0] = p[k0] + unzigzag(in.get(s_kept));
    d[k1] = p[k1] + unzigzag(in.get(s_kept));
    // unsigned: a corrupt stream's squares may pass 2^63, a valid one's stay below 2^60
    const uint64_t rest = uint64_t(r2) - uint64_t(d[k0]) * uint64_t(d[k0]) - uint64_t(d[k1]) * uint64_t(d[k1]);
    const int64_t r = isqrt_round(int64_t(rest)) + unzigzag(radial);
    d[drop] = neg ? -r : r;
    for (int c = 0; c < 3; c++) frame.cur[3 * i + c] = int32_t(frame.cur[3 * j + c] + d[c]);
}

// H2 - O = d on the circle |d|^2 = oh2, 2 d.h1 = oh2 + |h1|^2 - hh2, h1 = H1 - O. The coordinate of d on the axis
// where h1 is smallest leaves two points of the circle, one on each side of the plane through h1 and that axis; the
// side of d (intra: sent, else as a change from the side of the prediction p) picks one, and two residuals correct it.
// Integer products and a square root and a division of exact doubles, so the decoder finds the same point. Beyond the
// lattice sizes of a water (a molecule split by the boundary) the point is {coord on the axis, 0, 0}.
ALWAYS_INLINE int circle_axis(const int64_t h1[3]) {
    int axis = std::llabs(h1[1]) < std::llabs(h1[0]) ? 1 : 0;
    return std::llabs(h1[2]) < std::llabs(h1[axis]) ? 2 : axis;
}
// the side of v: 1 when (-h1[j], h1[i]) . (v[i], v[j]) < 0, for the other two axes i, j
ALWAYS_INLINE int circle_side(const int64_t h1[3], int axis, const int64_t v[3]) {
    const int i = (axis + 1) % 3, j = (axis + 2) % 3;
    return h1[i] * v[j] - h1[j] * v[i] < 0;
}
ALWAYS_INLINE void circle_point(const int64_t h1[3], int64_t oh2, int64_t hh2, int axis, int64_t coord, int side,
                                int64_t point[3]) {
    const int i = (axis + 1) % 3, j = (axis + 2) % 3;
    const int64_t u = h1[i], w = h1[j], lim = 16384;
    const bool ok = std::llabs(u) < lim && std::llabs(w) < lim && std::llabs(h1[axis]) < lim &&
                    std::llabs(coord) < lim && u * u + w * w > 0;
    const int64_t g2 = ok ? u * u + w * w : 1;
    const int64_t M = ok ? oh2 + g2 + h1[axis] * h1[axis] - hh2 - 2 * h1[axis] * coord : 0;  // 2 (u P_i + w P_j)
    const int64_t s = ok ? (side ? -1 : 1) * isqrt_round(4 * g2 * (oh2 - coord * coord) - M * M) : 0;
    point[axis] = coord;
    point[i] = ok ? round_half_away(double(M * u - s * w) / double(2 * g2)) : 0;
    point[j] = ok ? round_half_away(double(M * w + s * u) / double(2 * g2)) : 0;
}
// H2 - O of the previous frame, or 0 beyond the lattice sizes of a water, so that 2 zigzag(d - p) + 1 fits 32 bits
ALWAYS_INLINE void previous_h2(const FrameContext &frame, size_t o, int64_t p[3]) {
    displacement(frame.prev, o + 2, o, p);
    if (std::llabs(p[0]) >= 16384 || std::llabs(p[1]) >= 16384 || std::llabs(p[2]) >= 16384) p[0] = p[1] = p[2] = 0;
}
template <class Sink>
ALWAYS_INLINE void put_circle(const FrameContext &frame, int mode, size_t o, Sink &out) {
    int64_t h1[3], d[3], p[3], point[3];
    displacement(frame.cur, o + 1, o, h1);
    displacement(frame.cur, o + 2, o, d);
    const int axis = circle_axis(h1), side = circle_side(h1, axis, d);
    if (mode) {
        previous_h2(frame, o, p);
        out.put(S_WATER_H2_AXIS, zigzag(d[axis] - p[axis]) * 2 + uint32_t(side != circle_side(h1, axis, p)));
    } else {
        out.put(S_WATER_H2_AXIS, zigzag(d[axis]) * 2 + uint32_t(side));
    }
    circle_point(h1, frame.water_oh2, frame.water_hh2, axis, d[axis], side, point);
    for (int c = 0; c < 3; c++)
        if (c != axis) out.put(S_WATER_H2_RESIDUAL, zigzag(d[c] - point[c]));
}
ALWAYS_INLINE void get_circle(const FrameContext &frame, int mode, size_t o, SymbolReader &in) {
    int64_t h1[3], p[3] = {0, 0, 0}, point[3];
    displacement(frame.cur, o + 1, o, h1);
    const int axis = circle_axis(h1);
    if (mode) previous_h2(frame, o, p);
    const uint32_t v = in.get(S_WATER_H2_AXIS);
    const int side = int(v & 1) ^ (mode ? circle_side(h1, axis, p) : 0);
    circle_point(h1, frame.water_oh2, frame.water_hh2, axis, p[axis] + unzigzag(v >> 1), side, point);
    for (int c = 0; c < 3; c++)
        frame.cur[3 * (o + 2) + c] =
            int32_t(frame.cur[3 * o + c] + (c == axis ? point[c] : point[c] + unzigzag(in.get(S_WATER_H2_RESIDUAL))));
}

// A water O: from the previous frame, or (intra) its point in the box.
template <class Sink>
ALWAYS_INLINE void put_water_o(const FrameContext &frame, int mode, size_t i, Sink &out) {
    if (mode)
        for (int c = 0; c < 3; c++) out.put(S_WATER_O, zigzag(frame.cur[3 * i + c] - frame.prev[3 * i + c]));
    else
        put_box_point(frame, frame.cur + 3 * i, out);
}

// One atom: a water (O and the rest of its molecule), or another atom. Water H come with their O.
template <class Sink>
ALWAYS_INLINE void put_atom(const FrameContext &frame, const int *mode, size_t i, Sink &out) {
    const Layout &layout = *frame.layout;
    const int32_t *q = frame.cur;
    if (layout.kind[i] == K_WATER_O) {
        put_water_o(frame, mode[G_WATER_O], i, out);
        put_sphere(frame, mode[G_WATER_H], i + 1, i, frame.water_oh2, S_WATER_H1_FACE, S_WATER_H1_KEPT,
                   S_WATER_H1_RADIAL, out);
        put_circle(frame, mode[G_WATER_H], i, out);
    } else if (layout.kind[i] == K_OTHER) {
        const uint16_t bond = layout.bond[i];
        if (bond)
            put_sphere(frame, mode[G_BONDED], i, i - bond_offset(bond), frame.bond_r2[bond_class(bond)], S_BOND_FACE,
                       S_BOND_KEPT, S_BOND_RADIAL, out);
        else
            for (int c = 0; c < 3; c++)
                put_escaped(out, S_UNBONDED, zigzag(q[3 * i + c] - predict_unbonded(frame, mode[G_UNBONDED], i, c)));
    }
}
ALWAYS_INLINE void get_atom(const FrameContext &frame, const int *mode, size_t i, SymbolReader &in) {
    const Layout &layout = *frame.layout;
    int32_t *q = frame.cur;
    if (layout.kind[i] == K_WATER_O) {
        if (mode[G_WATER_O])
            for (int c = 0; c < 3; c++) q[3 * i + c] = int32_t(frame.prev[3 * i + c] + unzigzag(in.get(S_WATER_O)));
        else
            get_box_point(frame, q + 3 * i, in);
        get_sphere(frame, mode[G_WATER_H], i + 1, i, frame.water_oh2, S_WATER_H1_FACE, S_WATER_H1_KEPT,
                   S_WATER_H1_RADIAL, in);
        get_circle(frame, mode[G_WATER_H], i, in);
    } else if (layout.kind[i] == K_OTHER) {
        const uint16_t bond = layout.bond[i];
        if (bond)
            get_sphere(frame, mode[G_BONDED], i, i - bond_offset(bond), frame.bond_r2[bond_class(bond)], S_BOND_FACE,
                       S_BOND_KEPT, S_BOND_RADIAL, in);
        else
            for (int c = 0; c < 3; c++)
                q[3 * i + c] =
                    int32_t(predict_unbonded(frame, mode[G_UNBONDED], i, c) + unzigzag(get_escaped(in, S_UNBONDED)));
    }
}

// q = round(x * inv_step), inv_step the reciprocal of the step
template <class T>
inline void quantize(const T *x, size_t n, double inv_step, int32_t *q) {
    for (size_t i = 0; i < n; i++) {
        const double y = double(x[i]) * inv_step;
        q[i] = int32_t(y + std::copysign(0.5, y));
    }
}

// the box of the water O's of frame q
inline void water_box(const int32_t *q, const std::vector<uint32_t> &waters, int32_t lo[3], int32_t hi[3]) {
    for (int c = 0; c < 3; c++) lo[c] = hi[c] = waters.empty() ? 0 : q[3 * waters[0] + c];
    for (uint32_t i : waters)
        for (int c = 0; c < 3; c++) {
            lo[c] = std::min(lo[c], q[3 * i + c]);
            hi[c] = std::max(hi[c], q[3 * i + c]);
        }
}

// Per group, the mode (0 intra, 1 previous frame) that costs fewer bits on a sample of about 256 of its atoms.
inline void choose_modes(const FrameContext &frame, const std::vector<uint32_t> *const group_atoms[NUM_GROUPS],
                         int mode[NUM_GROUPS]) {
    for (int g = 0; g < NUM_GROUPS; g++) {
        const auto &sample = *group_atoms[g];
        double best = 1e300;
        for (int m = 0; m < 2 && !sample.empty(); m++) {
            CostEstimate cost;
            int trial[NUM_GROUPS] = {mode[0], mode[1], mode[2], mode[3]};
            trial[g] = m;
            const size_t stride = std::max<size_t>(1, sample.size() / 256);
            for (size_t k = 0; k < sample.size(); k += stride) {
                if (g == G_WATER_O)  // the O alone: put_atom would add its molecule's H
                    put_water_o(frame, m, sample[k], cost);
                else
                    put_atom(frame, trial, sample[k], cost);
            }
            if (cost.bits < best) best = cost.bits, mode[g] = m;
        }
    }
}

}  // namespace biomd

template <class T, uint N>
class SZBioMDDecomposition : public concepts::DecompositionInterface<T, int, N> {
   public:
    explicit SZBioMDDecomposition(const Config &conf)
        : frames_(N == 3 ? conf.dims[0] : 1), atoms_(N >= 2 ? conf.dims[N - 2] : 1), error_bound_(conf.absErrorBound) {
        if (N > 3 || conf.dims[N - 1] != 3) throw std::invalid_argument("SZ3 BioMD: data must be {frames, atoms, 3}");
        if (!std::is_floating_point<T>::value) throw std::invalid_argument("SZ3 BioMD: data must be float or double");
    }

    // The symbol streams of the chunk, one after the other as [count, symbols]; the rest goes to save().
    std::vector<int> compress(const Config & /*conf*/, T *data) override {
        using namespace biomd;
        const size_t frame_values = atoms_ * 3;
        find_fill_frames(data);
        choose_step(data);
        Layout layout;
        std::vector<uint32_t> waters, bonded, unbonded;  // water O's, other atoms
        detect_layout(data, layout, waters, bonded, unbonded);

        // a buffer per stream, as large as the symbols the chunk can put in it (every value is written before it is
        // read)
        const std::vector<uint32_t> *group_atoms[NUM_GROUPS] = {&waters, &waters, &bonded, &unbonded};
        std::unique_ptr<int[]> buffers[NUM_STREAMS];
        SymbolWriter writer;
        for (int s = 0; s < NUM_STREAMS; s++) {
            const size_t n = STREAM_GROUP[s] < 0 ? atoms_ : coded_frames_ * group_atoms[STREAM_GROUP[s]]->size();
            buffers[s].reset(new int[MAX_SYMBOLS_PER_UNIT[s] * n + 1]);
            writer.cursor[s] = buffers[s].get();
        }
        raw_bits_.clear();
        BitAppender raw(raw_bits_);
        writer.raw = &raw;
        for (size_t k = 0; k < waters.size(); k++) writer.put(S_WATER_GAP, waters[k] - (k ? waters[k - 1] : 0));
        for (size_t i = 0; i < atoms_; i++)
            if (layout.kind[i] == K_OTHER) writer.put(S_BOND_REF, layout.bond[i]);

        std::unique_ptr<int32_t[]> lattice(new int32_t[2 * frame_values]);  // written before read
        FrameContext frame = frame_context(lattice.get(), layout);
        const int intra[NUM_GROUPS] = {0, 0, 0, 0};
        std::fill(modes_, modes_ + NUM_GROUPS, 0);
        box_min_.clear();
        box_size_.clear();
        for (size_t t = 0; t < coded_frames_; t++) {
            quantize(data + t * frame_values, frame_values, 1.0 / step_, frame.cur);
            int32_t box_max[3];
            uint32_t box_size[3];
            water_box(frame.cur, waters, frame.box_min, box_max);
            for (int c = 0; c < 3; c++) box_size[c] = uint32_t(box_max[c] - frame.box_min[c]) + 1;
            set_water_box(frame, box_size);
            // per group, the predictor that is cheapest on a sample of frame 1, for every frame after the first
            if (t == 1) choose_modes(frame, group_atoms, modes_);
            if (needs_box(t)) {
                box_min_.insert(box_min_.end(), frame.box_min, frame.box_min + 3);
                box_size_.insert(box_size_.end(), box_size, box_size + 3);
            }
            for (size_t i = 0; i < atoms_; i++) put_atom(frame, t ? modes_ : intra, i, writer);
            std::swap(frame.cur, frame.prev);
        }
        raw.flush();

        std::vector<int> bins;
        size_t n = NUM_STREAMS;
        for (int s = 0; s < NUM_STREAMS; s++) {
            if (size_t(writer.cursor[s] - buffers[s].get()) > size_t(std::numeric_limits<int>::max()))
                throw Fallback("SZ3 BioMD: a stream of more symbols than an int counts");
            n += size_t(writer.cursor[s] - buffers[s].get());
        }
        bins.reserve(n);
        for (int s = 0; s < NUM_STREAMS; s++) {
            bins.push_back(int(writer.cursor[s] - buffers[s].get()));
            bins.insert(bins.end(), buffers[s].get(), writer.cursor[s]);
        }
        return bins;
    }

    // The chunk from its symbol streams, each as [count, symbols], after load().
    T *decompress(const Config & /*conf*/, std::vector<int> &quant_inds, T *dec_data) override {
        using namespace biomd;
        const size_t frame_values = atoms_ * 3;
        SymbolReader reader;
        for (size_t k = 0, s = 0; s < NUM_STREAMS; s++) {
            const size_t n = k < quant_inds.size() ? size_t(quant_inds[k]) : 0;
            if (k >= quant_inds.size() || n > quant_inds.size() - k - 1)
                throw std::runtime_error("SZ3 BioMD: corrupt stream");
            reader.cursor[s] = quant_inds.data() + k + 1;
            reader.end[s] = reader.cursor[s] + n;
            k += n + 1;
        }
        reader.raw = BitConsumer(raw_bits_.data(), raw_bits_.data() + raw_bits_.size());
        Layout layout;
        read_layout(reader, layout);
        std::unique_ptr<int32_t[]> lattice(new int32_t[2 * frame_values]);  // written before read
        FrameContext frame = frame_context(lattice.get(), layout);
        const int intra[NUM_GROUPS] = {0, 0, 0, 0};
        for (size_t t = 0, box = 0; t < coded_frames_; t++) {
            if (needs_box(t)) {
                for (int c = 0; c < 3; c++) frame.box_min[c] = box_min_[box * 3 + c];
                set_water_box(frame, &box_size_[box++ * 3]);
            }
            for (size_t i = 0; i < atoms_; i++) get_atom(frame, t ? modes_ : intra, i, reader);
            T *x = dec_data + t * frame_values;
            for (size_t i = 0; i < frame_values; i++) x[i] = T(double(frame.cur[i]) * step_);
            std::swap(frame.cur, frame.prev);
        }
        // a stream that holds more than the frames took is corrupt
        for (int s = 0; s < NUM_STREAMS; s++)
            if (reader.cursor[s] != reader.end[s]) throw std::runtime_error("SZ3 BioMD: corrupt stream");
        if (!reader.raw.at_end()) throw std::runtime_error("SZ3 BioMD: corrupt stream");
        std::fill(dec_data + coded_frames_ * frame_values, dec_data + frames_ * frame_values, fill_value_);
        return dec_data;
    }

    void save(uchar *&c) override {
        write(uint64_t(coded_frames_), c);
        write(fill_value_, c);
        write(step_, c);
        write(water_oh2_, c);
        write(water_hh2_, c);
        write(uint64_t(num_waters_), c);
        write(uint8_t(bond_r2_.size()), c);
        if (!bond_r2_.empty()) write(bond_r2_.data(), bond_r2_.size(), c);
        for (int m : modes_) write(uint8_t(m), c);
        if (!box_min_.empty()) {  // no box without water
            write(box_min_.data(), box_min_.size(), c);
            write(box_size_.data(), box_size_.size(), c);
        }
        write(uint64_t(raw_bits_.size()), c);
        if (!raw_bits_.empty()) write(raw_bits_.data(), raw_bits_.size(), c);
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        uint8_t classes;
        uint64_t waters, coded, raw_size = 0;
        read(coded, c, remaining_length);
        if (coded > frames_) throw std::runtime_error("SZ3 BioMD: corrupt stream");
        coded_frames_ = size_t(coded);
        read(fill_value_, c, remaining_length);
        read(step_, c, remaining_length);
        read(water_oh2_, c, remaining_length);
        read(water_hh2_, c, remaining_length);
        read(waters, c, remaining_length);
        num_waters_ = size_t(std::min<uint64_t>(waters, atoms_));
        read(classes, c, remaining_length);
        bond_r2_.resize(classes);
        if (classes) read(bond_r2_.data(), classes, c, remaining_length);
        for (int &m : modes_) {
            uint8_t v;
            read(v, c, remaining_length);
            m = v & 1;
        }
        size_t boxes = 0;
        for (size_t t = 0; t < coded_frames_; t++) boxes += needs_box(t);
        box_min_.resize(boxes * 3);
        box_size_.resize(boxes * 3);
        if (boxes) {
            read(box_min_.data(), box_min_.size(), c, remaining_length);
            read(box_size_.data(), box_size_.size(), c, remaining_length);
        }
        read(raw_size, c, remaining_length);
        if (raw_size > remaining_length) throw std::runtime_error("SZ3 BioMD: corrupt stream");
        raw_bits_.resize(size_t(raw_size));
        if (raw_size) read(raw_bits_.data(), raw_bits_.size(), c, remaining_length);
        bool ok = std::isfinite(step_) && step_ > 0 && waters <= atoms_ && water_oh2_ >= 0 &&
                  water_oh2_ < (int64_t(1) << 28) && water_hh2_ >= 0 && water_hh2_ < (int64_t(1) << 30) &&
                  classes <= biomd::MAX_BOND_CLASSES;
        for (int64_t r2 : bond_r2_)
            ok = ok && r2 >= 0;  // a class of MAX_BOND_R2 or more may be stored, but no atom refers to it
        for (uint32_t side : box_size_) ok = ok && side >= 1 && side <= (uint32_t(1) << 30);
        if (!ok) throw std::runtime_error("SZ3 BioMD: corrupt stream");
    }

    // a bound on what save() writes
    size_t size_est() override { return 64 + 8 * bond_r2_.size() + 24 * frames_ + raw_bits_.size(); }

    // every bin, a count or a symbol, is in [0, INT_MAX]
    std::pair<int, int> get_out_range() override { return {0, std::numeric_limits<int>::max()}; }

   private:
    // a frame whose water O are intra stores the box of the water O's
    bool needs_box(size_t t) const { return num_waters_ && (t == 0 || !modes_[biomd::G_WATER_O]); }

    // Trailing frames all of one value (the unwritten rest of a chunk) are stored as that value; a chunk may be all
    // fill (an unwritten chunk of NaN).
    void find_fill_frames(const T *data) {
        const size_t frame_values = atoms_ * 3;
        fill_value_ = data[(frames_ - 1) * frame_values];
        auto is_fill = [&](size_t t) {
            for (size_t i = t * frame_values; i < (t + 1) * frame_values; i++)
                if (memcmp(&data[i], &fill_value_, sizeof(T)) != 0) return false;
            return true;
        };
        for (coded_frames_ = frames_; coded_frames_ > 0 && is_fill(coded_frames_ - 1);) coded_frames_--;
    }

    // Lattice step: 2 (eb - m ulp) for the binade below B, so that q step rounded to T stays within eb for every
    // |q step| < B. The rounding of x * (1 / step) to q, done in double, is off by up to 2 double ulps of x, and q step
    // by 1 more: m = 1 for float, whose ulp is 2^29 times larger, and 4 for double.
    //  * B is the smallest power of two of at least 2^LATTICE_SPAN_BITS eb when m ulp < eb / 2 there and every |q| step
    //    of the chunk stays below B: the step then depends on the bound alone, and a decompressed chunk compressed
    //    again lands on the same lattice points.
    //  * Otherwise B is the power of two above max|x|, with the step eb when m ulp is eb / 2 or more: q eb rounded to T
    //    is then no further from q eb than x is, so within eb of x. Compressed again, the chunk may land on another
    //    lattice, as with any other algorithm.
    //  * A chunk of more than MAX_LATTICE lattice steps, or whose lattice would pass the largest value of T, goes to
    //    LORENZO_REG.
    void choose_step(const T *data) {
        if (!(error_bound_ > 0) || !std::isfinite(error_bound_))
            throw std::invalid_argument("SZ3 BioMD: the error bound must be positive and finite");
        using Bits = typename std::conditional<sizeof(T) == 4, uint32_t, uint64_t>::type;
        Bits max_bits = 0;  // the largest |x| as bits: NaN and inf are above every finite value
        for (size_t i = 0; i < coded_frames_ * atoms_ * 3; i++) {
            Bits b;
            memcpy(&b, &data[i], sizeof(T));
            b &= Bits(~Bits(0)) >> 1;
            max_bits = b > max_bits ? b : max_bits;
        }
        T max_abs;
        memcpy(&max_abs, &max_bits, sizeof(T));
        if (!(max_abs <= std::numeric_limits<T>::max()))
            throw biomd::Fallback("SZ3 BioMD: NaN or Inf in a frame that is not trailing fill");
        const double m = sizeof(T) == 4 ? 1 : 4;
        auto margin_below = [&](int exponent) {  // m ulp in the binade below 2^exponent
            return m * std::max(double(std::numeric_limits<T>::denorm_min()),
                                std::ldexp(1.0, exponent - std::numeric_limits<T>::digits));
        };
        int exponent;
        std::frexp(std::ldexp(error_bound_, biomd::LATTICE_SPAN_BITS), &exponent);
        double margin = margin_below(exponent);
        if (margin < 0.5 * error_bound_) {
            step_ = 2.0 * (error_bound_ - margin);
            const double max_q = std::fabs(double(biomd::round_half_away(double(max_abs) * (1.0 / step_))));
            if (max_q <= biomd::MAX_LATTICE && max_q * step_ < std::ldexp(1.0, exponent) &&
                std::ldexp(1.0, exponent) <= double(std::numeric_limits<T>::max()))
                return;
        }
        std::frexp(double(max_abs), &exponent);
        margin = margin_below(exponent);
        step_ = margin < 0.5 * error_bound_ ? 2.0 * (error_bound_ - margin) : error_bound_;
        // |q step| <= max|x| + step / 2 must stay finite in T
        if (!(step_ > 0) || !std::isfinite(step_) || !std::isfinite(1.0 / step_) ||
            max_abs / step_ > biomd::MAX_LATTICE || !(double(max_abs) + step_ <= double(std::numeric_limits<T>::max())))
            throw biomd::Fallback("SZ3 BioMD: coordinates beyond the lattice the error bound allows");
    }

    // The layout of the first frame, its geometry on the lattice, and the atoms of each group.
    void detect_layout(const T *data, biomd::Layout &layout, std::vector<uint32_t> &waters,
                       std::vector<uint32_t> &bonded, std::vector<uint32_t> &unbonded) {
        using namespace biomd;
        detect_water(data, atoms_, layout, step_);
        const double oh = layout.water_oh / step_, hh = layout.water_hh / step_;
        water_oh2_ = oh < 16384 ? round_half_away(oh * oh) : 0;  // squares that load() accepts
        water_hh2_ = hh < 32768 ? round_half_away(hh * hh) : 0;
        if (!(water_oh2_ > 0 && water_oh2_ < (int64_t(1) << 28) && water_hh2_ > 0 && water_hh2_ < (int64_t(1) << 30))) {
            water_oh2_ = water_hh2_ = 0;  // no water, or too many lattice steps across one for exact products
            std::fill(layout.kind.begin(), layout.kind.end(), K_OTHER);
        }
        detect_bonds(data, atoms_, layout);
        bond_r2_.clear();
        for (double b : layout.bond_lengths) bond_r2_.push_back(round_half_away((b / step_) * (b / step_)));
        for (size_t i = 0; i < atoms_; i++) {
            uint16_t &bond = layout.bond[i];
            if (bond && bond_r2_[bond_class(bond)] >= MAX_BOND_R2) bond = 0;
            if (layout.kind[i] == K_WATER_O) waters.push_back(uint32_t(i));
            if (layout.kind[i] == K_OTHER) (bond ? bonded : unbonded).push_back(uint32_t(i));
        }
        num_waters_ = waters.size();
    }

    // The layout from its streams: the water O's by their gaps, then the bond of each other atom.
    void read_layout(biomd::SymbolReader &reader, biomd::Layout &layout) const {
        using namespace biomd;
        layout.kind.assign(atoms_, K_OTHER);
        for (size_t k = 0, o = 0; k < num_waters_; k++) {
            const uint32_t gap = reader.get(S_WATER_GAP);
            if ((k && gap < 3) || gap > atoms_ || o + gap + 2 >= atoms_)
                throw std::runtime_error("SZ3 BioMD: corrupt stream");
            o += gap;
            layout.kind[o] = K_WATER_O;
            layout.kind[o + 1] = layout.kind[o + 2] = K_WATER_H;
        }
        layout.bond.assign(atoms_, 0);
        for (size_t i = 0; i < atoms_; i++)
            if (layout.kind[i] == K_OTHER) {
                const uint32_t bond = reader.get(S_BOND_REF);
                if (bond > uint32_t(MAX_BOND_OFFSET * 16 + MAX_BOND_CLASSES))
                    throw std::runtime_error("SZ3 BioMD: corrupt stream");
                if (bond && (bond_offset(uint16_t(bond)) == 0 || bond_offset(uint16_t(bond)) > i ||
                             bond_class(uint16_t(bond)) >= bond_r2_.size()))  // a class field of 0 gives SIZE_MAX
                    throw std::runtime_error("SZ3 BioMD: corrupt stream");
                layout.bond[i] = uint16_t(bond);
            }
    }

    // the lattice holds this frame and the one before
    biomd::FrameContext frame_context(int32_t *lattice, const biomd::Layout &layout) const {
        return {lattice, lattice + atoms_ * 3, &layout, water_oh2_, water_hh2_, bond_r2_.data(), {0, 0, 0}, {1, 1, 1},
                0};
    }

    size_t frames_, atoms_, coded_frames_ = 1;  // coded: the frames before the fill
    T fill_value_ = 0;
    double error_bound_, step_ = 0;
    int64_t water_oh2_ = 0, water_hh2_ = 0;  // squared O-H and H-H distances of the water (lattice units)
    size_t num_waters_ = 0;
    std::vector<int64_t> bond_r2_;                 // squared bond lengths per class (lattice units)
    int modes_[biomd::NUM_GROUPS] = {0, 0, 0, 0};  // per group: the predictor of the frames after the first
    std::vector<int32_t> box_min_;    // per frame whose water O are intra: the corner of the box of the water O's
    std::vector<uint32_t> box_size_;  // and its sides
    std::vector<uchar> raw_bits_;     // intra water O's and the low bytes of escaped values
};

template <class T, uint N>
SZBioMDDecomposition<T, N> make_decomposition_biomd(const Config &conf) {
    return SZBioMDDecomposition<T, N>(conf);
}

}  // namespace SZ3
#endif
