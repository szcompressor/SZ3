#ifndef SZ3_BIOMD_SIMD_HPP
#define SZ3_BIOMD_SIMD_HPP

// ALGO_BIOMD, part 2: AVX2 kernels for the intra-frame geometry, picked at run time on x86 (elsewhere, and on CPUs
// without AVX2, the scalar code runs). Every quantity is an exact integer in double (radii are limited so that
// products stay below 2^53) or a product rounded by the same round-to-nearest conversion as the scalar rnd_pred, so
// the results equal the scalar code bit for bit, with or without FMA.

#include "SZ3/compressor/specialized/biomd/BioMDCore.hpp"
#if defined(SZ3_BIOMD_X86_DISPATCH)
#include <immintrin.h>
#endif

namespace SZ3 {
namespace biomd {

constexpr size_t SIMD_PAD = 4;  // vector width, in doubles

// SoA input: d = H1 - O, h = H2 - O, mo = M - O (4-site only); output per molecule, before zigzag
struct WaterBatch {
    std::vector<double> dx, dy, dz, hx, hy, hz, mx, my, mz;
    std::vector<int32_t> face, a, b, e, v, r1, r2, side, m0, m1, m2;
    void resize(size_t n) {
        for (auto *p : {&dx, &dy, &dz, &hx, &hy, &hz, &mx, &my, &mz})
            if (p->size() < n + SIMD_PAD) p->resize(n + SIMD_PAD);
        for (auto *p : {&face, &a, &b, &e, &v, &r1, &r2, &side, &m0, &m1, &m2})
            if (p->size() < n + SIMD_PAD) p->resize(n + SIMD_PAD);
    }
};
// SoA input: d = atom - partner, r2 = squared bond length (lattice units); output as SphereCode
struct SphereBatch {
    std::vector<double> dx, dy, dz, r2;
    std::vector<int32_t> face, a, b, e;
    void resize(size_t n) {
        for (auto *p : {&dx, &dy, &dz, &r2})
            if (p->size() < n + SIMD_PAD) p->resize(n + SIMD_PAD);
        for (auto *p : {&face, &a, &b, &e})
            if (p->size() < n + SIMD_PAD) p->resize(n + SIMD_PAD);
    }
};

#if defined(SZ3_BIOMD_X86_DISPATCH)
namespace k_avx2 {
#define SZ3_BIOMD_TARGET __attribute__((target("avx2")))
#define SZ3_BIOMD_H __attribute__((target("avx2"), always_inline)) static inline
constexpr size_t NL = 4;
typedef __m256d VD;
typedef __m256d VM;
SZ3_BIOMD_H VD set1(double x) { return _mm256_set1_pd(x); }
SZ3_BIOMD_H VD ld(const double *p) { return _mm256_loadu_pd(p); }
SZ3_BIOMD_H VD add(VD a, VD b) { return _mm256_add_pd(a, b); }
SZ3_BIOMD_H VD sub(VD a, VD b) { return _mm256_sub_pd(a, b); }
SZ3_BIOMD_H VD mul(VD a, VD b) { return _mm256_mul_pd(a, b); }
SZ3_BIOMD_H VD dv(VD a, VD b) { return _mm256_div_pd(a, b); }
SZ3_BIOMD_H VD vabs(VD x) { return _mm256_andnot_pd(_mm256_set1_pd(-0.0), x); }
SZ3_BIOMD_H VM gt(VD a, VD b) { return _mm256_cmp_pd(a, b, _CMP_GT_OQ); }
SZ3_BIOMD_H VM lt(VD a, VD b) { return _mm256_cmp_pd(a, b, _CMP_LT_OQ); }
SZ3_BIOMD_H VM le(VD a, VD b) { return _mm256_cmp_pd(a, b, _CMP_LE_OQ); }
SZ3_BIOMD_H VM mor(VM a, VM b) { return _mm256_or_pd(a, b); }
SZ3_BIOMD_H VM mandn(VM a, VM b) { return _mm256_andnot_pd(a, b); }
SZ3_BIOMD_H VD sel(VM m, VD t, VD f) { return _mm256_blendv_pd(f, t, m); }
SZ3_BIOMD_H VD mone(VM m) { return _mm256_and_pd(m, _mm256_set1_pd(1.0)); }
SZ3_BIOMD_H VD rndi(VD y) { return _mm256_cvtepi32_pd(_mm256_cvtpd_epi32(y)); }
SZ3_BIOMD_H VD isq(VD n) {
    const VD zero = _mm256_setzero_pd();
    const VD r = _mm256_round_pd(_mm256_add_pd(_mm256_sqrt_pd(_mm256_max_pd(n, zero)), _mm256_set1_pd(0.5)),
                                 _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
    return sel(le(n, zero), zero, r);
}
SZ3_BIOMD_H void sti(int32_t *o, VD x) { _mm_storeu_si128(reinterpret_cast<__m128i *>(o), _mm256_cvtpd_epi32(x)); }
// intra-frame rigid-water geometry: sphere for H1, circle for H2, virtual site (as enc_WH, mode 0)
SZ3_BIOMD_TARGET static void water_intra(WaterBatch &B, size_t n, int64_t R2, double Rc2R, double invR2, bool four_site,
                                         double vs_a) {
    const VD zero = set1(0.0), one = set1(1.0);
    const VD vR2 = set1(double(R2)), v4R2 = set1(4.0 * double(R2));
    const VD vRc2R = set1(Rc2R), vinvR2 = set1(invR2), vvs = set1(vs_a);
    for (size_t w = 0; w < n; w += NL) {
        const VD d0 = ld(&B.dx[w]), d1 = ld(&B.dy[w]), d2 = ld(&B.dz[w]);
        const VD h0 = ld(&B.hx[w]), h1 = ld(&B.hy[w]), h2 = ld(&B.hz[w]);
        const VD a0 = vabs(d0), a1 = vabs(d1), a2 = vabs(d2);
        // --- H1 sphere: face f = argmax |d| (first on ties), kept coordinates (i, j) = (f+1, f+2) cyclic
        const VM c1 = gt(a1, a0);
        VD m = sel(c1, a1, a0);
        const VM c2 = gt(a2, m);
        m = sel(c2, a2, m);
        const VM f1 = mandn(c2, c1), f2 = c2;  // f == 1, f == 2
        const VD df = sel(f2, d2, sel(f1, d1, d0));
        const VD sa = sel(f2, d0, sel(f1, d2, d1)), sb = sel(f2, d1, sel(f1, d0, d2));
        const VD fnum = add(sel(f2, set1(4.0), sel(f1, set1(2.0), zero)), mone(lt(df, zero)));
        const VD e = sub(m, isq(sub(sub(vR2, mul(sa, sa)), mul(sb, sb))));
        sti(&B.face[w], fnum);
        sti(&B.a[w], sa);
        sti(&B.b[w], sb);
        sti(&B.e[w], e);
        // --- H2 circle: k = argmin |d| (first on ties), (i, j) = (k+1, k+2) cyclic
        const VD uu = add(add(mul(d0, d0), mul(d1, d1)), mul(d2, d2));
        const VM broken = gt(uu, v4R2);
        VM k1 = lt(a1, a0);
        const VD mk = sel(k1, a1, a0);
        VM k2 = lt(a2, mk);
        k1 = mandn(k2, k1);
        k1 = mandn(broken, k1);
        k2 = mandn(broken, k2);  // broken molecule: k = 0
        const VD P = sel(broken, zero, rndi(mul(vRc2R, add(uu, vR2))));
        const VD uk = sel(k2, d2, sel(k1, d1, d0));
        const VD al = sel(k2, d0, sel(k1, d2, d1)), be = sel(k2, d1, sel(k1, d0, d2));
        const VD ck = sel(broken, zero, rndi(mul(mul(P, uk), vinvR2)));
        const VD A = add(mul(al, al), mul(be, be));
        const VD hk = sel(k2, h2, sel(k1, h1, h0)), hi = sel(k2, h0, sel(k1, h2, h1)),
                 hj = sel(k2, h1, sel(k1, h0, h2));
        const VD invA = dv(one, A);
        const VM nosolve = mor(mor(broken, le(A, zero)), mor(gt(mul(hk, hk), v4R2), gt(A, v4R2)));
        const VD L = sub(P, mul(uk, hk));
        const VD S = sub(vR2, mul(hk, hk));
        const VD sq = isq(sub(mul(S, A), mul(L, L)));
        const VD la = mul(L, al), lb = mul(L, be), bs = mul(be, sq), as = mul(al, sq);
        const VD x0 = sel(nosolve, zero, rndi(mul(add(la, bs), invA)));
        const VD y0 = sel(nosolve, zero, rndi(mul(sub(lb, as), invA)));
        const VD x1 = sel(nosolve, zero, rndi(mul(sub(la, bs), invA)));
        const VD y1 = sel(nosolve, zero, rndi(mul(add(lb, as), invA)));
        const VD e0 = add(vabs(sub(hi, x0)), vabs(sub(hj, y0)));
        const VD e1 = add(vabs(sub(hi, x1)), vabs(sub(hj, y1)));
        const VM sg = lt(e1, e0);
        sti(&B.v[w], sub(hk, ck));
        sti(&B.side[w], mone(sg));
        sti(&B.r1[w], sub(hi, sel(sg, x1, x0)));
        sti(&B.r2[w], sub(hj, sel(sg, y1, y0)));
        if (four_site) {
            const VD q0 = ld(&B.mx[w]), q1 = ld(&B.my[w]), q2 = ld(&B.mz[w]);
            sti(&B.m0[w], sub(q0, rndi(mul(vvs, add(d0, h0)))));
            sti(&B.m1[w], sub(q1, rndi(mul(vvs, add(d1, h1)))));
            sti(&B.m2[w], sub(q2, rndi(mul(vvs, add(d2, h2)))));
        }
    }
}

// sphere code of bonded atoms around their partner, radius per atom (as sphere_encode)
SZ3_BIOMD_TARGET static void sphere_batch(SphereBatch &B, size_t n) {
    const VD zero = set1(0.0);
    for (size_t w = 0; w < n; w += NL) {
        const VD d0 = ld(&B.dx[w]), d1 = ld(&B.dy[w]), d2 = ld(&B.dz[w]), r2 = ld(&B.r2[w]);
        const VD a0 = vabs(d0), a1 = vabs(d1), a2 = vabs(d2);
        const VM c1 = gt(a1, a0);
        VD m = sel(c1, a1, a0);
        const VM c2 = gt(a2, m);
        m = sel(c2, a2, m);
        const VM f1 = mandn(c2, c1), f2 = c2;
        const VD df = sel(f2, d2, sel(f1, d1, d0));
        const VD sa = sel(f2, d0, sel(f1, d2, d1)), sb = sel(f2, d1, sel(f1, d0, d2));
        const VD fnum = add(sel(f2, set1(4.0), sel(f1, set1(2.0), zero)), mone(lt(df, zero)));
        sti(&B.face[w], fnum);
        sti(&B.a[w], sa);
        sti(&B.b[w], sb);
        sti(&B.e[w], sub(m, isq(sub(sub(r2, mul(sa, sa)), mul(sb, sb)))));
    }
}

}  // namespace k_avx2
#undef SZ3_BIOMD_TARGET
#undef SZ3_BIOMD_H
#endif

// the AVX2 kernels if avx2 and this CPU runs them, else null (use the scalar code)
struct Kernels {
    void (*water)(WaterBatch &, size_t, int64_t, double, double, bool, double) = nullptr;
    void (*sphere)(SphereBatch &, size_t) = nullptr;
};
inline Kernels simd_kernels(bool avx2) {
    Kernels k;
#if defined(SZ3_BIOMD_X86_DISPATCH)
    if (avx2 && cpu_avx2()) {
        k.water = k_avx2::water_intra;
        k.sphere = k_avx2::sphere_batch;
    }
#else
    (void)avx2;
#endif
    return k;
}

}  // namespace biomd
}  // namespace SZ3
#endif
