#ifndef SZ3_MDC_SIMD_HPP
#define SZ3_MDC_SIMD_HPP

// ALGO_MDC, part 2: SIMD kernels for SSE4.1, AVX2 and AVX-512 (picked at run time) and NEON (the aarch64 baseline).
// The kernel bodies are in MDCSimdKernels.inc, compiled once per instruction set in its own namespace over a small set
// of helpers; they equal the scalar code bit for bit. The scalar code is the fallback everywhere else.

#include "SZ3/compressor/specialized/mdc/MDCCore.hpp"
#if defined(SZ3_MDC_X86_DISPATCH)
#include <immintrin.h>
#endif
#if defined(__aarch64__) && defined(__ARM_NEON) && (defined(__GNUC__) || defined(__clang__))
#include <arm_neon.h>
#define SZ3_MDC_NEON_KERNELS 1
#endif

namespace SZ3 {
namespace mdc {

constexpr size_t SIMD_PAD = 8;  // widest vector, in doubles

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

// instruction-set levels: 0 scalar, 1 SSE4.1 or NEON (2 lanes), 2 AVX2 (4), 3 AVX-512 (8)
enum { SIMD_SCALAR = 0, SIMD_128 = 1, SIMD_AVX2 = 2, SIMD_AVX512 = 3 };
inline int simd_level() {
#if defined(SZ3_MDC_X86_DISPATCH)
    static const int lv = __builtin_cpu_supports("avx512f")  ? SIMD_AVX512
                          : __builtin_cpu_supports("avx2")   ? SIMD_AVX2
                          : __builtin_cpu_supports("sse4.1") ? SIMD_128
                                                             : SIMD_SCALAR;
    return lv;
#elif defined(SZ3_MDC_NEON_KERNELS)
    return SIMD_128;
#else
    return SIMD_SCALAR;
#endif
}

#ifdef SZ3_MDC_X86_DISPATCH
namespace k_sse41 {
#define SZ3_MDC_TARGET __attribute__((target("sse4.1")))
#define SZ3_MDC_H __attribute__((target("sse4.1"), always_inline)) static inline
constexpr size_t NL = 2;
typedef __m128d VD;
typedef __m128d VM;
SZ3_MDC_H VD set1(double x) { return _mm_set1_pd(x); }
SZ3_MDC_H VD ld(const double *p) { return _mm_loadu_pd(p); }
SZ3_MDC_H VD add(VD a, VD b) { return _mm_add_pd(a, b); }
SZ3_MDC_H VD sub(VD a, VD b) { return _mm_sub_pd(a, b); }
SZ3_MDC_H VD mul(VD a, VD b) { return _mm_mul_pd(a, b); }
SZ3_MDC_H VD dv(VD a, VD b) { return _mm_div_pd(a, b); }
SZ3_MDC_H VD vabs(VD x) { return _mm_andnot_pd(_mm_set1_pd(-0.0), x); }
SZ3_MDC_H VM gt(VD a, VD b) { return _mm_cmpgt_pd(a, b); }
SZ3_MDC_H VM lt(VD a, VD b) { return _mm_cmplt_pd(a, b); }
SZ3_MDC_H VM le(VD a, VD b) { return _mm_cmple_pd(a, b); }
SZ3_MDC_H VM mor(VM a, VM b) { return _mm_or_pd(a, b); }
SZ3_MDC_H VM mandn(VM a, VM b) { return _mm_andnot_pd(a, b); }  // ~a & b
SZ3_MDC_H VD sel(VM m, VD t, VD f) { return _mm_blendv_pd(f, t, m); }
SZ3_MDC_H VD mone(VM m) { return _mm_and_pd(m, _mm_set1_pd(1.0)); }
SZ3_MDC_H VD rndi(VD y) { return _mm_cvtepi32_pd(_mm_cvtpd_epi32(y)); }  // nearest, ties to even, as rnd_pred
SZ3_MDC_H VD isq(VD n) {  // n <= 0 ? 0 : trunc(sqrt(n) + 0.5), as isqrt_round
    const VD r = _mm_round_pd(_mm_add_pd(_mm_sqrt_pd(_mm_max_pd(n, _mm_setzero_pd())), _mm_set1_pd(0.5)),
                              _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
    return sel(le(n, _mm_setzero_pd()), _mm_setzero_pd(), r);
}
SZ3_MDC_H void sti(int32_t *o, VD x) { _mm_storel_epi64(reinterpret_cast<__m128i *>(o), _mm_cvtpd_epi32(x)); }
#include "SZ3/compressor/specialized/mdc/MDCSimdKernels.inc"
#undef SZ3_MDC_TARGET
#undef SZ3_MDC_H
}  // namespace k_sse41

namespace k_avx2 {
#define SZ3_MDC_TARGET __attribute__((target("avx2")))
#define SZ3_MDC_H __attribute__((target("avx2"), always_inline)) static inline
constexpr size_t NL = 4;
typedef __m256d VD;
typedef __m256d VM;
SZ3_MDC_H VD set1(double x) { return _mm256_set1_pd(x); }
SZ3_MDC_H VD ld(const double *p) { return _mm256_loadu_pd(p); }
SZ3_MDC_H VD add(VD a, VD b) { return _mm256_add_pd(a, b); }
SZ3_MDC_H VD sub(VD a, VD b) { return _mm256_sub_pd(a, b); }
SZ3_MDC_H VD mul(VD a, VD b) { return _mm256_mul_pd(a, b); }
SZ3_MDC_H VD dv(VD a, VD b) { return _mm256_div_pd(a, b); }
SZ3_MDC_H VD vabs(VD x) { return _mm256_andnot_pd(_mm256_set1_pd(-0.0), x); }
SZ3_MDC_H VM gt(VD a, VD b) { return _mm256_cmp_pd(a, b, _CMP_GT_OQ); }
SZ3_MDC_H VM lt(VD a, VD b) { return _mm256_cmp_pd(a, b, _CMP_LT_OQ); }
SZ3_MDC_H VM le(VD a, VD b) { return _mm256_cmp_pd(a, b, _CMP_LE_OQ); }
SZ3_MDC_H VM mor(VM a, VM b) { return _mm256_or_pd(a, b); }
SZ3_MDC_H VM mandn(VM a, VM b) { return _mm256_andnot_pd(a, b); }
SZ3_MDC_H VD sel(VM m, VD t, VD f) { return _mm256_blendv_pd(f, t, m); }
SZ3_MDC_H VD mone(VM m) { return _mm256_and_pd(m, _mm256_set1_pd(1.0)); }
SZ3_MDC_H VD rndi(VD y) { return _mm256_cvtepi32_pd(_mm256_cvtpd_epi32(y)); }
SZ3_MDC_H VD isq(VD n) {
    const VD zero = _mm256_setzero_pd();
    const VD r = _mm256_round_pd(_mm256_add_pd(_mm256_sqrt_pd(_mm256_max_pd(n, zero)), _mm256_set1_pd(0.5)),
                                 _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
    return sel(le(n, zero), zero, r);
}
SZ3_MDC_H void sti(int32_t *o, VD x) { _mm_storeu_si128(reinterpret_cast<__m128i *>(o), _mm256_cvtpd_epi32(x)); }
#include "SZ3/compressor/specialized/mdc/MDCSimdKernels.inc"
#undef SZ3_MDC_TARGET
#undef SZ3_MDC_H
}  // namespace k_avx2

namespace k_avx512 {
#define SZ3_MDC_TARGET __attribute__((target("avx512f")))
#define SZ3_MDC_H __attribute__((target("avx512f"), always_inline)) static inline
constexpr size_t NL = 8;
typedef __m512d VD;
typedef __mmask8 VM;
SZ3_MDC_H VD set1(double x) { return _mm512_set1_pd(x); }
SZ3_MDC_H VD ld(const double *p) { return _mm512_loadu_pd(p); }
SZ3_MDC_H VD add(VD a, VD b) { return _mm512_add_pd(a, b); }
SZ3_MDC_H VD sub(VD a, VD b) { return _mm512_sub_pd(a, b); }
SZ3_MDC_H VD mul(VD a, VD b) { return _mm512_mul_pd(a, b); }
SZ3_MDC_H VD dv(VD a, VD b) { return _mm512_div_pd(a, b); }
SZ3_MDC_H VD vabs(VD x) { return _mm512_abs_pd(x); }
SZ3_MDC_H VM gt(VD a, VD b) { return _mm512_cmp_pd_mask(a, b, _CMP_GT_OQ); }
SZ3_MDC_H VM lt(VD a, VD b) { return _mm512_cmp_pd_mask(a, b, _CMP_LT_OQ); }
SZ3_MDC_H VM le(VD a, VD b) { return _mm512_cmp_pd_mask(a, b, _CMP_LE_OQ); }
SZ3_MDC_H VM mor(VM a, VM b) { return VM(a | b); }
SZ3_MDC_H VM mandn(VM a, VM b) { return VM(~a & b); }
SZ3_MDC_H VD sel(VM m, VD t, VD f) { return _mm512_mask_blend_pd(m, f, t); }
SZ3_MDC_H VD mone(VM m) { return _mm512_maskz_mov_pd(m, _mm512_set1_pd(1.0)); }
SZ3_MDC_H VD rndi(VD y) { return _mm512_cvtepi32_pd(_mm512_cvtpd_epi32(y)); }
SZ3_MDC_H VD isq(VD n) {
    const VD zero = _mm512_setzero_pd();
    const VD r = _mm512_roundscale_pd(_mm512_add_pd(_mm512_sqrt_pd(_mm512_max_pd(n, zero)), _mm512_set1_pd(0.5)),
                                      _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
    return sel(le(n, zero), zero, r);
}
SZ3_MDC_H void sti(int32_t *o, VD x) { _mm256_storeu_si256(reinterpret_cast<__m256i *>(o), _mm512_cvtpd_epi32(x)); }
#include "SZ3/compressor/specialized/mdc/MDCSimdKernels.inc"
#undef SZ3_MDC_TARGET
#undef SZ3_MDC_H
}  // namespace k_avx512
#endif

#ifdef SZ3_MDC_NEON_KERNELS
namespace k_neon {
#define SZ3_MDC_TARGET
#define SZ3_MDC_H __attribute__((always_inline)) static inline
constexpr size_t NL = 2;
typedef float64x2_t VD;
typedef uint64x2_t VM;
SZ3_MDC_H VD set1(double x) { return vdupq_n_f64(x); }
SZ3_MDC_H VD ld(const double *p) { return vld1q_f64(p); }
SZ3_MDC_H VD add(VD a, VD b) { return vaddq_f64(a, b); }
SZ3_MDC_H VD sub(VD a, VD b) { return vsubq_f64(a, b); }
SZ3_MDC_H VD mul(VD a, VD b) { return vmulq_f64(a, b); }
SZ3_MDC_H VD dv(VD a, VD b) { return vdivq_f64(a, b); }
SZ3_MDC_H VD vabs(VD x) { return vabsq_f64(x); }
SZ3_MDC_H VM gt(VD a, VD b) { return vcgtq_f64(a, b); }
SZ3_MDC_H VM lt(VD a, VD b) { return vcltq_f64(a, b); }
SZ3_MDC_H VM le(VD a, VD b) { return vcleq_f64(a, b); }
SZ3_MDC_H VM mor(VM a, VM b) { return vorrq_u64(a, b); }
SZ3_MDC_H VM mandn(VM a, VM b) { return vbicq_u64(b, a); }  // b & ~a
SZ3_MDC_H VD sel(VM m, VD t, VD f) { return vbslq_f64(m, t, f); }
SZ3_MDC_H VD mone(VM m) { return vreinterpretq_f64_u64(vandq_u64(m, vreinterpretq_u64_f64(vdupq_n_f64(1.0)))); }
SZ3_MDC_H VD rndi(VD y) { return vcvtq_f64_s64(vcvtnq_s64_f64(y)); }  // nearest, ties to even, as rnd_pred
SZ3_MDC_H VD isq(VD n) {
    const VD zero = vdupq_n_f64(0.0);
    const VD r = vrndq_f64(vaddq_f64(vsqrtq_f64(vmaxq_f64(n, zero)), vdupq_n_f64(0.5)));
    return sel(le(n, zero), zero, r);
}
SZ3_MDC_H void sti(int32_t *o, VD x) { vst1_s32(o, vmovn_s64(vcvtq_s64_f64(x))); }
#include "SZ3/compressor/specialized/mdc/MDCSimdKernels.inc"
#undef SZ3_MDC_TARGET
#undef SZ3_MDC_H
}  // namespace k_neon
#endif

// the kernels of the widest instruction set up to max_level that this CPU runs; null: use the scalar code
struct Kernels {
    int level = SIMD_SCALAR;
    void (*water)(WaterBatch &, size_t, int64_t, double, double, bool, double) = nullptr;
    void (*sphere)(SphereBatch &, size_t) = nullptr;
};
inline Kernels simd_kernels(int max_level) {
    Kernels k;
    const int lv = std::min(max_level, simd_level());
#if defined(SZ3_MDC_X86_DISPATCH)
    if (lv >= SIMD_AVX512) {
        k.level = SIMD_AVX512;
        k.water = k_avx512::water_intra;
        k.sphere = k_avx512::sphere_batch;
    } else if (lv == SIMD_AVX2) {
        k.level = SIMD_AVX2;
        k.water = k_avx2::water_intra;
        k.sphere = k_avx2::sphere_batch;
    } else if (lv == SIMD_128) {
        k.level = SIMD_128;
        k.water = k_sse41::water_intra;
        k.sphere = k_sse41::sphere_batch;
    }
#elif defined(SZ3_MDC_NEON_KERNELS)
    if (lv >= SIMD_128) {
        k.level = SIMD_128;
        k.water = k_neon::water_intra;
        k.sphere = k_neon::sphere_batch;
    }
#else
    (void)lv;
#endif
    return k;
}

}  // namespace mdc
}  // namespace SZ3
#endif
