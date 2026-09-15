/**
 * @file zfp_codec.hpp
 * @brief zfp CODEC 5's block codec, instantiated so the decorrelating transform and the
 *        embedded coder can be driven as two separate stages.
 *
 * Upstream expands its block codec through macros -- one translation unit per (scalar
 * type, dimensionality) -- and keeps every stage function `static`. The public API
 * therefore begins at `zfp_encode_block_*`, which carries a block all the way from
 * floating point to bitstream, past the point a composable pipeline has to cut. Rather
 * than strip `static` from upstream, this header runs the same macro expansion once per
 * pair inside its own namespace and adds wrappers at the stage boundary.
 *
 * The seam sits where zfp already divides the work: everything up to and including
 * `fwd_order` carries no error bound, and the bound enters at `encode_ints` by way of
 * the per-block precision. That holds for ZFP_ROUND_NEVER, which is upstream's default;
 * with ZFP_ROUND_FIRST the transform would need the precision too and the seam moves.
 */

#ifndef SZ3_ZFP_CODEC_HPP
#define SZ3_ZFP_CODEC_HPP

// Upstream C compiled as C++, against SZ3's warning flags.
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC system_header
#endif

// Every name undefined at the end of this header is saved here first: an installed zfp
// defines most of them, and undefining them would break it for whoever includes it next.
#pragma push_macro("SZ3_ZFP_OPS")
#pragma push_macro("MIN")
#pragma push_macro("MAX")
#pragma push_macro("inline_")
#pragma push_macro("restrict_")
#pragma push_macro("extern_")
#pragma push_macro("cache_align_")
#pragma push_macro("ZFP_MIN_BITS")
#pragma push_macro("ZFP_MAX_BITS")
#pragma push_macro("ZFP_MAX_PREC")
#pragma push_macro("ZFP_MIN_EXP")
#pragma push_macro("ZFP_ROUNDING_MODE")
#pragma push_macro("ZFP_ROUND_FIRST")
#pragma push_macro("ZFP_ROUND_NEVER")
#pragma push_macro("ZFP_ROUND_LAST")
#pragma push_macro("_t1")
#pragma push_macro("_t2")
#pragma push_macro("_cat2")
#pragma push_macro("_cat3")
#pragma push_macro("wsize")
#pragma push_macro("ZFP_CACHE_LINE_SIZE")
#pragma push_macro("INT64C")
#pragma push_macro("UINT64C")
#pragma push_macro("INT64PRId")
#pragma push_macro("INT64PRIi")
#pragma push_macro("UINT64PRIo")
#pragma push_macro("UINT64PRIu")
#pragma push_macro("UINT64PRIx")
#pragma push_macro("INT64SCNd")
#pragma push_macro("INT64SCNi")
#pragma push_macro("UINT64SCNo")
#pragma push_macro("UINT64SCNu")
#pragma push_macro("UINT64SCNx")

// The vendored sources include these from inside the namespace below; pull them in at
// global scope first so their include guards make those directives inert.
#include <cfloat>
#include <cinttypes>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <float.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

namespace SZ3 {
namespace ZFP {

// inline.h keys off __STDC_VERSION__, which C++ never defines, so `inline_` would become
// plain `static` and every unused helper would warn.
#define inline_ static inline
#include "include/zfp/internal/zfp/inline.h"
#include "include/zfp/internal/zfp/types.h"
#include "include/zfp/internal/zfp/macros.h"
#include "include/zfp/bitstream.h"
#include "include/zfp/bitstream.inl"
#include "src/template/template.h"

// zfp's rounding modes. These live in zfp.h, which is not vendored here; without them
// the `#if ZFP_ROUNDING_MODE == ZFP_ROUND_LAST` guards inside the codec compare two
// undefined identifiers, both of which the preprocessor reads as 0, and the rounding
// they gate is applied when upstream's default leaves it off.
#define ZFP_ROUND_FIRST (-1)
#define ZFP_ROUND_NEVER 0
#define ZFP_ROUND_LAST  1
#ifndef ZFP_ROUNDING_MODE
  #define ZFP_ROUNDING_MODE ZFP_ROUND_NEVER
#endif

// zfp's stream limits, which codec.h and the stage functions below expect.
#define ZFP_MIN_BITS     1
#define ZFP_MAX_BITS 16658
#define ZFP_MAX_PREC    64
#define ZFP_MIN_EXP  -1074

/* The block codec reads exactly these five fields off a zfp_stream. SZ3 drives the
 * stages directly and never constructs one; this shim exists only so that upstream's
 * block entry points still compile without vendoring the whole C API. */
typedef struct {
    uint minbits;
    uint maxbits;
    uint maxprec;
    int minexp;
    bitstream* stream;
} zfp_stream;

namespace f1 {
namespace {  // upstream's public block entry points are not static, and this header
             // reaches many translation units
#define SZ3_ZFP_BLOCK      "src/block1.h"
#define SZ3_ZFP_TRAITS     "src/traitsf.h"
#define SZ3_ZFP_CODEC_N    "src/template/codec1.c"
#define SZ3_ZFP_ENCODE_N   "src/template/encode1.c"
#define SZ3_ZFP_DECODE_N   "src/template/decode1.c"
#define SZ3_ZFP_REVENCODE_N "src/template/revencode1.c"
#define SZ3_ZFP_REVDECODE_N "src/template/revdecode1.c"
#include "zfp_inst.inc"
}  // anonymous
}  // namespace f1

namespace d1 {
namespace {  // upstream's public block entry points are not static, and this header
             // reaches many translation units
#define SZ3_ZFP_BLOCK      "src/block1.h"
#define SZ3_ZFP_TRAITS     "src/traitsd.h"
#define SZ3_ZFP_CODEC_N    "src/template/codec1.c"
#define SZ3_ZFP_ENCODE_N   "src/template/encode1.c"
#define SZ3_ZFP_DECODE_N   "src/template/decode1.c"
#define SZ3_ZFP_REVENCODE_N "src/template/revencode1.c"
#define SZ3_ZFP_REVDECODE_N "src/template/revdecode1.c"
#include "zfp_inst.inc"
}  // anonymous
}  // namespace d1

namespace f2 {
namespace {  // upstream's public block entry points are not static, and this header
             // reaches many translation units
#define SZ3_ZFP_BLOCK      "src/block2.h"
#define SZ3_ZFP_TRAITS     "src/traitsf.h"
#define SZ3_ZFP_CODEC_N    "src/template/codec2.c"
#define SZ3_ZFP_ENCODE_N   "src/template/encode2.c"
#define SZ3_ZFP_DECODE_N   "src/template/decode2.c"
#define SZ3_ZFP_REVENCODE_N "src/template/revencode2.c"
#define SZ3_ZFP_REVDECODE_N "src/template/revdecode2.c"
#include "zfp_inst.inc"
}  // anonymous
}  // namespace f2

namespace d2 {
namespace {  // upstream's public block entry points are not static, and this header
             // reaches many translation units
#define SZ3_ZFP_BLOCK      "src/block2.h"
#define SZ3_ZFP_TRAITS     "src/traitsd.h"
#define SZ3_ZFP_CODEC_N    "src/template/codec2.c"
#define SZ3_ZFP_ENCODE_N   "src/template/encode2.c"
#define SZ3_ZFP_DECODE_N   "src/template/decode2.c"
#define SZ3_ZFP_REVENCODE_N "src/template/revencode2.c"
#define SZ3_ZFP_REVDECODE_N "src/template/revdecode2.c"
#include "zfp_inst.inc"
}  // anonymous
}  // namespace d2

namespace f3 {
namespace {  // upstream's public block entry points are not static, and this header
             // reaches many translation units
#define SZ3_ZFP_BLOCK      "src/block3.h"
#define SZ3_ZFP_TRAITS     "src/traitsf.h"
#define SZ3_ZFP_CODEC_N    "src/template/codec3.c"
#define SZ3_ZFP_ENCODE_N   "src/template/encode3.c"
#define SZ3_ZFP_DECODE_N   "src/template/decode3.c"
#define SZ3_ZFP_REVENCODE_N "src/template/revencode3.c"
#define SZ3_ZFP_REVDECODE_N "src/template/revdecode3.c"
#include "zfp_inst.inc"
}  // anonymous
}  // namespace f3

namespace d3 {
namespace {  // upstream's public block entry points are not static, and this header
             // reaches many translation units
#define SZ3_ZFP_BLOCK      "src/block3.h"
#define SZ3_ZFP_TRAITS     "src/traitsd.h"
#define SZ3_ZFP_CODEC_N    "src/template/codec3.c"
#define SZ3_ZFP_ENCODE_N   "src/template/encode3.c"
#define SZ3_ZFP_DECODE_N   "src/template/decode3.c"
#define SZ3_ZFP_REVENCODE_N "src/template/revencode3.c"
#define SZ3_ZFP_REVDECODE_N "src/template/revdecode3.c"
#include "zfp_inst.inc"
}  // anonymous
}  // namespace d3

// These call the instantiations above, which have internal linkage because upstream's stage
// functions are static. So these must have it too: with external linkage, several TUs would
// define one symbol over different copies of those functions. The templates that call these
// are external and so formally name a per-TU entity; no diagnostic is required and the
// definitions are identical, every TU having compiled this same header.
namespace {

/**
 * @brief Dispatches a scalar type to its instantiation, since a namespace cannot be a
 *        template argument. 1D, 2D and 3D in float and double; 4D would be the same
 *        pattern over upstream's block4 sources.
 */
template <class T, uint N>
struct ops;

#define SZ3_ZFP_OPS(T, D, NS)                                                                                \
    template <>                                                                                              \
    struct ops<T, D> {                                                                                       \
        using uint_type = NS::uint_type;                                                                     \
        static constexpr int block_size = NS::block_size;                                                    \
        static constexpr int ebias = NS::ebias;                                                              \
        static int fwd(const T *p, const ptrdiff_t *s, const size_t *n, uint_type *ub) {                     \
            return NS::fwd_block(p, s, n, ub);                                                               \
        }                                                                                                    \
        static void inv(T *p, const ptrdiff_t *s, const size_t *n, int emax, const uint_type *ub) {          \
            NS::inv_block(p, s, n, emax, ub);                                                                \
        }                                                                                                    \
        static unsigned enc(bitstream *s, int emax, const uint_type *ub, int minexp) {                       \
            return NS::enc_block(s, emax, ub, ZFP_MAX_PREC, minexp, ZFP_MIN_BITS, ZFP_MAX_BITS);             \
        }                                                                                                    \
        static unsigned dec(bitstream *s, int *emax, uint_type *ub, int minexp) {                            \
            return NS::dec_block(s, emax, ub, ZFP_MAX_PREC, minexp, ZFP_MIN_BITS, ZFP_MAX_BITS);             \
        }                                                                                                    \
    }

SZ3_ZFP_OPS(float, 1, f1);
SZ3_ZFP_OPS(double, 1, d1);
SZ3_ZFP_OPS(float, 2, f2);
SZ3_ZFP_OPS(double, 2, d2);
SZ3_ZFP_OPS(float, 3, f3);
SZ3_ZFP_OPS(double, 3, d3);
#undef SZ3_ZFP_OPS

/// The exponent floor that realises an absolute error bound, as zfp_stream_set_accuracy does.
inline int minexp_for(double tolerance) {
    int emin = ZFP_MIN_EXP;
    if (tolerance > 0) {
        frexp(tolerance, &emin);
        emin--;
    }
    return emin;
}

// Macros have no namespace, so anything still defined here leaks. The per-instantiation
// ones (Scalar, Int, DIMS, BLOCK_SIZE, PERM) are undefined by zfp_inst.inc.
#undef MIN
#undef MAX
#undef inline_
#undef restrict_
#undef extern_
#undef cache_align_
#undef ZFP_MIN_BITS
#undef ZFP_MAX_BITS
#undef ZFP_MAX_PREC
#undef ZFP_MIN_EXP
#undef ZFP_ROUNDING_MODE
#undef ZFP_ROUND_FIRST
#undef ZFP_ROUND_NEVER
#undef ZFP_ROUND_LAST
#undef _t1
#undef _t2
#undef _cat2
#undef _cat3
#undef wsize
#undef ZFP_CACHE_LINE_SIZE
#undef INT64C
#undef UINT64C
#undef INT64PRId
#undef INT64PRIi
#undef UINT64PRIo
#undef UINT64PRIu
#undef UINT64PRIx
#undef INT64SCNd
#undef INT64SCNi
#undef UINT64SCNo
#undef UINT64SCNu
#undef UINT64SCNx

// Restore whatever the user had.
#pragma pop_macro("SZ3_ZFP_OPS")
#pragma pop_macro("MIN")
#pragma pop_macro("MAX")
#pragma pop_macro("inline_")
#pragma pop_macro("restrict_")
#pragma pop_macro("extern_")
#pragma pop_macro("cache_align_")
#pragma pop_macro("ZFP_MIN_BITS")
#pragma pop_macro("ZFP_MAX_BITS")
#pragma pop_macro("ZFP_MAX_PREC")
#pragma pop_macro("ZFP_MIN_EXP")
#pragma pop_macro("ZFP_ROUNDING_MODE")
#pragma pop_macro("ZFP_ROUND_FIRST")
#pragma pop_macro("ZFP_ROUND_NEVER")
#pragma pop_macro("ZFP_ROUND_LAST")
#pragma pop_macro("_t1")
#pragma pop_macro("_t2")
#pragma pop_macro("_cat2")
#pragma pop_macro("_cat3")
#pragma pop_macro("wsize")
#pragma pop_macro("ZFP_CACHE_LINE_SIZE")
#pragma pop_macro("INT64C")
#pragma pop_macro("UINT64C")
#pragma pop_macro("INT64PRId")
#pragma pop_macro("INT64PRIi")
#pragma pop_macro("UINT64PRIo")
#pragma pop_macro("UINT64PRIu")
#pragma pop_macro("UINT64PRIx")
#pragma pop_macro("INT64SCNd")
#pragma pop_macro("INT64SCNi")
#pragma pop_macro("UINT64SCNo")
#pragma pop_macro("UINT64SCNu")
#pragma pop_macro("UINT64SCNx")

}  // anonymous

}  // namespace ZFP
}  // namespace SZ3
#endif
