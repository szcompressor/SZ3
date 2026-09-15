/* Instantiate zfp's block codec for one (Scalar, DIMS) pair and expose the stage
 * boundary between the decorrelating transform and the embedded coder.
 *
 * Include once per namespace with SZ3_ZFP_TRAITS and SZ3_ZFP_BLOCK set. Upstream's template
 * sources are included unmodified; the wrappers at the bottom are the only additions.
 */
#include SZ3_ZFP_BLOCK
#include SZ3_ZFP_TRAITS
#include "src/template/codec.h"
#include "src/template/codec.c"
#include "src/template/codecf.c"
#include SZ3_ZFP_CODEC_N
#include "src/template/encode.c"
#include "src/template/encodef.c"
#include SZ3_ZFP_ENCODE_N
#include "src/template/revcodecf.c"
#include "src/template/revencode.c"
#include "src/template/revencodef.c"
#include SZ3_ZFP_REVENCODE_N
#include "src/template/decode.c"
#include "src/template/decodef.c"
#include SZ3_ZFP_DECODE_N
#include "src/template/revdecode.c"
#include "src/template/revdecodef.c"
#include SZ3_ZFP_REVDECODE_N

/* ---------------- the seam ---------------- */

typedef Scalar scalar_type;
typedef Int    int_type;
typedef UInt   uint_type;

enum {
    block_size = BLOCK_SIZE,
    ebits      = EBITS,
    ebias      = EBIAS,
    dims       = DIMS
};

/* Transform stage: gather a (possibly partial) block, align it to its own exponent, run
 * the decorrelating transform, reorder and map to negabinary. Carries no error bound --
 * that enters downstream, where each block's precision is derived. */
inline int fwd_block(const Scalar* p, const ptrdiff_t* s, const size_t* n, UInt* ublock) {
    Scalar fblock[BLOCK_SIZE];
#if DIMS == 1
    if (n[0] == 4) _t2(gather, Scalar, 1)(fblock, p, s[0]);
    else           _t2(gather_partial, Scalar, 1)(fblock, p, n[0], s[0]);
#elif DIMS == 2
    if (n[0] == 4 && n[1] == 4) _t2(gather, Scalar, 2)(fblock, p, s[0], s[1]);
    else                        _t2(gather_partial, Scalar, 2)(fblock, p, n[0], n[1], s[0], s[1]);
#else
    if (n[0] == 4 && n[1] == 4 && n[2] == 4) _t2(gather, Scalar, 3)(fblock, p, s[0], s[1], s[2]);
    else _t2(gather_partial, Scalar, 3)(fblock, p, n[0], n[1], n[2], s[0], s[1], s[2]);
#endif
    int emax = _t1(exponent_block, Scalar)(fblock, BLOCK_SIZE);
    /* fwd_cast scales by quantize(1, emax) = LDEXP(1, CHAR_BIT*sizeof(Scalar) - 2 - emax), which
       is finite only while that exponent stays within EBIAS. Below this emax it is infinity, and
       infinity times a coefficient is either a NaN or an out-of-range product, neither of which
       has a defined conversion to Int. Upstream reaches the same code but guards it with
       `maxprec != 0` (encode.c), which needs the error bound -- the thing this stage deliberately
       does not carry. The exponent alone decides the overflow, so guard on that instead. It
       covers an all-zero block, whose emax is -EBIAS. enc_block writes such a block from emax
       alone and never reads ublock, so zeroing it here changes nothing it can observe. */
    if (emax < (int)(CHAR_BIT * sizeof(Scalar)) - 2 - EBIAS) {
        uint i;
        for (i = 0; i < BLOCK_SIZE; i++) ublock[i] = 0;
        return emax;
    }
    Int iblock[BLOCK_SIZE];
    _t1(fwd_cast, Scalar)(iblock, fblock, BLOCK_SIZE, emax);
    _t2(fwd_xform, Int, DIMS)(iblock);
    _t1(fwd_order, Int)(ublock, iblock, PERM, BLOCK_SIZE);
    return emax;
}

/* Inverse transform stage. */
inline void inv_block(Scalar* p, const ptrdiff_t* s, const size_t* n, int emax, const UInt* ublock) {
    Scalar fblock[BLOCK_SIZE];
    Int iblock[BLOCK_SIZE];
    _t1(inv_order, Int)(ublock, iblock, PERM, BLOCK_SIZE);
    _t2(inv_xform, Int, DIMS)(iblock);
    _t1(inv_cast, Scalar)(iblock, fblock, BLOCK_SIZE, emax);
#if DIMS == 1
    if (n[0] == 4) _t2(scatter, Scalar, 1)(fblock, p, s[0]);
    else           _t2(scatter_partial, Scalar, 1)(fblock, p, n[0], s[0]);
#elif DIMS == 2
    if (n[0] == 4 && n[1] == 4) _t2(scatter, Scalar, 2)(fblock, p, s[0], s[1]);
    else                        _t2(scatter_partial, Scalar, 2)(fblock, p, n[0], n[1], s[0], s[1]);
#else
    if (n[0] == 4 && n[1] == 4 && n[2] == 4) _t2(scatter, Scalar, 3)(fblock, p, s[0], s[1], s[2]);
    else _t2(scatter_partial, Scalar, 3)(fblock, p, n[0], n[1], n[2], s[0], s[1], s[2]);
#endif
}

/* Coder stage: the error bound enters here, as minexp. Mirrors upstream's
 * encode_block(Scalar, DIMS) with the transform already done. */
inline unsigned enc_block(bitstream* s, int emax, const UInt* ublock,
                          unsigned maxprec_cap, int minexp,
                          unsigned zminbits, unsigned zmaxbits) {
    unsigned bits = 1;
    unsigned maxprec = precision(emax, maxprec_cap, minexp, DIMS);
    unsigned e = maxprec ? (unsigned)(emax + EBIAS) : 0;
    if (e) {
        bits += EBITS;
        stream_write_bits(s, 2 * e + 1, bits);
        unsigned mb = zminbits - MIN(bits, zminbits);
        unsigned Mb = zmaxbits - bits;
        unsigned b = _t1(encode_ints, UInt)(s, Mb, maxprec, ublock, BLOCK_SIZE);
        if (b < mb) { stream_pad(s, mb - b); b = mb; }
        bits += b;
    } else {
        stream_write_bit(s, 0);
        if (zminbits > bits) { stream_pad(s, zminbits - bits); bits = zminbits; }
    }
    return bits;
}

/* Inverse coder stage. A zero block comes back as emax = -EBIAS with zero coefficients,
 * which the inverse transform turns into zeros. */
inline unsigned dec_block(bitstream* s, int* emax_out, UInt* ublock,
                          unsigned maxprec_cap, int minexp,
                          unsigned zminbits, unsigned zmaxbits) {
    unsigned bits = 1;
    if (stream_read_bit(s)) {
        bits += EBITS;
        int emax = (int)stream_read_bits(s, EBITS) - EBIAS;
        unsigned maxprec = precision(emax, maxprec_cap, minexp, DIMS);
        unsigned mb = zminbits - MIN(bits, zminbits);
        unsigned Mb = zmaxbits - bits;
        unsigned b = _t1(decode_ints, UInt)(s, Mb, maxprec, ublock, BLOCK_SIZE);
        if (b < mb) { stream_skip(s, mb - b); b = mb; }
        bits += b;
        *emax_out = emax;
    } else {
        unsigned i;
        for (i = 0; i < BLOCK_SIZE; i++) ublock[i] = 0;
        *emax_out = -EBIAS;
        if (zminbits > bits) { stream_skip(s, zminbits - bits); bits = zminbits; }
    }
    return bits;
}

#undef Scalar
#undef Int
#undef UInt
#undef EBITS
#undef PBITS
#undef NBMASK
#undef TCMASK
#undef SCALAR_MIN
#undef FABS
#undef FREXP
#undef LDEXP
#undef DIMS
#undef PERM
#undef BLOCK_SIZE
#undef EBIAS
#undef REVERSIBLE
#undef SZ3_ZFP_TRAITS
#undef SZ3_ZFP_BLOCK
#undef SZ3_ZFP_CODEC_N
#undef SZ3_ZFP_ENCODE_N
#undef SZ3_ZFP_DECODE_N
#undef SZ3_ZFP_REVENCODE_N
#undef SZ3_ZFP_REVDECODE_N
