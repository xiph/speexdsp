#ifndef SPEEXDSP_TESTS_CHECKASM_RESAMPLE_WRAP_H
#define SPEEXDSP_TESTS_CHECKASM_RESAMPLE_WRAP_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include "config.h"
#include "arch.h"
#include "speex/speex_resampler.h" /* opaque SpeexResamplerState type */

/* These tests exercise the four *full* resampler kernels in
 * libspeexdsp/resample.c -- resampler_basic_{direct,interpolate}_{single,double}
 * -- end to end (sinc-table indexing + main loop + saturation + the inner
 * product), comparing a pure-C build of each against the SIMD build.
 *
 * The C-vs-SIMD difference is internal: each resampler_basic_* function selects
 * its inner-product kernel via #ifndef OVERRIDE_*. To get both a C and a SIMD
 * build of the same static function in one binary, resample.c is #included once
 * per kernel variant in its own translation unit (wrap_resample_c.c forces the
 * C kernels; wrap_resample_{sse,sse2,neon}.c keep the native SIMD kernels), each
 * exposing a uniquely named non-static wrapper. */

/* ------------- Per-routine SIMD-availability gates -------------
 * Today resample_neon.h overrides only inner_product_single, so under NEON only
 * resampler_basic_direct_single differs from the C reference. Flip more of these
 * on when resample_neon.h grows OVERRIDE_INTERPOLATE_PRODUCT_SINGLE / *_DOUBLE. */
#ifdef USE_NEON
#  define HAVE_NEON_DIRECT_SINGLE 1
/* #  define HAVE_NEON_INTERPOLATE_SINGLE 1 */
#endif

/* RVV overrides both single kernels (both modes) and both double kernels (float only). */
#ifdef USE_RVV
#  define HAVE_RVV_DIRECT_SINGLE 1
#  define HAVE_RVV_INTERPOLATE_SINGLE 1
#  ifndef FIXED_POINT
#    define HAVE_RVV_DIRECT_DOUBLE 1
#    define HAVE_RVV_INTERPOLATE_DOUBLE 1
#  endif
#endif

/* ------------- Resampler "kind" = (method, precision) -------------
 * update_filter() picks one of the four functions based on the sample-rate
 * ratio (direct vs interpolate) and, in floating point, quality>8 (single vs
 * double). Tests build a state then assert it landed in the expected kind. */
enum resample_kind {
    RESAMPLE_KIND_OTHER = 0,
    RESAMPLE_KIND_DIRECT_SINGLE,
    RESAMPLE_KIND_INTERPOLATE_SINGLE,
#ifndef FIXED_POINT
    RESAMPLE_KIND_DIRECT_DOUBLE,
    RESAMPLE_KIND_INTERPOLATE_DOUBLE,
#endif
};

/* ------------- State helpers (defined once, in wrap_resample_c.c) -------------
 * Only that TU has both the four static function addresses and the public init
 * API in scope, so it is the sole owner of state construction/inspection. The
 * harness treats SpeexResamplerState as opaque. */
SpeexResamplerState *resample_make_state(unsigned in_rate, unsigned out_rate, int quality);
/* Multi-channel state for the interleaved integration path (see below). */
SpeexResamplerState *resample_make_state_ch(unsigned in_rate, unsigned out_rate,
        int quality, unsigned channels);
void resample_destroy_state(SpeexResamplerState *st);
enum resample_kind resample_kind(const SpeexResamplerState *st);
unsigned resample_filt_len(const SpeexResamplerState *st);
unsigned resample_den_rate(const SpeexResamplerState *st);
unsigned resample_oversample(const SpeexResamplerState *st);

/* ------------- Functions under test -------------
 * Signature matches resampler_basic_func in resample.c. Each wrapper resets the
 * per-channel cursor (st->last_sample/samp_frac_num) before calling, so repeated
 * invocations -- in particular checkasm_bench_new's call loop -- re-do identical
 * full work rather than starting past the end of the input. */
int resampler_basic_direct_single_c(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len);
int resampler_basic_interpolate_single_c(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len);
#ifndef FIXED_POINT
int resampler_basic_direct_double_c(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len);
int resampler_basic_interpolate_double_c(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len);
#endif

#if defined(USE_SSE) && !defined(FIXED_POINT)
int resampler_basic_direct_single_sse(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len);
int resampler_basic_interpolate_single_sse(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len);
#endif
#if defined(USE_SSE2) && !defined(FIXED_POINT)
int resampler_basic_direct_double_sse2(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len);
int resampler_basic_interpolate_double_sse2(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len);
#endif
#ifdef HAVE_NEON_DIRECT_SINGLE
int resampler_basic_direct_single_neon(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len);
#endif
#ifdef USE_RVV
int resampler_basic_direct_single_rvv(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len);
int resampler_basic_interpolate_single_rvv(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len);
#ifndef FIXED_POINT
int resampler_basic_direct_double_rvv(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len);
int resampler_basic_interpolate_double_rvv(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len);
#endif
#endif

/* ------------- Integration benchmark: full resampler pipeline -------------
 * These drive the *public* speex_resampler_process_float / process_int end to
 * end, so they measure a whole sample-rate conversion (sinc-table build happens
 * at make-state time, outside the timed call). The kernel used inside the
 * pipeline is fixed in the state by update_filter at init time, so each build
 * needs its own state constructor: a C-built state runs the C kernels; a NEON-
 * or SSE-built state runs that ISA's kernels. (The SSE translation unit keeps
 * both USE_SSE and USE_SSE2, so its state runs SSE single-precision and SSE2
 * double-precision kernels -- i.e. the full library pipeline.)
 *
 * The _int variants additionally exercise the int16 output path, i.e. the NEON
 * WORD2INT override (saturate_float_to_16bit), which process_float and the
 * kernel tests never touch.
 *
 * The _il variants drive the *interleaved* multi-channel API
 * (speex_resampler_process_interleaved_*), exercising the per-channel loop and
 * (for int) the deinterleave/reinterleave that the single-channel path skips.
 * Their in_len/out_len are PER CHANNEL; the return is the per-channel output
 * count. State for these comes from resample_make_state*_ch(.., channels).
 *
 * Lengths are passed BY VALUE (process_* overwrites its in_len/out_len, so a
 * by-pointer signature would shrink the work on every benchmark iteration). Each
 * wrapper zeroes the filter memory and resets every channel's cursor first, so
 * every call -- the correctness pair and each benchmark iteration -- redoes
 * identical, deterministic work. Returns the number of output samples produced. */
#ifndef DISABLE_FLOAT_API
int resample_process_c(SpeexResamplerState *st, const float *in,
        spx_uint32_t in_len, float *out, spx_uint32_t out_len);
int resample_process_int_c(SpeexResamplerState *st, const spx_int16_t *in,
        spx_uint32_t in_len, spx_int16_t *out, spx_uint32_t out_len);
int resample_process_il_c(SpeexResamplerState *st, const float *in,
        spx_uint32_t in_len, float *out, spx_uint32_t out_len);
int resample_process_int_il_c(SpeexResamplerState *st, const spx_int16_t *in,
        spx_uint32_t in_len, spx_int16_t *out, spx_uint32_t out_len);
#ifdef USE_NEON
SpeexResamplerState *resample_make_state_neon(unsigned in_rate, unsigned out_rate, int quality);
SpeexResamplerState *resample_make_state_neon_ch(unsigned in_rate, unsigned out_rate,
        int quality, unsigned channels);
int resample_process_neon(SpeexResamplerState *st, const float *in,
        spx_uint32_t in_len, float *out, spx_uint32_t out_len);
int resample_process_int_neon(SpeexResamplerState *st, const spx_int16_t *in,
        spx_uint32_t in_len, spx_int16_t *out, spx_uint32_t out_len);
int resample_process_il_neon(SpeexResamplerState *st, const float *in,
        spx_uint32_t in_len, float *out, spx_uint32_t out_len);
int resample_process_int_il_neon(SpeexResamplerState *st, const spx_int16_t *in,
        spx_uint32_t in_len, spx_int16_t *out, spx_uint32_t out_len);
#endif
#if defined(USE_SSE) && !defined(FIXED_POINT)
SpeexResamplerState *resample_make_state_sse(unsigned in_rate, unsigned out_rate, int quality);
SpeexResamplerState *resample_make_state_sse_ch(unsigned in_rate, unsigned out_rate,
        int quality, unsigned channels);
int resample_process_sse(SpeexResamplerState *st, const float *in,
        spx_uint32_t in_len, float *out, spx_uint32_t out_len);
int resample_process_int_sse(SpeexResamplerState *st, const spx_int16_t *in,
        spx_uint32_t in_len, spx_int16_t *out, spx_uint32_t out_len);
int resample_process_il_sse(SpeexResamplerState *st, const float *in,
        spx_uint32_t in_len, float *out, spx_uint32_t out_len);
int resample_process_int_il_sse(SpeexResamplerState *st, const spx_int16_t *in,
        spx_uint32_t in_len, spx_int16_t *out, spx_uint32_t out_len);
#endif
#ifdef USE_RVV
SpeexResamplerState *resample_make_state_rvv(unsigned in_rate, unsigned out_rate, int quality);
SpeexResamplerState *resample_make_state_rvv_ch(unsigned in_rate, unsigned out_rate,
        int quality, unsigned channels);
int resample_process_rvv(SpeexResamplerState *st, const float *in,
        spx_uint32_t in_len, float *out, spx_uint32_t out_len);
int resample_process_int_rvv(SpeexResamplerState *st, const spx_int16_t *in,
        spx_uint32_t in_len, spx_int16_t *out, spx_uint32_t out_len);
int resample_process_il_rvv(SpeexResamplerState *st, const float *in,
        spx_uint32_t in_len, float *out, spx_uint32_t out_len);
int resample_process_int_il_rvv(SpeexResamplerState *st, const spx_int16_t *in,
        spx_uint32_t in_len, spx_int16_t *out, spx_uint32_t out_len);
#endif
#endif /* !DISABLE_FLOAT_API */

/* ------------- Test-input fill -------------
 * Reuses checkasm's randomizers. Inputs feed the real sinc-table dot product:
 * the table has ~unity gain (coefficients sum to ~1.0 / ~32768), so even
 * full-range int16 input keeps the int32 accumulator well within range. */
#include <checkasm/utils.h>

/* checkasm_init fills raw bytes; reinterpreted as float that yields the full
 * IEEE-754 exponent range plus NaN/Inf, which would make a long dot product
 * almost certainly NaN/Inf. Use checkasm's float-aware uniform helper and shift
 * to symmetric [-1, 1). */
static inline void resample_fill_float(float *buf, unsigned n)
{
    checkasm_randomize_rangef(buf, (int) n, 2.0f);
    for (unsigned i = 0; i < n; i++)
        buf[i] -= 1.0f;
}

static inline void resample_fill_input(spx_word16_t *buf, unsigned n)
{
#ifdef FIXED_POINT
    checkasm_init(buf, n * sizeof *buf);
#else
    resample_fill_float(buf, n);
#endif
}

/* ------------- Output comparison -------------
 * All four functions write spx_word16_t out[] (float in floating point). The
 * C and SIMD builds differ only in the inner-product summation, so:
 *  - float:  the buffer agrees to a relative tolerance, normalized by the
 *            buffer's peak |ref| (L-infinity). Normalizing by the peak rather
 *            than per-sample |ref| avoids the cancellation flakiness that a
 *            near-zero output sample would otherwise cause (see Higham, backward
 *            error of a dot product).
 *  - word16 (FIXED_POINT): C uses SATURATE32PSHR (clamps to +/-32767) while NEON
 *            uses sqrshrn/vqrshrn (reaches -32768), so at most one LSB diverges
 *            per output sample -> absolute 1-LSB tolerance. The RVV kernels
 *            instead keep C's wrapping int32 accumulate + scalar saturation,
 *            which is bit-exact -> they use resample_buffer_bitexact below so
 *            a tail-policy or combine regression can't hide in the LSB slack.
 *
 * Float tolerance pairings for the RVV kernels:
 *  - direct single + interpolate single: vfmacc is a true FMA and lane-split
 *    reduction reorders the sums; both sit far inside the 1e-4 tolerance.
 *  - direct double + interpolate double: vfwmacc forms exact f64 products where
 *    the C reference rounds each product to f32, so both pair with
 *    RESAMPLE_DOUBLE_FLOATMUL_REL_TOL (1e-6), not 1e-9. (SSE2 direct double
 *    multiplies in f32 like the C reference and keeps the tight 1e-9.) */
#define RESAMPLE_FLOAT_REL_TOL           1e-4f
#define RESAMPLE_DOUBLE_REL_TOL          1e-9   /* both sides multiply in float */
#define RESAMPLE_DOUBLE_FLOATMUL_REL_TOL 1e-6   /* C rounds products to f32; RVV's vfwmacc forms exact f64 */
#define RESAMPLE_WORD16_LSB_TOL          1

#ifndef FIXED_POINT
static inline double resample_peak(const spx_word16_t *a, unsigned n)
{
    double p = 0.0;
    for (unsigned i = 0; i < n; i++) {
        double v = fabs((double) a[i]);
        if (v > p) p = v;
    }
    return p;
}
#endif

/* Returns 1 if res agrees with ref over all n samples, else 0 (with a
 * diagnostic to stderr). rel_tol is ignored in FIXED_POINT. */
static inline int resample_buffer_within_tol(const spx_word16_t *ref,
        const spx_word16_t *res, unsigned n, double rel_tol)
{
#ifdef FIXED_POINT
    (void) rel_tol;
    for (unsigned i = 0; i < n; i++) {
        int32_t r = (int32_t) ref[i], s = (int32_t) res[i];
        int32_t diff = r > s ? r - s : s - r;
        if (diff > RESAMPLE_WORD16_LSB_TOL) {
            fprintf(stderr, "FAILED: sample %u ref=%" PRId32 " res=%" PRId32
                    " diff=%" PRId32 " (tol %d)\n",
                    i, r, s, diff, RESAMPLE_WORD16_LSB_TOL);
            return 0;
        }
    }
    return 1;
#else
    double peak = resample_peak(ref, n);
    for (unsigned i = 0; i < n; i++) {
        double diff = fabs((double) ref[i] - (double) res[i]);
        double rel  = peak > 0.0 ? diff / peak : diff;
        /* NaN-safe: rel > rel_tol would be false for NaN and wrongly pass. */
        if (!(rel <= rel_tol)) {
            fprintf(stderr, "FAILED: sample %u ref=%g res=%g diff=%g peak=%g "
                    "rel=%.2e (tol %g)\n",
                    i, (double) ref[i], (double) res[i], diff, peak, rel, rel_tol);
            return 0;
        }
    }
    return 1;
#endif
}

#ifdef FIXED_POINT
/* Exact-equality comparator for the RVV fixed-point branches: those kernels
 * reproduce the C kernels bit for bit (exact widening products, wrapping int32
 * accumulation -- associative under any lane split -- and the same scalar
 * MULT16_32_Q15 / SATURATE32PSHR combine), so any difference is a real bug. */
static inline int resample_buffer_bitexact(const spx_word16_t *ref,
        const spx_word16_t *res, unsigned n)
{
    for (unsigned i = 0; i < n; i++) {
        if (ref[i] != res[i]) {
            fprintf(stderr, "FAILED: sample %u ref=%" PRId32 " res=%" PRId32
                    " (bit-exact required)\n",
                    i, (int32_t) ref[i], (int32_t) res[i]);
            return 0;
        }
    }
    return 1;
}
#endif

/* ------------- Integration output comparison -------------
 * speex_resampler_process_float always writes plain float (even in FIXED_POINT,
 * where it converts), so the integration test cannot reuse the spx_word16_t
 * helper above. Compare relative to the buffer peak (L-infinity), like the
 * float path of resample_buffer_within_tol. */
#ifndef DISABLE_FLOAT_API
#define RESAMPLE_PROCESS_REL_TOL 1e-3   /* full-pipeline C-vs-SIMD, peak-relative */

static inline int resample_float_within_tol(const float *ref, const float *res,
        unsigned n, double rel_tol)
{
    double peak = 0.0;
    for (unsigned i = 0; i < n; i++) {
        double v = fabs((double) ref[i]);
        if (v > peak) peak = v;
    }
    for (unsigned i = 0; i < n; i++) {
        double diff = fabs((double) ref[i] - (double) res[i]);
        double rel  = peak > 0.0 ? diff / peak : diff;
        /* NaN-safe: rel > rel_tol would be false for NaN and wrongly pass. */
        if (!(rel <= rel_tol)) {
            fprintf(stderr, "FAILED: out[%u] ref=%g res=%g diff=%g peak=%g "
                    "rel=%.2e (tol %g)\n",
                    i, (double) ref[i], (double) res[i], diff, peak, rel, rel_tol);
            return 0;
        }
    }
    return 1;
}

/* int16 output path (process_int). Random full-range int16 input naturally
 * drives some samples into saturation, exercising WORD2INT. Two bounded effects
 * make C and SIMD diverge: the float resampling kernel itself differs by
 * ~1e-4 relative (different accumulation order), and the WORD2INT rounding
 * differs (NEON's fcvtas rounds to nearest-ties-away vs the C floor(.5+x)) by at
 * most 1 LSB. Both are small, so the int16 outputs differ by only a few LSB.
 *
 * We therefore gate on an absolute LSB bound rather than a peak-relative one: at
 * full scale the old 1e-3 peak-relative tolerance worked out to ~33 LSB of
 * slack -- loose enough to let a real WORD2INT/saturate regression through while
 * still passing. Measured worst case is 1 LSB on x86 over 30 seeds x 11
 * conversions (x86 shares the C WORD2INT, so only the kernel divergence shows);
 * aarch64 adds at most ~1 LSB more from fcvtas ties-away vs the C floor(.5+x).
 * RESAMPLE_PROCESS_INT_MAX_LSB sits just above that floor: tight enough to trip
 * on a rounding/saturation regression, with a little headroom for the unmeasured
 * NEON path. */
#define RESAMPLE_PROCESS_INT_MAX_LSB 4   /* absolute int16 LSB, full-pipeline C-vs-SIMD */

static inline void resample_fill_int16(spx_int16_t *buf, unsigned n)
{
    checkasm_init(buf, (size_t) n * sizeof *buf);
}

static inline int resample_int16_within_lsb(const spx_int16_t *ref,
        const spx_int16_t *res, unsigned n, int max_lsb)
{
    for (unsigned i = 0; i < n; i++) {
        int diff = abs((int) ref[i] - (int) res[i]);
        if (diff > max_lsb) {
            fprintf(stderr, "FAILED: out[%u] ref=%d res=%d diff=%d LSB (max %d)\n",
                    i, (int) ref[i], (int) res[i], diff, max_lsb);
            return 0;
        }
    }
    return 1;
}
#endif /* !DISABLE_FLOAT_API */

/* ------------- Per-routine test entry points -------------
 * Each is always defined; the body returns early when no SIMD build of that
 * routine exists for the current (mode, arch), so calling it is always safe. */
void test_resampler_basic_direct_single(void);
void test_resampler_basic_direct_double(void);
void test_resampler_basic_interpolate_single(void);
void test_resampler_basic_interpolate_double(void);

#endif /* SPEEXDSP_TESTS_CHECKASM_RESAMPLE_WRAP_H */
