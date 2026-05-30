#ifndef SPEEXDSP_TESTS_CHECKASM_RESAMPLE_WRAP_H
#define SPEEXDSP_TESTS_CHECKASM_RESAMPLE_WRAP_H

#include <stdio.h>
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
 *            per output sample -> absolute 1-LSB tolerance. */
#define RESAMPLE_FLOAT_REL_TOL           1e-4f
#define RESAMPLE_DOUBLE_REL_TOL          1e-9   /* both sides multiply in float */
#define RESAMPLE_DOUBLE_FLOATMUL_REL_TOL 1e-6   /* C multiplies in double, SIMD in float */
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
        if (rel > rel_tol) {
            fprintf(stderr, "FAILED: sample %u ref=%g res=%g diff=%g peak=%g "
                    "rel=%.2e (tol %g)\n",
                    i, (double) ref[i], (double) res[i], diff, peak, rel, rel_tol);
            return 0;
        }
    }
    return 1;
#endif
}

/* ------------- Per-routine test entry points -------------
 * Each is always defined; the body returns early when no SIMD build of that
 * routine exists for the current (mode, arch), so calling it is always safe. */
void test_resampler_basic_direct_single(void);
void test_resampler_basic_direct_double(void);
void test_resampler_basic_interpolate_single(void);
void test_resampler_basic_interpolate_double(void);

#endif /* SPEEXDSP_TESTS_CHECKASM_RESAMPLE_WRAP_H */
