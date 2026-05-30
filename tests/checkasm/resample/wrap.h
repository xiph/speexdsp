#ifndef SPEEXDSP_TESTS_CHECKASM_RESAMPLE_WRAP_H
#define SPEEXDSP_TESTS_CHECKASM_RESAMPLE_WRAP_H

#include <stdio.h>
#include <math.h>
#include <inttypes.h>
#include "config.h"
#include "arch.h"

/* Per-routine NEON-availability gates. wrap_neon.c only provides a wrapper
 * when libspeexdsp's resample_neon.h defines the matching OVERRIDE_* macro.
 * Today, only inner_product_single has a NEON impl. Flip these on when
 * resample_neon.h grows new OVERRIDE_* macros. */
#ifdef USE_NEON
#  define HAVE_NEON_INNER_PRODUCT_SINGLE       1
/* #  define HAVE_NEON_INNER_PRODUCT_DOUBLE       1 */
/* #  define HAVE_NEON_INTERPOLATE_PRODUCT_SINGLE 1 */
/* #  define HAVE_NEON_INTERPOLATE_PRODUCT_DOUBLE 1 */
#endif

/* ------------- inner_product_single -------------
 * Available in both FIXED_POINT and FLOATING_POINT.
 * SIMD: NEON (both modes), SSE (FP only).
 */
spx_word32_t inner_product_single_c(const spx_word16_t *a,
                                    const spx_word16_t *b,
                                    unsigned int len);
#ifdef HAVE_NEON_INNER_PRODUCT_SINGLE
spx_word32_t inner_product_single_neon(const spx_word16_t *a,
                                       const spx_word16_t *b,
                                       unsigned int len);
#endif
#if defined(USE_SSE) && !defined(FIXED_POINT)
spx_word32_t inner_product_single_sse(const spx_word16_t *a,
                                      const spx_word16_t *b,
                                      unsigned int len);
#endif

/* ------------- inner_product_double -------------
 * FLOATING_POINT only.
 * SIMD: SSE2 (not SSE).
 * TODO: NEON
 */
#ifndef FIXED_POINT
double inner_product_double_c(const spx_word16_t *a,
                              const spx_word16_t *b,
                              unsigned int len);
#ifdef USE_SSE2
double inner_product_double_sse2(const spx_word16_t *a,
                                 const spx_word16_t *b,
                                 unsigned int len);
#endif
#ifdef HAVE_NEON_INNER_PRODUCT_DOUBLE
double inner_product_double_neon(const spx_word16_t *a,
                                 const spx_word16_t *b,
                                 unsigned int len);
#endif
#endif

/* ------------- interpolate_product_single -------------
 * Available in both FIXED_POINT and FLOATING_POINT.
 * SIMD: SSE (FP only).
 * TODO: NEON
 */
spx_word32_t interpolate_product_single_c(const spx_word16_t *a,
                                          const spx_word16_t *b,
                                          unsigned int len,
                                          spx_uint32_t oversample,
                                          const spx_word16_t *interp);
#if defined(USE_SSE) && !defined(FIXED_POINT)
spx_word32_t interpolate_product_single_sse(const spx_word16_t *a,
                                            const spx_word16_t *b,
                                            unsigned int len,
                                            spx_uint32_t oversample,
                                            const spx_word16_t *interp);
#endif
#ifdef HAVE_NEON_INTERPOLATE_PRODUCT_SINGLE
spx_word32_t interpolate_product_single_neon(const spx_word16_t *a,
                                             const spx_word16_t *b,
                                             unsigned int len,
                                             spx_uint32_t oversample,
                                             const spx_word16_t *interp);
#endif

/* ------------- interpolate_product_double -------------
 * FLOATING_POINT only.
 * SIMD: SSE2 (not SSE).
 * TODO: NEON
 */
#ifndef FIXED_POINT
double interpolate_product_double_c(const spx_word16_t *a,
                                    const spx_word16_t *b,
                                    unsigned int len,
                                    spx_uint32_t oversample,
                                    const spx_word16_t *interp);
#ifdef USE_SSE2
double interpolate_product_double_sse2(const spx_word16_t *a,
                                       const spx_word16_t *b,
                                       unsigned int len,
                                       spx_uint32_t oversample,
                                       const spx_word16_t *interp);
#endif
#ifdef HAVE_NEON_INTERPOLATE_PRODUCT_DOUBLE
double interpolate_product_double_neon(const spx_word16_t *a,
                                       const spx_word16_t *b,
                                       unsigned int len,
                                       spx_uint32_t oversample,
                                       const spx_word16_t *interp);
#endif
#endif

/* ------------- Test-input buffer fillers -------------
 * Single-buffer form because some tests (interpolate_product_*) have
 * differently-sized a and b. */
#include <checkasm/utils.h>

#ifdef FIXED_POINT
/* Random bytes are always valid int16, but a dot product of int16^2 over
 * hundreds of terms overflows int32. Shift right by 6 to cap |a*b| at 2^18,
 * keeping comfortable int32 headroom. The bound also makes the two NEON
 * variants bit-identical to C: aarch64 (saddlv -> sqxtn -> sqrshrn #15) and
 * aarch32 (vaddl.s32 + vadd.s64 -> vqmovn.s64 -> vqrshrn.s32 #15) perform
 * the same width/saturation sequence on inputs in this range. */
static inline void resample_fill_word16(spx_word16_t *buf, unsigned n)
{
    checkasm_init(buf, n * sizeof *buf);
    for (unsigned i = 0; i < n; i++)
        buf[i] >>= 6;
}
#endif

/* checkasm_init fills raw bytes; reinterpreted as float that yields the full
 * IEEE-754 exponent range plus NaN/Inf, so a long dot product is almost
 * certainly NaN/Inf and the comparison trivially "passes" (NaN > x is false).
 * Use checkasm's float-aware uniform helper and shift to symmetric [-1, 1). */
static inline void resample_fill_float(float *buf, unsigned n)
{
    checkasm_randomize_rangef(buf, (int) n, 2.0f);
    for (unsigned i = 0; i < n; i++)
        buf[i] -= 1.0f;
}

/* ------------- Result tolerance checks -------------
 * Per-precision relative tolerances shared by all resample tests. Each helper
 * returns 1 if the SIMD result agrees with the C reference within tolerance,
 * 0 otherwise (with a diagnostic emitted to stderr). Caller invokes
 * checkasm_fail() on a 0 return.
 *
 * The error is normalized by the accumulation scale S = sum |a_j * b_j| (the
 * caller passes it in, e.g. via resample_abs_dot_scale below) rather than by
 * |ref|. This is the standard backward-error denominator for a dot product
 * (Higham, "Accuracy and Stability of Numerical Algorithms", section 3.1):
 * two correct float dot products differing only in summation order satisfy
 * |res - ref| <~ N*eps*S regardless of the signed result. Normalizing by
 * |ref| instead breaks down under cancellation -- when the random inputs
 * happen to sum near zero, a perfectly normal absolute rounding difference
 * divided by a tiny |ref| spuriously trips the tolerance (a flaky failure
 * dependent on the checkasm seed). S stays ~N/4 even when ref ~ 0, so the
 * normalized error tracks the actual conditioning of the computation.
 *
 * Tolerances (now applied to diff/S):
 *   float:  1e-4 -- generous over N*eps_float (~1.5e-5 at N=256).
 *   double: 1e-9 -- ~10^4 x headroom over N*eps_double (~2.8e-14 at N=256).
 *   word16: 1 LSB absolute -- C uses SATURATE32PSHR (clamps to +/-32767),
 *           NEON uses sqrshrn/vqrshrn (reaches -32768). At most one LSB
 *           divergence at the negative endpoint. */
#define RESAMPLE_FLOAT_REL_TOL  1e-4f
#define RESAMPLE_DOUBLE_REL_TOL 1e-9
#define RESAMPLE_WORD16_LSB_TOL 1

/* Accumulation scale S = sum_j |a_j * b_j| for a length-len dot product, used
 * as the denominator in the float/double tolerance checks. Accumulated in
 * double so the scale itself is effectively exact relative to the float work
 * being judged. */
static inline double resample_abs_dot_scale(const float *a, const float *b, unsigned len)
{
    double s = 0.0;
    for (unsigned j = 0; j < len; j++)
        s += fabs((double)a[j] * (double)b[j]);
    return s;
}

/* Accumulation scale for the interpolate routines, whose result is
 * sum_k interp[k] * (sum_j a[j] * b[j*oversample + k]) (resample.h:82-118).
 * The matching abs-conditioning scale is
 * sum_k |interp[k]| * sum_j |a[j] * b[j*oversample + k]|. */
static inline double resample_abs_interp_scale(const float *a, const float *b,
        unsigned len, unsigned oversample, const float *interp)
{
    double accum[4] = {0, 0, 0, 0};
    for (unsigned j = 0; j < len; j++)
        for (int k = 0; k < 4; k++)
            accum[k] += fabs((double)a[j] * (double)b[j * oversample + k]);
    double s = 0.0;
    for (int k = 0; k < 4; k++)
        s += fabs((double)interp[k]) * accum[k];
    return s;
}

static inline int is_resample_result_within_tolerance_float(float ref, float res, double scale, unsigned len)
{
    double diff = fabs((double)ref - (double)res);
    /* scale <= 0 means an all-zero input; fall back to an absolute check so a
     * genuine divergence is still caught and we never divide by zero. */
    double rel  = scale > 0.0 ? diff / scale : diff;
    if (rel > RESAMPLE_FLOAT_REL_TOL) {
        fprintf(stderr,
            "FAILED: len=%u ref=%g res=%g diff=%g scale=%g rel=%.2e (tol %g)\n",
            len, ref, res, diff, scale, rel, (double)RESAMPLE_FLOAT_REL_TOL);
        return 0;
    }
    return 1;
}

static inline int is_resample_result_within_tolerance_double(double ref, double res, double scale, unsigned len)
{
    double diff = fabs(ref - res);
    double rel  = scale > 0.0 ? diff / scale : diff;
    if (rel > RESAMPLE_DOUBLE_REL_TOL) {
        fprintf(stderr,
            "FAILED: len=%u ref=%g res=%g diff=%g scale=%g rel=%.2e (tol %g)\n",
            len, ref, res, diff, scale, rel, RESAMPLE_DOUBLE_REL_TOL);
        return 0;
    }
    return 1;
}

#ifdef FIXED_POINT
static inline int is_resample_result_within_tolerance_word16(int32_t ref, int32_t res, unsigned len)
{
    int32_t diff = ref > res ? ref - res : res - ref;
    if (diff > RESAMPLE_WORD16_LSB_TOL) {
        fprintf(stderr,
            "FAILED: len=%u ref=%" PRId32 " res=%" PRId32 " diff=%" PRId32 " (tol %d)\n",
            len, ref, res, diff, RESAMPLE_WORD16_LSB_TOL);
        return 0;
    }
    return 1;
}
#endif

/* Per-routine test entry points. Each is always defined; the body returns
 * early when no SIMD implementation exists for the current (mode, arch)
 * build, so calling them is always safe. */
void test_inner_product_single(void);
void test_inner_product_double(void);
void test_interpolate_product_single(void);
void test_interpolate_product_double(void);

#endif /* SPEEXDSP_TESTS_CHECKASM_RESAMPLE_WRAP_H */
