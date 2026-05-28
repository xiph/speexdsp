#ifndef SPEEXDSP_TESTS_CHECKASM_RESAMPLE_WRAP_H
#define SPEEXDSP_TESTS_CHECKASM_RESAMPLE_WRAP_H

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

/* Per-routine test entry points. Each is always defined; the body is empty
 * when no SIMD implementation exists for the current (mode, arch) build. */
void test_inner_product_single(void);
void test_inner_product_double(void);
void test_interpolate_product_single(void);
void test_interpolate_product_double(void);

#endif /* SPEEXDSP_TESTS_CHECKASM_RESAMPLE_WRAP_H */
