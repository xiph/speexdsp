#include "wrap.h"
#include "../../../libspeexdsp/resample_neon.h"

#ifdef HAVE_NEON_INNER_PRODUCT_SINGLE
spx_word32_t inner_product_single_neon(const spx_word16_t *a,
                                       const spx_word16_t *b,
                                       unsigned int len)
{
    return inner_product_single(a, b, len);
}
#endif

/* TODO: when resample_neon.h grows OVERRIDE_INNER_PRODUCT_DOUBLE, flip
 * HAVE_NEON_INNER_PRODUCT_DOUBLE on in wrap.h and uncomment this block. */
#ifdef HAVE_NEON_INNER_PRODUCT_DOUBLE
#ifndef FIXED_POINT
double inner_product_double_neon(const spx_word16_t *a,
                                 const spx_word16_t *b,
                                 unsigned int len)
{
    return inner_product_double(a, b, len);
}
#endif
#endif

/* TODO: when resample_neon.h grows OVERRIDE_INTERPOLATE_PRODUCT_SINGLE, flip
 * HAVE_NEON_INTERPOLATE_PRODUCT_SINGLE on in wrap.h and uncomment this block. */
#ifdef HAVE_NEON_INTERPOLATE_PRODUCT_SINGLE
spx_word32_t interpolate_product_single_neon(const spx_word16_t *a,
                                             const spx_word16_t *b,
                                             unsigned int len,
                                             spx_uint32_t oversample,
                                             const spx_word16_t *interp)
{
    /* If the NEON inline takes a non-const last arg (matching the SSE form),
     * cast through. If it takes const, the cast is a no-op. */
    return interpolate_product_single(a, b, len, oversample,
                                      (spx_word16_t *)interp);
}
#endif

/* TODO: when resample_neon.h grows OVERRIDE_INTERPOLATE_PRODUCT_DOUBLE, flip
 * HAVE_NEON_INTERPOLATE_PRODUCT_DOUBLE on in wrap.h and uncomment this block. */
#ifdef HAVE_NEON_INTERPOLATE_PRODUCT_DOUBLE
#ifndef FIXED_POINT
double interpolate_product_double_neon(const spx_word16_t *a,
                                       const spx_word16_t *b,
                                       unsigned int len,
                                       spx_uint32_t oversample,
                                       const spx_word16_t *interp)
{
    return interpolate_product_double(a, b, len, oversample,
                                      (spx_word16_t *)interp);
}
#endif
#endif
