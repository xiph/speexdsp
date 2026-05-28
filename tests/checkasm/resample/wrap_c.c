/* Always compiled. Provides the *_c reference implementations against which
 * the SIMD wrappers are checked. We include only resample.h here -- crucially
 * NOT resample_sse.h or resample_neon.h -- so OVERRIDE_INNER_PRODUCT_SINGLE
 * (and friends) stay undefined and the generic C fallbacks at resample.h:43
 * and below are what end up compiled into this TU. The SIMD wrappers live in
 * their own TUs (wrap_sse.c, wrap_sse2.c, wrap_neon.c), each pulling in the
 * matching resample_*.h, so the four variants coexist in the same binary
 * without symbol collision (each routine is static inline, file-scoped). */
#include "wrap.h"
#include "../../../libspeexdsp/resample.h"

spx_word32_t inner_product_single_c(const spx_word16_t *a,
                                    const spx_word16_t *b,
                                    unsigned int len)
{
    return inner_product_single(a, b, len);
}

#ifndef FIXED_POINT
double inner_product_double_c(const spx_word16_t *a,
                              const spx_word16_t *b,
                              unsigned int len)
{
    return inner_product_double(a, b, len);
}
#endif

spx_word32_t interpolate_product_single_c(const spx_word16_t *a,
                                          const spx_word16_t *b,
                                          unsigned int len,
                                          spx_uint32_t oversample,
                                          const spx_word16_t *interp)
{
    return interpolate_product_single(a, b, len, oversample, interp);
}

#ifndef FIXED_POINT
double interpolate_product_double_c(const spx_word16_t *a,
                                    const spx_word16_t *b,
                                    unsigned int len,
                                    spx_uint32_t oversample,
                                    const spx_word16_t *interp)
{
    return interpolate_product_double(a, b, len, oversample, interp);
}
#endif
