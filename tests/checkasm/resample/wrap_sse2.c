/* Compiled only when has_sse2 && !FIXED_POINT.
 * SSE2 additionally provides inner_product_double and interpolate_product_double
 * in floating-point mode (see libspeexdsp/resample_sse.h, USE_SSE2 block). */
#include "wrap.h"
#include "../../../libspeexdsp/resample_sse.h"

double inner_product_double_sse2(const spx_word16_t *a,
                                 const spx_word16_t *b,
                                 unsigned int len)
{
    return inner_product_double(a, b, len);
}

double interpolate_product_double_sse2(const spx_word16_t *a,
                                       const spx_word16_t *b,
                                       unsigned int len,
                                       spx_uint32_t oversample,
                                       const spx_word16_t *interp)
{
    return interpolate_product_double(a, b, len, oversample,
                                      (spx_word16_t *)interp);
}
