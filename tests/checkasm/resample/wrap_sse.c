/* Compiled only when has_sse && !FIXED_POINT.
 * SSE provides inner_product_single and interpolate_product_single in
 * floating-point mode (see libspeexdsp/resample_sse.h). */
#include "wrap.h"
#include "../../../libspeexdsp/resample_sse.h"

spx_word32_t inner_product_single_sse(const spx_word16_t *a,
                                      const spx_word16_t *b,
                                      unsigned int len)
{
    return inner_product_single(a, b, len);
}

spx_word32_t interpolate_product_single_sse(const spx_word16_t *a,
                                            const spx_word16_t *b,
                                            unsigned int len,
                                            spx_uint32_t oversample,
                                            const spx_word16_t *interp)
{
    /* SSE inline takes (float *frac) — non-const. The wrapper exposes the
     * canonical const signature; the cast is safe because the SSE body only
     * reads from interp. */
    return interpolate_product_single(a, b, len, oversample,
                                      (spx_word16_t *)interp);
}
