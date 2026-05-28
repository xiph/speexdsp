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
