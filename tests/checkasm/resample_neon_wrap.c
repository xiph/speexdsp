#include "config.h"
#define OVERRIDE_INNER_PRODUCT_SINGLE
#include "../../../libspeexdsp/resample_neon.h"

// Explicit declaration needed because the header defines it as static inline
#ifdef FIXED_POINT
int32_t inner_product_single(const int16_t *a, const int16_t *b, unsigned int len);
#else
float inner_product_single(const float *a, const float *b, unsigned int len);
#endif

#ifdef FIXED_POINT
int32_t inner_product_single_neon(const int16_t *a, const int16_t *b, unsigned int len) {
    return inner_product_single(a, b, len);
}
#else
float inner_product_single_neon(const float *a, const float *b, unsigned int len) {
    return inner_product_single(a, b, len);
}
#endif
