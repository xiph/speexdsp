#include "config.h"
#include "../../../libspeexdsp/resample.h"

#ifdef FIXED_POINT
int32_t inner_product_single_c(const int16_t *a, const int16_t *b, unsigned int len) {
    return inner_product_single(a, b, len);
}
#else
float inner_product_single_c(const float *a, const float *b, unsigned int len) {
    return inner_product_single(a, b, len);
}
#endif
