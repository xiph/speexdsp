/* Rename kiss_fft.c's public API to a per-TU-unique prefix so each
 * wrap_fft_*.c can #include the whole file without link-time collisions --
 * same scheme as ../resample/wrap_resample_rename.h. #define CKA_PREFIX and
 * include this BEFORE wrap.h. */
#include "../wrap_common.h"
#define CKA_MODULE_SRC "kiss_fft.c"

#define kiss_fft_alloc  CKA_CAT(CKA_PREFIX, kiss_fft_alloc)
#define kiss_fft_stride CKA_CAT(CKA_PREFIX, kiss_fft_stride)
#define kiss_fft        CKA_CAT(CKA_PREFIX, kiss_fft)
