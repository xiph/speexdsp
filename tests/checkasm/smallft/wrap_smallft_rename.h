/* Rename smallft.c's public API to a per-TU-unique prefix so each
 * wrap_smallft_*.c can #include the whole file without link-time
 * collisions -- same scheme as ../fft/wrap_fft_rename.h. #define
 * CKA_PREFIX and include this BEFORE wrap.h. */
#ifndef CKA_PREFIX
#  error "define CKA_PREFIX (a unique token) before including wrap_smallft_rename.h"
#endif

#define CKA_CAT2(a, b) a ## b
#define CKA_CAT(a, b)  CKA_CAT2(a, b)

#define spx_drft_forward  CKA_CAT(CKA_PREFIX, spx_drft_forward)
#define spx_drft_backward CKA_CAT(CKA_PREFIX, spx_drft_backward)
#define spx_drft_init     CKA_CAT(CKA_PREFIX, spx_drft_init)
#define spx_drft_clear    CKA_CAT(CKA_PREFIX, spx_drft_clear)
