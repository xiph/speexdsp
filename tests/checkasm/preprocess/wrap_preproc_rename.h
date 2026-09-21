/* Rename preprocess.c's public API (and the fftwrap/filterbank/mdf
 * symbols it references, which each TU stubs out) to a per-TU-unique
 * prefix so each wrap_preproc_*.c can #include the whole file without
 * link-time collisions -- same scheme as ../mdf/wrap_mdf_rename.h.
 * #define CKA_PREFIX and include this BEFORE wrap.h; it also names the
 * module source for ../wrap_impl.h. */
#include "../wrap_common.h"
#define CKA_MODULE_SRC "preprocess.c"

#define speex_preprocess_state_init      CKA_CAT(CKA_PREFIX, speex_preprocess_state_init)
#define speex_preprocess_state_destroy   CKA_CAT(CKA_PREFIX, speex_preprocess_state_destroy)
#define speex_preprocess                 CKA_CAT(CKA_PREFIX, speex_preprocess)
#define speex_preprocess_run             CKA_CAT(CKA_PREFIX, speex_preprocess_run)
#define speex_preprocess_estimate_update CKA_CAT(CKA_PREFIX, speex_preprocess_estimate_update)
#define speex_preprocess_ctl             CKA_CAT(CKA_PREFIX, speex_preprocess_ctl)

/* fftwrap/filterbank/mdf symbols preprocess.c links against; stubbed in
 * wrap_preproc_stubs.h. */
#define spx_fft_init                     CKA_CAT(CKA_PREFIX, spx_fft_init)
#define spx_fft_destroy                  CKA_CAT(CKA_PREFIX, spx_fft_destroy)
#define spx_fft                          CKA_CAT(CKA_PREFIX, spx_fft)
#define spx_ifft                         CKA_CAT(CKA_PREFIX, spx_ifft)
#define filterbank_new                   CKA_CAT(CKA_PREFIX, filterbank_new)
#define filterbank_destroy               CKA_CAT(CKA_PREFIX, filterbank_destroy)
#define filterbank_compute_bank32        CKA_CAT(CKA_PREFIX, filterbank_compute_bank32)
#define filterbank_compute_psd16         CKA_CAT(CKA_PREFIX, filterbank_compute_psd16)
#define filterbank_compute_bank          CKA_CAT(CKA_PREFIX, filterbank_compute_bank)
#define filterbank_compute_psd           CKA_CAT(CKA_PREFIX, filterbank_compute_psd)
#define speex_echo_get_residual          CKA_CAT(CKA_PREFIX, speex_echo_get_residual)
