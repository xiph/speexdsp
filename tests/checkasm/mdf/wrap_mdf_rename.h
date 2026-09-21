/* Rename mdf.c's public API (and the fftwrap symbols it references, which
 * each TU stubs out) to a per-TU-unique prefix so each wrap_mdf_*.c can
 * #include the whole file without link-time collisions -- same scheme as
 * ../fft/wrap_fft_rename.h. #define CKA_PREFIX and include this BEFORE
 * wrap.h; it also names the module source for ../wrap_impl.h. */
#include "../wrap_common.h"
#define CKA_MODULE_SRC "mdf.c"

#define speex_echo_state_init      CKA_CAT(CKA_PREFIX, speex_echo_state_init)
#define speex_echo_state_init_mc   CKA_CAT(CKA_PREFIX, speex_echo_state_init_mc)
#define speex_echo_state_reset     CKA_CAT(CKA_PREFIX, speex_echo_state_reset)
#define speex_echo_state_destroy   CKA_CAT(CKA_PREFIX, speex_echo_state_destroy)
#define speex_echo_capture         CKA_CAT(CKA_PREFIX, speex_echo_capture)
#define speex_echo_playback        CKA_CAT(CKA_PREFIX, speex_echo_playback)
#define speex_echo_cancel          CKA_CAT(CKA_PREFIX, speex_echo_cancel)
#define speex_echo_cancellation    CKA_CAT(CKA_PREFIX, speex_echo_cancellation)
#define speex_echo_ctl             CKA_CAT(CKA_PREFIX, speex_echo_ctl)
#define speex_echo_get_residual    CKA_CAT(CKA_PREFIX, speex_echo_get_residual)

/* fftwrap symbols mdf.c links against; stubbed in wrap_mdf_stubs.h. */
#define spx_fft_init               CKA_CAT(CKA_PREFIX, spx_fft_init)
#define spx_fft_destroy            CKA_CAT(CKA_PREFIX, spx_fft_destroy)
#define spx_fft                    CKA_CAT(CKA_PREFIX, spx_fft)
#define spx_ifft                   CKA_CAT(CKA_PREFIX, spx_ifft)
