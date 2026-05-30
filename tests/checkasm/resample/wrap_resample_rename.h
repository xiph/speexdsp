/* Rename resample.c's public speex_resampler_* API to a per-TU-unique prefix.
 *
 * Each wrap_resample_*.c #includes the whole resample.c so it can build a C or
 * SIMD variant of the four static resampler_basic_* functions. resample.c's
 * internal forward-references to the public API (e.g. speex_resampler_init_frac
 * calls speex_resampler_set_rate_frac, which is defined later) resolve through
 * the external prototypes in speex_resampler.h, so the public functions must
 * keep external linkage -- a `static` redefinition would both clash with those
 * prototypes and break the forward calls. Instead we give each TU its own unique
 * prefix so the (external) public symbols don't collide at link time. Everything
 * else in resample.c is already file-scope static.
 *
 * The caller must #define CKA_PREFIX to a unique token and include this BEFORE
 * wrap.h, so speex_resampler.h's prototypes are renamed to match the (renamed)
 * definitions and forward-call sites in resample.c. */
#ifndef CKA_PREFIX
#  error "define CKA_PREFIX (a unique token) before including wrap_resample_rename.h"
#endif

#define CKA_CAT2(a, b) a ## b
#define CKA_CAT(a, b)  CKA_CAT2(a, b)

#define speex_resampler_init                      CKA_CAT(CKA_PREFIX, speex_resampler_init)
#define speex_resampler_init_frac                 CKA_CAT(CKA_PREFIX, speex_resampler_init_frac)
#define speex_resampler_destroy                   CKA_CAT(CKA_PREFIX, speex_resampler_destroy)
#define speex_resampler_process_int               CKA_CAT(CKA_PREFIX, speex_resampler_process_int)
#define speex_resampler_process_float             CKA_CAT(CKA_PREFIX, speex_resampler_process_float)
#define speex_resampler_process_interleaved_float CKA_CAT(CKA_PREFIX, speex_resampler_process_interleaved_float)
#define speex_resampler_process_interleaved_int   CKA_CAT(CKA_PREFIX, speex_resampler_process_interleaved_int)
#define speex_resampler_set_rate                  CKA_CAT(CKA_PREFIX, speex_resampler_set_rate)
#define speex_resampler_get_rate                  CKA_CAT(CKA_PREFIX, speex_resampler_get_rate)
#define speex_resampler_set_rate_frac             CKA_CAT(CKA_PREFIX, speex_resampler_set_rate_frac)
#define speex_resampler_get_ratio                 CKA_CAT(CKA_PREFIX, speex_resampler_get_ratio)
#define speex_resampler_set_quality               CKA_CAT(CKA_PREFIX, speex_resampler_set_quality)
#define speex_resampler_get_quality               CKA_CAT(CKA_PREFIX, speex_resampler_get_quality)
#define speex_resampler_set_input_stride          CKA_CAT(CKA_PREFIX, speex_resampler_set_input_stride)
#define speex_resampler_get_input_stride          CKA_CAT(CKA_PREFIX, speex_resampler_get_input_stride)
#define speex_resampler_set_output_stride         CKA_CAT(CKA_PREFIX, speex_resampler_set_output_stride)
#define speex_resampler_get_output_stride         CKA_CAT(CKA_PREFIX, speex_resampler_get_output_stride)
#define speex_resampler_get_input_latency         CKA_CAT(CKA_PREFIX, speex_resampler_get_input_latency)
#define speex_resampler_get_output_latency        CKA_CAT(CKA_PREFIX, speex_resampler_get_output_latency)
#define speex_resampler_skip_zeros                CKA_CAT(CKA_PREFIX, speex_resampler_skip_zeros)
#define speex_resampler_reset_mem                 CKA_CAT(CKA_PREFIX, speex_resampler_reset_mem)
#define speex_resampler_strerror                  CKA_CAT(CKA_PREFIX, speex_resampler_strerror)
