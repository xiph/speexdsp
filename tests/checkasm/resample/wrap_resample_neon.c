/* NEON build of the resampler functions. Compiled only when has_neon. Keeps the
 * native USE_NEON so resample.c pulls in resample_neon.h. Today that header only
 * overrides inner_product_single (OVERRIDE_INNER_PRODUCT_SINGLE), so only
 * resampler_basic_direct_single actually differs from the C reference under
 * NEON; the other three would be byte-identical to C and are not exposed.
 * See wrap_resample_impl.h for the EXPORT/HAVE_CONFIG_H/rename mechanics. */
#define CKA_PREFIX ckaneon_
#include "wrap_resample_rename.h"
#include "wrap.h"

#ifdef USE_NEON
#  if !defined(__aarch64__) && !defined(__ARM_NEON)
#    error "wrap_resample_neon.c built without aarch64 / ARM NEON ISA; check toolchain flags"
#  endif
#endif

#include "wrap_resample_impl.h"

#ifdef HAVE_NEON_DIRECT_SINGLE
int resampler_basic_direct_single_neon(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len)
{
    st->last_sample[channel_index]   = 0;
    st->samp_frac_num[channel_index] = 0;
    return resampler_basic_direct_single(st, channel_index, in, in_len, out, out_len);
}
#endif

/* When resample_neon.h grows OVERRIDE_INTERPOLATE_PRODUCT_SINGLE / *_DOUBLE,
 * flip the matching HAVE_NEON_* macro on in wrap.h and add the wrapper(s) here,
 * mirroring resampler_basic_direct_single_neon above. */
