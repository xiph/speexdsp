/* SSE2 build of the double-precision resampler functions. Compiled only when
 * has_sse2 && !fixed-point. The native USE_SSE2 makes resample_sse.h define
 * OVERRIDE_INNER_PRODUCT_DOUBLE / OVERRIDE_INTERPOLATE_PRODUCT_DOUBLE, so the
 * double kernels use the SSE2 intrinsics. See wrap_resample_impl.h for the
 * EXPORT/HAVE_CONFIG_H/rename mechanics. */
#define CKA_PREFIX ckasse2_
#include "wrap_resample_rename.h"
#include "wrap.h"
#include "wrap_resample_impl.h"

int resampler_basic_direct_double_sse2(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len)
{
    st->last_sample[channel_index]   = 0;
    st->samp_frac_num[channel_index] = 0;
    return resampler_basic_direct_double(st, channel_index, in, in_len, out, out_len);
}

int resampler_basic_interpolate_double_sse2(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len)
{
    st->last_sample[channel_index]   = 0;
    st->samp_frac_num[channel_index] = 0;
    return resampler_basic_interpolate_double(st, channel_index, in, in_len, out, out_len);
}
