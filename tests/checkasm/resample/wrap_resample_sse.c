/* SSE build of the single-precision resampler functions. Compiled only when
 * has_sse && !fixed-point. Keeps the native USE_SSE so resample.c pulls in
 * resample_sse.h, overriding inner_product_single / interpolate_product_single.
 * See wrap_resample_impl.h for the EXPORT/HAVE_CONFIG_H/rename mechanics. */
#define CKA_PREFIX ckasse_
#include "wrap_resample_rename.h"
#include "wrap.h"
#include "wrap_resample_impl.h"

int resampler_basic_direct_single_sse(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len)
{
    st->last_sample[channel_index]   = 0;
    st->samp_frac_num[channel_index] = 0;
    return resampler_basic_direct_single(st, channel_index, in, in_len, out, out_len);
}

int resampler_basic_interpolate_single_sse(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len)
{
    st->last_sample[channel_index]   = 0;
    st->samp_frac_num[channel_index] = 0;
    return resampler_basic_interpolate_single(st, channel_index, in, in_len, out, out_len);
}
