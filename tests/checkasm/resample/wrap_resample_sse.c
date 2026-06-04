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

/* ------------- Integration: full-pipeline wrapper (SSE/SSE2 kernels) ----------
 * On x86-64 this TU keeps both USE_SSE and USE_SSE2, so the state's resampler_ptr
 * lands on the SSE single / SSE2 double kernels as appropriate -- the full
 * library pipeline. See wrap.h for the by-value / zero-state contract. */
#ifndef DISABLE_FLOAT_API
SpeexResamplerState *resample_make_state_sse(unsigned in_rate, unsigned out_rate, int quality)
{
    int err = RESAMPLER_ERR_SUCCESS;
    return speex_resampler_init(1, in_rate, out_rate, quality, &err);
}

int resample_process_sse(SpeexResamplerState *st, const float *in,
        spx_uint32_t in_len, float *out, spx_uint32_t out_len)
{
    spx_uint32_t il = in_len, ol = out_len;
    memset(st->mem, 0, (size_t) st->mem_alloc_size * st->nb_channels * sizeof(spx_word16_t));
    st->last_sample[0]   = 0;
    st->samp_frac_num[0] = 0;
    st->magic_samples[0] = 0;
    st->started          = 0;
    speex_resampler_process_float(st, 0, in, &il, out, &ol);
    return (int) ol;
}

int resample_process_int_sse(SpeexResamplerState *st, const spx_int16_t *in,
        spx_uint32_t in_len, spx_int16_t *out, spx_uint32_t out_len)
{
    spx_uint32_t il = in_len, ol = out_len;
    memset(st->mem, 0, (size_t) st->mem_alloc_size * st->nb_channels * sizeof(spx_word16_t));
    st->last_sample[0]   = 0;
    st->samp_frac_num[0] = 0;
    st->magic_samples[0] = 0;
    st->started          = 0;
    speex_resampler_process_int(st, 0, in, &il, out, &ol);
    return (int) ol;
}
#endif
