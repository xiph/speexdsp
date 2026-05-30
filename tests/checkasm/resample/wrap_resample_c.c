/* Pure-C reference build of the four resampler_basic_* functions, plus the
 * shared state helpers. This is the only TU that owns SpeexResamplerState
 * construction/inspection (it has both the four static function addresses and
 * the public init API in scope after #including resample.c).
 *
 * We #undef the native USE_SSE/USE_SSE2/USE_NEON before pulling in resample.c so
 * its #ifdef USE_SSE/USE_NEON skip the SIMD headers and the inner-product
 * kernels stay the generic C fallbacks. wrap_resample_impl.h then handles the
 * EXPORT/HAVE_CONFIG_H/rename mechanics and includes resample.c. */
#define CKA_PREFIX ckac_
#include "wrap_resample_rename.h"
#include "wrap.h"

#undef USE_SSE
#undef USE_SSE2
#undef USE_NEON

#include "wrap_resample_impl.h"

/* ------------- State helpers ------------- */

SpeexResamplerState *resample_make_state(unsigned in_rate, unsigned out_rate, int quality)
{
    int err = RESAMPLER_ERR_SUCCESS;
    return speex_resampler_init(1, in_rate, out_rate, quality, &err);
}

void resample_destroy_state(SpeexResamplerState *st)
{
    speex_resampler_destroy(st);
}

enum resample_kind resample_kind(const SpeexResamplerState *st)
{
    if (st->resampler_ptr == resampler_basic_direct_single)
        return RESAMPLE_KIND_DIRECT_SINGLE;
    if (st->resampler_ptr == resampler_basic_interpolate_single)
        return RESAMPLE_KIND_INTERPOLATE_SINGLE;
#ifndef FIXED_POINT
    if (st->resampler_ptr == resampler_basic_direct_double)
        return RESAMPLE_KIND_DIRECT_DOUBLE;
    if (st->resampler_ptr == resampler_basic_interpolate_double)
        return RESAMPLE_KIND_INTERPOLATE_DOUBLE;
#endif
    return RESAMPLE_KIND_OTHER;
}

unsigned resample_filt_len(const SpeexResamplerState *st)   { return st->filt_len; }
unsigned resample_den_rate(const SpeexResamplerState *st)   { return st->den_rate; }
unsigned resample_oversample(const SpeexResamplerState *st) { return st->oversample; }

/* ------------- C reference wrappers -------------
 * Reset the per-channel cursor so every call (notably checkasm_bench_new's
 * loop) re-does identical full work. */
int resampler_basic_direct_single_c(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len)
{
    st->last_sample[channel_index]   = 0;
    st->samp_frac_num[channel_index] = 0;
    return resampler_basic_direct_single(st, channel_index, in, in_len, out, out_len);
}

int resampler_basic_interpolate_single_c(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len)
{
    st->last_sample[channel_index]   = 0;
    st->samp_frac_num[channel_index] = 0;
    return resampler_basic_interpolate_single(st, channel_index, in, in_len, out, out_len);
}

#ifndef FIXED_POINT
int resampler_basic_direct_double_c(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len)
{
    st->last_sample[channel_index]   = 0;
    st->samp_frac_num[channel_index] = 0;
    return resampler_basic_direct_double(st, channel_index, in, in_len, out, out_len);
}

int resampler_basic_interpolate_double_c(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len)
{
    st->last_sample[channel_index]   = 0;
    st->samp_frac_num[channel_index] = 0;
    return resampler_basic_interpolate_double(st, channel_index, in, in_len, out, out_len);
}
#endif
