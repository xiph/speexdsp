/* RVV build of the resampler kernels for checkasm. Base-ISA C; the V
 * instructions live only in resample_rvv_asm.S (linked in alongside). Keeps the
 * native USE_RVV so resample.c pulls in resample_rvv.h.
 *
 * The library dispatches at runtime via spx_rvv_enabled; the harness instead
 * sets RESAMPLE_RVV_FORCE_ON (=> SPX_RVV_ON == 1) so checkasm deterministically
 * tests the asm kernel, not whatever the build host's getauxval reports.
 * See wrap_resample_impl.h for the rename mechanics. */
#define CKA_PREFIX ckarvv_
#include "wrap_resample_rename.h"
#include "wrap.h"

#ifdef USE_RVV
#  define RESAMPLE_RVV_FORCE_ON
#endif

#include "wrap_resample_impl.h"

#ifdef USE_RVV

int resampler_basic_direct_single_rvv(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len)
{
    st->last_sample[channel_index]   = 0;
    st->samp_frac_num[channel_index] = 0;
    return resampler_basic_direct_single(st, channel_index, in, in_len, out, out_len);
}

int resampler_basic_interpolate_single_rvv(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len)
{
    st->last_sample[channel_index]   = 0;
    st->samp_frac_num[channel_index] = 0;
    return resampler_basic_interpolate_single(st, channel_index, in, in_len, out, out_len);
}

#ifndef FIXED_POINT
int resampler_basic_direct_double_rvv(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len)
{
    st->last_sample[channel_index]   = 0;
    st->samp_frac_num[channel_index] = 0;
    return resampler_basic_direct_double(st, channel_index, in, in_len, out, out_len);
}

int resampler_basic_interpolate_double_rvv(SpeexResamplerState *st, spx_uint32_t channel_index,
        const spx_word16_t *in, spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len)
{
    st->last_sample[channel_index]   = 0;
    st->samp_frac_num[channel_index] = 0;
    return resampler_basic_interpolate_double(st, channel_index, in, in_len, out, out_len);
}
#endif

/* Full-pipeline wrappers: update_filter in this TU points resampler_ptr at the
 * RVV kernels, so the public process functions run the RVV path. */
#ifndef DISABLE_FLOAT_API
SpeexResamplerState *resample_make_state_rvv(unsigned in_rate, unsigned out_rate, int quality)
{
    int err = RESAMPLER_ERR_SUCCESS;
    return speex_resampler_init(1, in_rate, out_rate, quality, &err);
}

SpeexResamplerState *resample_make_state_rvv_ch(unsigned in_rate, unsigned out_rate,
        int quality, unsigned channels)
{
    int err = RESAMPLER_ERR_SUCCESS;
    return speex_resampler_init(channels, in_rate, out_rate, quality, &err);
}

int resample_process_rvv(SpeexResamplerState *st, const float *in,
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

int resample_process_int_rvv(SpeexResamplerState *st, const spx_int16_t *in,
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

/* Interleaved multi-channel paths (RVV kernels). See wrap_resample_c.c. */
static void resample_reset_all_rvv(SpeexResamplerState *st)
{
    spx_uint32_t c;
    memset(st->mem, 0, (size_t) st->mem_alloc_size * st->nb_channels * sizeof(spx_word16_t));
    for (c = 0; c < st->nb_channels; c++) {
        st->last_sample[c]   = 0;
        st->samp_frac_num[c] = 0;
        st->magic_samples[c] = 0;
    }
    st->started = 0;
}

int resample_process_il_rvv(SpeexResamplerState *st, const float *in,
        spx_uint32_t in_len, float *out, spx_uint32_t out_len)
{
    spx_uint32_t il = in_len, ol = out_len;
    resample_reset_all_rvv(st);
    speex_resampler_process_interleaved_float(st, in, &il, out, &ol);
    return (int) ol;
}

int resample_process_int_il_rvv(SpeexResamplerState *st, const spx_int16_t *in,
        spx_uint32_t in_len, spx_int16_t *out, spx_uint32_t out_len)
{
    spx_uint32_t il = in_len, ol = out_len;
    resample_reset_all_rvv(st);
    speex_resampler_process_interleaved_int(st, in, &il, out, &ol);
    return (int) ol;
}
#endif /* !DISABLE_FLOAT_API */

#endif /* USE_RVV */
