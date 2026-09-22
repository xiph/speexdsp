/* External-linkage shims shared by the wrap_resample_<variant>.c TUs,
 * with CKA_VARIANT naming the TU (foo -> foo_c / foo_sse / ...); wrap.h
 * carries the matching prototypes.
 *
 * CKA_RESAMPLER_BASIC_SHIM(kind) defines the shim for one of the four
 * static resampler_basic_<kind> kernels: same signature as
 * resampler_basic_func, resetting the per-channel cursor first so repeated
 * invocations -- in particular checkasm_bench_new's call loop -- redo
 * identical full work rather than starting past the end of the input.
 * Each TU instantiates it for the kernels its build overrides.
 *
 * The full-pipeline wrappers below drive the public process API; every
 * variant TU but CKA_KERNEL_SHIMS_ONLY ones defines them, since update_filter in that TU points
 * resampler_ptr at that TU's kernels (see wrap.h for the conventions:
 * lengths by value, filter memory zeroed and cursors reset per call). */

#define CKA_RESAMPLER_BASIC_SHIM(kind)                                        \
int CKA_SHIM(resampler_basic_##kind)(SpeexResamplerState *st,                 \
        spx_uint32_t channel_index, const spx_word16_t *in,                   \
        spx_uint32_t *in_len, spx_word16_t *out, spx_uint32_t *out_len)       \
{                                                                             \
    st->last_sample[channel_index]   = 0;                                     \
    st->samp_frac_num[channel_index] = 0;                                     \
    return resampler_basic_##kind(st, channel_index, in, in_len, out, out_len); \
}

#if !defined(CKA_KERNEL_SHIMS_ONLY)
SpeexResamplerState *CKA_SHIM(resample_make_state)(unsigned in_rate, unsigned out_rate, int quality)
{
    int err = RESAMPLER_ERR_SUCCESS;
    return speex_resampler_init(1, in_rate, out_rate, quality, &err);
}

SpeexResamplerState *CKA_CAT(CKA_SHIM(resample_make_state), _ch)(unsigned in_rate, unsigned out_rate,
        int quality, unsigned channels)
{
    int err = RESAMPLER_ERR_SUCCESS;
    return speex_resampler_init(channels, in_rate, out_rate, quality, &err);
}
#endif /* !CKA_KERNEL_SHIMS_ONLY */

#if !defined(DISABLE_FLOAT_API) && !defined(CKA_KERNEL_SHIMS_ONLY)
int CKA_SHIM(resample_process)(SpeexResamplerState *st, const float *in,
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

int CKA_SHIM(resample_process_int)(SpeexResamplerState *st, const spx_int16_t *in,
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

/* Interleaved multi-channel paths. in_len/out_len are per channel; the buffers
 * hold nb_channels * len. Reset every channel so each call redoes equal work. */
static void CKA_SHIM(resample_reset_all)(SpeexResamplerState *st)
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

int CKA_SHIM(resample_process_il)(SpeexResamplerState *st, const float *in,
        spx_uint32_t in_len, float *out, spx_uint32_t out_len)
{
    spx_uint32_t il = in_len, ol = out_len;
    CKA_SHIM(resample_reset_all)(st);
    speex_resampler_process_interleaved_float(st, in, &il, out, &ol);
    return (int) ol;
}

int CKA_SHIM(resample_process_int_il)(SpeexResamplerState *st, const spx_int16_t *in,
        spx_uint32_t in_len, spx_int16_t *out, spx_uint32_t out_len)
{
    spx_uint32_t il = in_len, ol = out_len;
    CKA_SHIM(resample_reset_all)(st);
    speex_resampler_process_interleaved_int(st, in, &il, out, &ol);
    return (int) ol;
}
#endif /* !DISABLE_FLOAT_API && !CKA_KERNEL_SHIMS_ONLY */
