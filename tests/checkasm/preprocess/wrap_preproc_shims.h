/* External-linkage shims around preprocess.c's static inline kernels, one per
 * kernel under test. Included by wrap_preproc_c.c and wrap_preproc_rvv.c after the
 * module, with CKA_VARIANT naming the TU (foo -> foo_c / foo_rvv); wrap.h
 * carries the matching prototypes. beta[0]/beta[1] = beta/beta_1, by
 * pointer for the same trampoline reason as the mdf scalars. */

void CKA_SHIM(preproc_window)(float *frame, const float *window, int len)
{
    preproc_window(frame, window, len);
}

void CKA_SHIM(preproc_power_spectrum)(const float *ft, float *ps, int N)
{
    preproc_power_spectrum(ft, ps, N);
}

void CKA_SHIM(preproc_smooth_spectrum)(float *S, const float *ps, int N)
{
    preproc_smooth_spectrum(S, ps, N);
}

void CKA_SHIM(preproc_min_track_swap)(float *Smin, float *Stmp, const float *S, int N)
{
    preproc_min_track_swap(Smin, Stmp, S, N);
}

void CKA_SHIM(preproc_min_track)(float *Smin, float *Stmp, const float *S, int N)
{
    preproc_min_track(Smin, Stmp, S, N);
}

void CKA_SHIM(preproc_update_prob)(const float *S, const float *Smin, int *update_prob, int N)
{
    preproc_update_prob(S, Smin, update_prob, N);
}

void CKA_SHIM(preproc_noise_update)(const int *update_prob, const float *ps, float *noise,
                                    const float *beta, int N)
{
    preproc_noise_update(update_prob, ps, noise, beta[0], beta[1], N);
}

void CKA_SHIM(preproc_snr_update)(const float *ps, const float *noise, const float *echo_noise,
                                  const float *reverb_estimate, const float *old_ps,
                                  float *post, float *prior, int len)
{
    preproc_snr_update(ps, noise, echo_noise, reverb_estimate, old_ps, post, prior, len);
}

void CKA_SHIM(preproc_zeta_smooth)(float *zeta, const float *prior, int N, int M)
{
    preproc_zeta_smooth(zeta, prior, N, M);
}

void CKA_SHIM(preproc_em_gain)(const float *prior, const float *post, const float *ps,
                               const float *gain_floor, float *gain, float *gain2,
                               float *old_ps, int N)
{
    preproc_em_gain(prior, post, ps, gain_floor, gain, gain2, old_ps, N);
}

void CKA_SHIM(preproc_apply_gain)(const float *gain2, float *ft, int N)
{
    preproc_apply_gain(gain2, ft, N);
}

void CKA_SHIM(preproc_overlap_output)(spx_int16_t *x, const float *outbuf,
                                      const float *frame, int len)
{
    preproc_overlap_output(x, outbuf, frame, len);
}
