#ifndef SPEEXDSP_TESTS_CHECKASM_PREPROCESS_WRAP_H
#define SPEEXDSP_TESTS_CHECKASM_PREPROCESS_WRAP_H

#include <string.h>
#include "config.h"
#include "arch.h"
#include "../compare.h"

/* Tests for preprocess.c's per-bin kernels (windowing, power spectrum,
 * noise tracking, SNR update, Ephraim-Malah gain, gain application), C vs
 * SIMD. Like the mdf tests, each variant #includes preprocess.c with its
 * public API renamed (wrap_preproc_rename.h): wrap_preproc_c.c forces the
 * scalar loops, wrap_preproc_rvv.c pins the RVV dispatch on. */

/* ------------- Per-ISA availability gates -------------
 * Mirrors preprocess_rvv.h's gate: float only, on an FP ABI. */
#if defined(USE_RVV) && !defined(FIXED_POINT) && defined(__riscv_float_abi_double)
#  define HAVE_RVV_PREPROC 1
#endif

/* ------------- Functions under test -------------
 * Thin shims around preprocess.c's static inline kernels; one per
 * variant TU. Only meaningful for float builds (the RVV kernels are
 * float-only), so the shims use plain float/int signatures. */
#ifndef FIXED_POINT

void preproc_window_c(float *frame, const float *window, int len);
void preproc_power_spectrum_c(const float *ft, float *ps, int N);
void preproc_smooth_spectrum_c(float *S, const float *ps, int N);
void preproc_min_track_swap_c(float *Smin, float *Stmp, const float *S, int N);
void preproc_min_track_c(float *Smin, float *Stmp, const float *S, int N);
void preproc_update_prob_c(const float *S, const float *Smin, int *update_prob, int N);
/* beta[0] = beta, beta[1] = beta_1; passed by pointer so the checked-call
 * trampoline never has to marshal float register args (same convention as
 * the mdf wsmul p argument). */
void preproc_noise_update_c(const int *update_prob, const float *ps, float *noise,
                            const float *beta, int N);
void preproc_snr_update_c(const float *ps, const float *noise, const float *echo_noise,
                          const float *reverb_estimate, const float *old_ps,
                          float *post, float *prior, int len);
void preproc_zeta_smooth_c(float *zeta, const float *prior, int N, int M);
void preproc_em_gain_c(const float *prior, const float *post, const float *ps,
                       const float *gain_floor, float *gain, float *gain2,
                       float *old_ps, int N);
void preproc_apply_gain_c(const float *gain2, float *ft, int N);
void preproc_overlap_output_c(spx_int16_t *x, const float *outbuf,
                              const float *frame, int len);

#ifdef HAVE_RVV_PREPROC
void preproc_window_rvv(float *frame, const float *window, int len);
void preproc_power_spectrum_rvv(const float *ft, float *ps, int N);
void preproc_smooth_spectrum_rvv(float *S, const float *ps, int N);
void preproc_min_track_swap_rvv(float *Smin, float *Stmp, const float *S, int N);
void preproc_min_track_rvv(float *Smin, float *Stmp, const float *S, int N);
void preproc_update_prob_rvv(const float *S, const float *Smin, int *update_prob, int N);
void preproc_noise_update_rvv(const int *update_prob, const float *ps, float *noise,
                              const float *beta, int N);
void preproc_snr_update_rvv(const float *ps, const float *noise, const float *echo_noise,
                            const float *reverb_estimate, const float *old_ps,
                            float *post, float *prior, int len);
void preproc_zeta_smooth_rvv(float *zeta, const float *prior, int N, int M);
void preproc_em_gain_rvv(const float *prior, const float *post, const float *ps,
                         const float *gain_floor, float *gain, float *gain2,
                         float *old_ps, int N);
void preproc_apply_gain_rvv(const float *gain2, float *ft, int N);
void preproc_overlap_output_rvv(spx_int16_t *x, const float *outbuf,
                                const float *frame, int len);
#endif

/* ------------- Test-input fill -------------
 * The kernels assume the value ranges the algorithm maintains (power
 * spectra are non-negative, SNRs are clamped at 100, gains sit in
 * [0, 1]); raw random bits would push them into states preprocess.c can
 * never reach (NaNs, negative powers), so fill accordingly. */
#define preproc_fill_signal checkasm_fill_symmetric_f32   /* [-1, 1) */

static inline void preproc_fill_power(float *buf, int n) /* [0, 1e6) */
{
    checkasm_randomize_rangef(buf, n, 1e6f);
}

static inline void preproc_fill_snr(float *buf, int n) /* [0, 100] incl. clamp value */
{
    checkasm_randomize_rangef(buf, n, 101.0f);
    for (int i = 0; i < n; i++)
        if (buf[i] > 100.0f)
            buf[i] = 100.0f;
}

static inline void preproc_fill_gain(float *buf, int n) /* [0, 1) */
{
    checkasm_randomize_rangef(buf, n, 1.0f);
}

static inline void preproc_fill_prob(int *buf, int n) /* 0/1 */
{
    for (int i = 0; i < n; i++)
        buf[i] = checkasm_rand() & 1;
}

/* ------------- Output comparison -------------
 * The RVV kernels use FMAs, a shared 1/x, and reordered sums, so compare
 * with checkasm_f32_within_tol relative to the buffer peak. em_gain also
 * crosses the .333*g > gain branch, whose two sides differ by 1e-3 at the
 * threshold (3*.333 = .999), hence its looser bound. The int16 output of
 * preproc_overlap_output and the update_prob flags are bit-exact
 * (checkasm_i16_bitexact / checkasm_i32_bitexact). */
#define PREPROC_ELTWISE_F32_REL_TOL 1e-6
#define PREPROC_SNR_F32_REL_TOL     1e-5
#define PREPROC_EM_F32_REL_TOL      2e-3

#endif /* !FIXED_POINT */

#endif /* SPEEXDSP_TESTS_CHECKASM_PREPROCESS_WRAP_H */
