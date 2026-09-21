/* Copyright (C) 2003 Epic Games (written by Jean-Marc Valin)
 * Copyright (C) 2004-2006 Epic Games
 * Copyright (C) 2026 Tristan Matthews
 */
/**
   @file preprocess_rvv.h
   @brief Preprocessor (denoiser) per-bin kernels (RISC-V Vector extension,
          runtime-dispatched)
*/

/* Runtime-dispatched RVV kernels for preprocess.c's per-bin loops. The
 * vector code lives out-of-line in preprocess_rvv_asm.S, so this header is
 * plain C and preprocess.c stays base-ISA. preprocess.c includes it right
 * after its scalar kernels: each <kernel>_dispatch wrapper below calls the
 * asm when SPX_PREPROC_RVV_ON (and the frame is at least a few bins) and
 * the scalar kernel otherwise, and `#define <kernel> <kernel>_dispatch`
 * then points every later call site at the wrapper.
 *
 * Float only: the float loops map to plain vector arithmetic (vfdiv,
 * vfsqrt, and an indexed-load gather for the hypergeom_gain table), while
 * their fixed-point twins are built from branchy range-reduced divides
 * (DIV32_16_Q8/Q15) and the spx_exp/spx_sqrt polynomials, which have no
 * bit-exact vector mapping worth the maintenance. All float kernels pair
 * with a small checkasm tolerance (FMAs, 1/x refactoring, reordered
 * rounding). checkasm defines PREPROC_RVV_FORCE_ON to test the asm
 * unconditionally. */

#ifndef PREPROC_RVV_H
#define PREPROC_RVV_H

#include "arch.h"

#if !defined(FIXED_POINT) && defined(__riscv_float_abi_double)

#ifdef PREPROC_RVV_FORCE_ON
#  define SPX_PREPROC_RVV_ON 1
#else
#include "rvv_cpu.h"
extern int spx_preproc_rvv_enabled;   /* defined in preprocess.c: -1 until its first init probes, then 0/1 */
#  define SPX_PREPROC_RVV_ON (spx_preproc_rvv_enabled > 0)
#  define PREPROC_RVV_RUNTIME 1      /* tells preprocess.c to define+detect the flag */
#endif

void spx_preproc_rvv_window_f32(float *frame, const float *window, int len);
void spx_preproc_rvv_power_spectrum_f32(const float *ft, float *ps, int N);
void spx_preproc_rvv_smooth_spectrum_f32(float *S, const float *ps, int N);
void spx_preproc_rvv_min_track_swap_f32(float *Smin, float *Stmp, const float *S, int N);
void spx_preproc_rvv_min_track_f32(float *Smin, float *Stmp, const float *S, int N);
void spx_preproc_rvv_update_prob_f32(const float *S, const float *Smin, int *update_prob, int N);
void spx_preproc_rvv_noise_update_f32(const int *update_prob, const float *ps, float *noise,
                                      float beta, float beta_1, int N);
void spx_preproc_rvv_snr_f32(const float *ps, const float *noise, const float *echo_noise,
                             const float *reverb_estimate, const float *old_ps,
                             float *post, float *prior, int len);
void spx_preproc_rvv_zeta_f32(float *zeta, const float *prior, int N);
void spx_preproc_rvv_em_gain_f32(const float *prior, const float *post, const float *ps,
                                 const float *gain_floor, float *gain, float *gain2,
                                 float *old_ps, int N);
void spx_preproc_rvv_apply_gain_f32(const float *gain2, float *ft, int N);
void spx_preproc_rvv_word2int_sum_f32(spx_int16_t *x, const float *a,
                                      const float *b, int len);

static inline void preproc_window_dispatch(spx_word16_t *frame, const spx_word16_t *window, int len)
{
   if (SPX_PREPROC_RVV_ON && len >= 4)
      spx_preproc_rvv_window_f32(frame, window, len);
   else
      preproc_window(frame, window, len);
}
#define preproc_window preproc_window_dispatch

static inline void preproc_power_spectrum_dispatch(const spx_word16_t *ft, spx_word32_t *ps, int N)
{
   if (SPX_PREPROC_RVV_ON && N >= 4)
   {
      ps[0]=MULT16_16(ft[0],ft[0]);
      /* interior complex bins: ps[1..N-1] */
      spx_preproc_rvv_power_spectrum_f32(ft, ps, N);
   }
   else
      preproc_power_spectrum(ft, ps, N);
}
#define preproc_power_spectrum preproc_power_spectrum_dispatch

static inline void preproc_smooth_spectrum_dispatch(spx_word32_t *S, const spx_word32_t *ps, int N)
{
   if (SPX_PREPROC_RVV_ON && N >= 4)
   {
      /* the kernel updates the interior bins S[1..N-2] */
      spx_preproc_rvv_smooth_spectrum_f32(S, ps, N);
      S[0] =  MULT16_32_Q15(QCONST16(.8f,15),S[0]) + MULT16_32_Q15(QCONST16(.2f,15),ps[0]);
      S[N-1] =  MULT16_32_Q15(QCONST16(.8f,15),S[N-1]) + MULT16_32_Q15(QCONST16(.2f,15),ps[N-1]);
   }
   else
      preproc_smooth_spectrum(S, ps, N);
}
#define preproc_smooth_spectrum preproc_smooth_spectrum_dispatch

static inline void preproc_min_track_swap_dispatch(spx_word32_t *Smin, spx_word32_t *Stmp, const spx_word32_t *S, int N)
{
   if (SPX_PREPROC_RVV_ON && N >= 4)
      spx_preproc_rvv_min_track_swap_f32(Smin, Stmp, S, N);
   else
      preproc_min_track_swap(Smin, Stmp, S, N);
}
#define preproc_min_track_swap preproc_min_track_swap_dispatch

static inline void preproc_min_track_dispatch(spx_word32_t *Smin, spx_word32_t *Stmp, const spx_word32_t *S, int N)
{
   if (SPX_PREPROC_RVV_ON && N >= 4)
      spx_preproc_rvv_min_track_f32(Smin, Stmp, S, N);
   else
      preproc_min_track(Smin, Stmp, S, N);
}
#define preproc_min_track preproc_min_track_dispatch

static inline void preproc_update_prob_dispatch(const spx_word32_t *S, const spx_word32_t *Smin, int *update_prob, int N)
{
   if (SPX_PREPROC_RVV_ON && N >= 4)
      spx_preproc_rvv_update_prob_f32(S, Smin, update_prob, N);
   else
      preproc_update_prob(S, Smin, update_prob, N);
}
#define preproc_update_prob preproc_update_prob_dispatch

static inline void preproc_noise_update_dispatch(const int *update_prob, const spx_word32_t *ps, spx_word32_t *noise, spx_word16_t beta, spx_word16_t beta_1, int N)
{
   if (SPX_PREPROC_RVV_ON && N >= 4)
      spx_preproc_rvv_noise_update_f32(update_prob, ps, noise, beta, beta_1, N);
   else
      preproc_noise_update(update_prob, ps, noise, beta, beta_1, N);
}
#define preproc_noise_update preproc_noise_update_dispatch

static inline void preproc_snr_update_dispatch(const spx_word32_t *ps, const spx_word32_t *noise, const spx_word32_t *echo_noise, const spx_word32_t *reverb_estimate, const spx_word32_t *old_ps, spx_word16_t *post, spx_word16_t *prior, int len)
{
   if (SPX_PREPROC_RVV_ON && len >= 4)
      spx_preproc_rvv_snr_f32(ps, noise, echo_noise, reverb_estimate, old_ps, post, prior, len);
   else
      preproc_snr_update(ps, noise, echo_noise, reverb_estimate, old_ps, post, prior, len);
}
#define preproc_snr_update preproc_snr_update_dispatch

static inline void preproc_zeta_smooth_dispatch(spx_word16_t *zeta, const spx_word16_t *prior, int N, int M)
{
   if (SPX_PREPROC_RVV_ON && N >= 4)
   {
      int i;
      zeta[0] = PSHR32(ADD32(MULT16_16(QCONST16(.7f,15),zeta[0]), MULT16_16(QCONST16(.3f,15),prior[0])),15);
      /* the kernel updates the interior bins zeta[1..N-2] */
      spx_preproc_rvv_zeta_f32(zeta, prior, N);
      for (i=N-1;i<N+M;i++)
         zeta[i] = PSHR32(ADD32(MULT16_16(QCONST16(.7f,15),zeta[i]), MULT16_16(QCONST16(.3f,15),prior[i])),15);
   }
   else
      preproc_zeta_smooth(zeta, prior, N, M);
}
#define preproc_zeta_smooth preproc_zeta_smooth_dispatch

static inline void preproc_em_gain_dispatch(const spx_word16_t *prior, const spx_word16_t *post, const spx_word32_t *ps, const spx_word16_t *gain_floor, spx_word16_t *gain, spx_word16_t *gain2, spx_word32_t *old_ps, int N)
{
   if (SPX_PREPROC_RVV_ON && N >= 4)
      spx_preproc_rvv_em_gain_f32(prior, post, ps, gain_floor, gain, gain2, old_ps, N);
   else
      preproc_em_gain(prior, post, ps, gain_floor, gain, gain2, old_ps, N);
}
#define preproc_em_gain preproc_em_gain_dispatch

static inline void preproc_apply_gain_dispatch(const spx_word16_t *gain2, spx_word16_t *ft, int N)
{
   if (SPX_PREPROC_RVV_ON && N >= 4)
   {
      /* the kernel scales the interior complex bins ft[1..2*N-2] */
      spx_preproc_rvv_apply_gain_f32(gain2, ft, N);
      ft[0] = MULT16_16_P15(gain2[0],ft[0]);
      ft[2*N-1] = MULT16_16_P15(gain2[N-1],ft[2*N-1]);
   }
   else
      preproc_apply_gain(gain2, ft, N);
}
#define preproc_apply_gain preproc_apply_gain_dispatch

/* Bit-exact vs the scalar loop: the kernel adds a+b under the ambient
 * rounding mode (matching C's +), then reproduces WORD2INT's
 * floor(.5+x) with a clamp and an frm=RDN add+convert (equivalence
 * proven exhaustively over all floats). */
static inline void preproc_overlap_output_dispatch(spx_int16_t *x, const spx_word16_t *outbuf, const spx_word16_t *frame, int len)
{
   if (SPX_PREPROC_RVV_ON && len >= 4)
      spx_preproc_rvv_word2int_sum_f32(x, outbuf, frame, len);
   else
      preproc_overlap_output(x, outbuf, frame, len);
}
#define preproc_overlap_output preproc_overlap_output_dispatch

#endif /* !FIXED_POINT && __riscv_float_abi_double */

#endif /* PREPROC_RVV_H */
