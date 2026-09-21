/* Copyright (C) 2003 Epic Games (written by Jean-Marc Valin)
 * Copyright (C) 2004-2006 Epic Games
 * Copyright (C) 2026 Tristan Matthews
 */
/**
   @file preprocess_rvv.h
   @brief Preprocessor (denoiser) per-bin kernels (RISC-V Vector extension,
          runtime-dispatched)
*/
/*
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

   - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   - Neither the name of the Xiph.org Foundation nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE FOUNDATION OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/* Runtime-dispatched RVV kernels for preprocess.c's per-bin loops. The
 * vector code lives out-of-line in preprocess_rvv_asm.S, so this header is
 * plain C and preprocess.c stays base-ISA; each OVERRIDE_PREPROC_* wrapper
 * calls a kernel only when SPX_PREPROC_RVV_ON and falls back to the scalar
 * loop otherwise.
 *
 * Float only: the float loops map to plain vector arithmetic (vfdiv,
 * vfsqrt, and an indexed-load gather for the hypergeom_gain table), while
 * their fixed-point twins are built from branchy range-reduced divides
 * (DIV32_16_Q8/Q15) and the spx_exp/spx_sqrt polynomials, which have no
 * bit-exact vector mapping worth the maintenance. All float kernels pair
 * with a small checkasm tolerance (FMAs, 1/x refactoring, reordered
 * rounding).
 *
 * Include from preprocess.c AFTER the SNR_SHIFT/EXPIN_SHIFT macros and the
 * hypergeom_gain definition (the scalar fallbacks reference them).
 * checkasm defines PREPROC_RVV_FORCE_ON to test the asm unconditionally. */

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

#define OVERRIDE_PREPROC_WINDOW
static inline void preproc_window(spx_word16_t *frame, const spx_word16_t *window, int len)
{
   int i;
   if (SPX_PREPROC_RVV_ON && len >= 4)
   {
      spx_preproc_rvv_window_f32(frame, window, len);
      return;
   }
   for (i=0;i<len;i++)
      frame[i] = MULT16_16_Q15(frame[i], window[i]);
}

#define OVERRIDE_PREPROC_POWER_SPECTRUM
static inline void preproc_power_spectrum(const spx_word16_t *ft, spx_word32_t *ps, int N)
{
   int i;
   ps[0]=MULT16_16(ft[0],ft[0]);
   if (SPX_PREPROC_RVV_ON && N >= 4)
   {
      /* interior complex bins: ps[1..N-1] */
      spx_preproc_rvv_power_spectrum_f32(ft, ps, N);
      return;
   }
   for (i=1;i<N;i++)
      ps[i]=MULT16_16(ft[2*i-1],ft[2*i-1]) + MULT16_16(ft[2*i],ft[2*i]);
}

#define OVERRIDE_PREPROC_SMOOTH_SPECTRUM
static inline void preproc_smooth_spectrum(spx_word32_t *S, const spx_word32_t *ps, int N)
{
   int i;
   if (SPX_PREPROC_RVV_ON && N >= 4)
   {
      /* the kernel updates the interior bins S[1..N-2] */
      spx_preproc_rvv_smooth_spectrum_f32(S, ps, N);
   } else {
      for (i=1;i<N-1;i++)
         S[i] =  MULT16_32_Q15(QCONST16(.8f,15),S[i]) + MULT16_32_Q15(QCONST16(.05f,15),ps[i-1])
                         + MULT16_32_Q15(QCONST16(.1f,15),ps[i]) + MULT16_32_Q15(QCONST16(.05f,15),ps[i+1]);
   }
   S[0] =  MULT16_32_Q15(QCONST16(.8f,15),S[0]) + MULT16_32_Q15(QCONST16(.2f,15),ps[0]);
   S[N-1] =  MULT16_32_Q15(QCONST16(.8f,15),S[N-1]) + MULT16_32_Q15(QCONST16(.2f,15),ps[N-1]);
}

#define OVERRIDE_PREPROC_MIN_TRACK_SWAP
static inline void preproc_min_track_swap(spx_word32_t *Smin, spx_word32_t *Stmp, const spx_word32_t *S, int N)
{
   int i;
   if (SPX_PREPROC_RVV_ON && N >= 4)
   {
      spx_preproc_rvv_min_track_swap_f32(Smin, Stmp, S, N);
      return;
   }
   for (i=0;i<N;i++)
   {
      Smin[i] = MIN32(Stmp[i], S[i]);
      Stmp[i] = S[i];
   }
}

#define OVERRIDE_PREPROC_MIN_TRACK
static inline void preproc_min_track(spx_word32_t *Smin, spx_word32_t *Stmp, const spx_word32_t *S, int N)
{
   int i;
   if (SPX_PREPROC_RVV_ON && N >= 4)
   {
      spx_preproc_rvv_min_track_f32(Smin, Stmp, S, N);
      return;
   }
   for (i=0;i<N;i++)
   {
      Smin[i] = MIN32(Smin[i], S[i]);
      Stmp[i] = MIN32(Stmp[i], S[i]);
   }
}

#define OVERRIDE_PREPROC_UPDATE_PROB
static inline void preproc_update_prob(const spx_word32_t *S, const spx_word32_t *Smin, int *update_prob, int N)
{
   int i;
   if (SPX_PREPROC_RVV_ON && N >= 4)
   {
      spx_preproc_rvv_update_prob_f32(S, Smin, update_prob, N);
      return;
   }
   for (i=0;i<N;i++)
   {
      if (MULT16_32_Q15(QCONST16(.4f,15),S[i]) > Smin[i])
         update_prob[i] = 1;
      else
         update_prob[i] = 0;
   }
}

#define OVERRIDE_PREPROC_NOISE_UPDATE
static inline void preproc_noise_update(const int *update_prob, const spx_word32_t *ps, spx_word32_t *noise, spx_word16_t beta, spx_word16_t beta_1, int N)
{
   int i;
   if (SPX_PREPROC_RVV_ON && N >= 4)
   {
      spx_preproc_rvv_noise_update_f32(update_prob, ps, noise, beta, beta_1, N);
      return;
   }
   for (i=0;i<N;i++)
   {
      if (!update_prob[i] || ps[i] < PSHR32(noise[i], NOISE_SHIFT))
         noise[i] = MAX32(EXTEND32(0),MULT16_32_Q15(beta_1,noise[i]) + MULT16_32_Q15(beta,SHL32(ps[i],NOISE_SHIFT)));
   }
}

#define OVERRIDE_PREPROC_SNR_UPDATE
static inline void preproc_snr_update(const spx_word32_t *ps, const spx_word32_t *noise, const spx_word32_t *echo_noise, const spx_word32_t *reverb_estimate, const spx_word32_t *old_ps, spx_word16_t *post, spx_word16_t *prior, int len)
{
   int i;
   if (SPX_PREPROC_RVV_ON && len >= 4)
   {
      spx_preproc_rvv_snr_f32(ps, noise, echo_noise, reverb_estimate, old_ps, post, prior, len);
      return;
   }
   for (i=0;i<len;i++)
   {
      spx_word16_t gamma;
      spx_word32_t tot_noise = ADD32(ADD32(ADD32(EXTEND32(1), PSHR32(noise[i],NOISE_SHIFT)) , echo_noise[i]) , reverb_estimate[i]);
      post[i] = SUB16(DIV32_16_Q8(ps[i],tot_noise), QCONST16(1.f,SNR_SHIFT));
      post[i]=MIN16(post[i], QCONST16(100.f,SNR_SHIFT));
      gamma = QCONST16(.1f,15)+MULT16_16_Q15(QCONST16(.89f,15),SQR16_Q15(DIV32_16_Q15(old_ps[i],ADD32(old_ps[i],tot_noise))));
      prior[i] = EXTRACT16(PSHR32(ADD32(MULT16_16(gamma,MAX16(0,post[i])), MULT16_16(Q15_ONE-gamma,DIV32_16_Q8(old_ps[i],tot_noise))), 15));
      prior[i]=MIN16(prior[i], QCONST16(100.f,SNR_SHIFT));
   }
}

#define OVERRIDE_PREPROC_ZETA_SMOOTH
static inline void preproc_zeta_smooth(spx_word16_t *zeta, const spx_word16_t *prior, int N, int M)
{
   int i;
   zeta[0] = PSHR32(ADD32(MULT16_16(QCONST16(.7f,15),zeta[0]), MULT16_16(QCONST16(.3f,15),prior[0])),15);
   if (SPX_PREPROC_RVV_ON && N >= 4)
   {
      /* the kernel updates the interior bins zeta[1..N-2] */
      spx_preproc_rvv_zeta_f32(zeta, prior, N);
   } else {
      for (i=1;i<N-1;i++)
         zeta[i] = PSHR32(ADD32(ADD32(ADD32(MULT16_16(QCONST16(.7f,15),zeta[i]), MULT16_16(QCONST16(.15f,15),prior[i])),
                              MULT16_16(QCONST16(.075f,15),prior[i-1])), MULT16_16(QCONST16(.075f,15),prior[i+1])),15);
   }
   for (i=N-1;i<N+M;i++)
      zeta[i] = PSHR32(ADD32(MULT16_16(QCONST16(.7f,15),zeta[i]), MULT16_16(QCONST16(.3f,15),prior[i])),15);
}

#define OVERRIDE_PREPROC_EM_GAIN
static inline void preproc_em_gain(const spx_word16_t *prior, const spx_word16_t *post, const spx_word32_t *ps, const spx_word16_t *gain_floor, spx_word16_t *gain, spx_word16_t *gain2, spx_word32_t *old_ps, int N)
{
   int i;
   if (SPX_PREPROC_RVV_ON && N >= 4)
   {
      spx_preproc_rvv_em_gain_f32(prior, post, ps, gain_floor, gain, gain2, old_ps, N);
      return;
   }
   for (i=0;i<N;i++)
   {
      spx_word32_t MM;
      spx_word32_t theta;
      spx_word16_t prior_ratio;
      spx_word16_t tmp;
      spx_word16_t p;
      spx_word16_t g;

      prior_ratio = PDIV32_16(SHL32(EXTEND32(prior[i]), 15), ADD16(prior[i], SHL32(1,SNR_SHIFT)));
      theta = MULT16_32_P15(prior_ratio, QCONST32(1.f,EXPIN_SHIFT)+SHL32(EXTEND32(post[i]),EXPIN_SHIFT-SNR_SHIFT));

      MM = hypergeom_gain(theta);
      g = EXTRACT16(MIN32(Q15_ONE, MULT16_32_Q15(prior_ratio, MM)));
      p = gain2[i];

      if (MULT16_16_Q15(QCONST16(.333f,15),g) > gain[i])
         g = MULT16_16(3,gain[i]);
      gain[i] = g;

      old_ps[i] = MULT16_32_P15(QCONST16(.2f,15),old_ps[i]) + MULT16_32_P15(MULT16_16_P15(QCONST16(.8f,15),SQR16_Q15(gain[i])),ps[i]);

      if (gain[i] < gain_floor[i])
         gain[i] = gain_floor[i];

      tmp = MULT16_16_P15(p,spx_sqrt(SHL32(EXTEND32(gain[i]),15))) + MULT16_16_P15(SUB16(Q15_ONE,p),spx_sqrt(SHL32(EXTEND32(gain_floor[i]),15)));
      gain2[i]=SQR16_Q15(tmp);
   }
}

#define OVERRIDE_PREPROC_APPLY_GAIN
static inline void preproc_apply_gain(const spx_word16_t *gain2, spx_word16_t *ft, int N)
{
   int i;
   if (SPX_PREPROC_RVV_ON && N >= 4)
   {
      /* the kernel scales the interior complex bins ft[1..2*N-2] */
      spx_preproc_rvv_apply_gain_f32(gain2, ft, N);
   } else {
      for (i=1;i<N;i++)
      {
         ft[2*i-1] = MULT16_16_P15(gain2[i],ft[2*i-1]);
         ft[2*i] = MULT16_16_P15(gain2[i],ft[2*i]);
      }
   }
   ft[0] = MULT16_16_P15(gain2[0],ft[0]);
   ft[2*N-1] = MULT16_16_P15(gain2[N-1],ft[2*N-1]);
}

/* Bit-exact vs the scalar loop: the kernel adds a+b under the ambient
 * rounding mode (matching C's +), then reproduces WORD2INT's
 * floor(.5+x) with a clamp and an frm=RDN add+convert (equivalence
 * proven exhaustively over all floats). */
#define OVERRIDE_PREPROC_OVERLAP_OUTPUT
static inline void preproc_overlap_output(spx_int16_t *x, const spx_word16_t *outbuf, const spx_word16_t *frame, int len)
{
   int i;
   if (SPX_PREPROC_RVV_ON && len >= 4)
   {
      spx_preproc_rvv_word2int_sum_f32(x, outbuf, frame, len);
      return;
   }
   for (i=0;i<len;i++)
      x[i] = WORD2INT(ADD32(EXTEND32(outbuf[i]), EXTEND32(frame[i])));
}

#endif /* !FIXED_POINT && __riscv_float_abi_double */

#endif /* PREPROC_RVV_H */
