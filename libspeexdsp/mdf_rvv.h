/* Copyright (C) 2003-2008 Jean-Marc Valin
 * Copyright (C) 2026 Tristan Matthews
 */
/**
   @file mdf_rvv.h
   @brief MDF echo canceller spectral kernels (RISC-V Vector extension,
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

/* Runtime-dispatched RVV kernels for mdf.c's O(M*N) spectral loops. The
 * vector code lives out-of-line in mdf_rvv_asm.S, so this header is plain
 * C and mdf.c stays base-ISA; each OVERRIDE_MDF_* wrapper calls a kernel
 * only when SPX_MDF_RVV_ON and falls back to the scalar loop otherwise.
 * The asm handles the interior complex bins of the packed spectrum; the
 * two real-only edge bins (0 and N-1) stay in C here. Fixed point is
 * bit-exact vs C; float pairs with a checkasm tolerance (FMAs, reordered
 * sums). Include from mdf.c AFTER the WEIGHT_SHIFT definition. checkasm
 * defines MDF_RVV_FORCE_ON to test the asm unconditionally. */

#ifndef MDF_RVV_H
#define MDF_RVV_H

#include "arch.h"

#if defined(FIXED_POINT) || defined(__riscv_float_abi_double)

#ifdef MDF_RVV_FORCE_ON
#  define SPX_MDF_RVV_ON 1
#else
#include "rvv_cpu.h"
extern int spx_mdf_rvv_enabled;   /* defined in mdf.c: -1 until its first init probes, then 0/1 */
#  define SPX_MDF_RVV_ON (spx_mdf_rvv_enabled > 0)
#  define MDF_RVV_RUNTIME 1          /* tells mdf.c to define+detect the flag */
#endif

#ifdef FIXED_POINT

void spx_mdf_rvv_smul_accum_i16(const spx_int16_t *X, const spx_int32_t *Y,
                                spx_int16_t *acc, int N, int M, int shift);
void spx_mdf_rvv_smul_accum16_i16(const spx_int16_t *X, const spx_int16_t *Y,
                                  spx_int16_t *acc, int N, int M, int shift);
void spx_mdf_rvv_power_spectrum_i16(const spx_int16_t *X, spx_int32_t *ps, int N);
void spx_mdf_rvv_res_window_i16(spx_int16_t *y, const spx_int16_t *window,
                                const spx_int16_t *last_y, int N);
void spx_mdf_rvv_res_scale_i32(spx_int32_t *residual_echo, int leak2, int len);
void spx_mdf_rvv_power_spectrum_accum_i16(const spx_int16_t *X, spx_int32_t *ps, int N);
spx_int32_t spx_mdf_rvv_inner_prod_i16(const spx_int16_t *x, const spx_int16_t *y,
                                       int len, int shift);
spx_int32_t spx_mdf_rvv_prop_sumsq_i16(const spx_int32_t *W, int len, int shift);
void spx_mdf_rvv_weight_update_i32(spx_int32_t *w, const spx_int32_t *phi, int N);

/* Kernels behind the build-independent wrappers below. */
#define SPX_MDF_RVV_WEIGHT_UPDATE         spx_mdf_rvv_weight_update_i32
#define SPX_MDF_RVV_RES_WINDOW            spx_mdf_rvv_res_window_i16
#define SPX_MDF_RVV_RES_SCALE             spx_mdf_rvv_res_scale_i32
#define SPX_MDF_RVV_POWER_SPECTRUM        spx_mdf_rvv_power_spectrum_i16
#define SPX_MDF_RVV_POWER_SPECTRUM_ACCUM  spx_mdf_rvv_power_spectrum_accum_i16

#define OVERRIDE_MDF_SPECTRAL_MUL_ACCUM
static inline void spectral_mul_accum(const spx_word16_t *X, const spx_word32_t *Y, spx_word16_t *acc, int N, int M)
{
   int i,j;
   spx_word32_t tmp1=0,tmp2=0;
   if (SPX_MDF_RVV_ON && N > 2 && !(N & 1) && M > 0)
   {
      spx_mdf_rvv_smul_accum_i16(X, Y, acc, N, M, WEIGHT_SHIFT);
      for (j=0;j<M;j++)
      {
         tmp1 = MAC16_16(tmp1, X[j*N],TOP16(Y[j*N]));
         tmp2 = MAC16_16(tmp2, X[(j+1)*N-1],TOP16(Y[(j+1)*N-1]));
      }
      acc[0] = PSHR32(tmp1,WEIGHT_SHIFT);
      acc[N-1] = PSHR32(tmp2,WEIGHT_SHIFT);
      return;
   }
   for (j=0;j<M;j++)
   {
      tmp1 = MAC16_16(tmp1, X[j*N],TOP16(Y[j*N]));
   }
   acc[0] = PSHR32(tmp1,WEIGHT_SHIFT);
   for (i=1;i<N-1;i+=2)
   {
      tmp1 = tmp2 = 0;
      for (j=0;j<M;j++)
      {
         tmp1 = SUB32(MAC16_16(tmp1, X[j*N+i],TOP16(Y[j*N+i])), MULT16_16(X[j*N+i+1],TOP16(Y[j*N+i+1])));
         tmp2 = MAC16_16(MAC16_16(tmp2, X[j*N+i+1],TOP16(Y[j*N+i])), X[j*N+i], TOP16(Y[j*N+i+1]));
      }
      acc[i] = PSHR32(tmp1,WEIGHT_SHIFT);
      acc[i+1] = PSHR32(tmp2,WEIGHT_SHIFT);
   }
   tmp1 = tmp2 = 0;
   for (j=0;j<M;j++)
   {
      tmp1 = MAC16_16(tmp1, X[(j+1)*N-1],TOP16(Y[(j+1)*N-1]));
   }
   acc[N-1] = PSHR32(tmp1,WEIGHT_SHIFT);
}

#define OVERRIDE_MDF_SPECTRAL_MUL_ACCUM16
static inline void spectral_mul_accum16(const spx_word16_t *X, const spx_word16_t *Y, spx_word16_t *acc, int N, int M)
{
   int i,j;
   spx_word32_t tmp1=0,tmp2=0;
   if (SPX_MDF_RVV_ON && N > 2 && !(N & 1) && M > 0)
   {
      spx_mdf_rvv_smul_accum16_i16(X, Y, acc, N, M, WEIGHT_SHIFT);
      for (j=0;j<M;j++)
      {
         tmp1 = MAC16_16(tmp1, X[j*N],Y[j*N]);
         tmp2 = MAC16_16(tmp2, X[(j+1)*N-1],Y[(j+1)*N-1]);
      }
      acc[0] = PSHR32(tmp1,WEIGHT_SHIFT);
      acc[N-1] = PSHR32(tmp2,WEIGHT_SHIFT);
      return;
   }
   for (j=0;j<M;j++)
   {
      tmp1 = MAC16_16(tmp1, X[j*N],Y[j*N]);
   }
   acc[0] = PSHR32(tmp1,WEIGHT_SHIFT);
   for (i=1;i<N-1;i+=2)
   {
      tmp1 = tmp2 = 0;
      for (j=0;j<M;j++)
      {
         tmp1 = SUB32(MAC16_16(tmp1, X[j*N+i],Y[j*N+i]), MULT16_16(X[j*N+i+1],Y[j*N+i+1]));
         tmp2 = MAC16_16(MAC16_16(tmp2, X[j*N+i+1],Y[j*N+i]), X[j*N+i], Y[j*N+i+1]);
      }
      acc[i] = PSHR32(tmp1,WEIGHT_SHIFT);
      acc[i+1] = PSHR32(tmp2,WEIGHT_SHIFT);
   }
   tmp1 = tmp2 = 0;
   for (j=0;j<M;j++)
   {
      tmp1 = MAC16_16(tmp1, X[(j+1)*N-1],Y[(j+1)*N-1]);
   }
   acc[N-1] = PSHR32(tmp1,WEIGHT_SHIFT);
}

#define OVERRIDE_MDF_INNER_PROD
static inline spx_word32_t mdf_inner_prod(const spx_word16_t *x, const spx_word16_t *y, int len)
{
   spx_word32_t sum=0;
   /* below ~32 elements the vsetvli/reduction overhead beats the gain */
   if (SPX_MDF_RVV_ON && len >= 32)
      return spx_mdf_rvv_inner_prod_i16(x, y, len, 6);
   len >>= 1;
   while(len--)
   {
      spx_word32_t part=0;
      part = MAC16_16(part,*x++,*y++);
      part = MAC16_16(part,*x++,*y++);
      /* HINT: If you had a 40-bit accumulator, you could shift only at the end */
      sum = ADD32(sum,SHR32(part,6));
   }
   return sum;
}

#define OVERRIDE_MDF_ADJUST_PROP
static inline void mdf_adjust_prop(const spx_word32_t *W, int N, int M, int P, spx_word16_t *prop)
{
   int i, j, p;
   spx_word16_t max_sum = 1;
   spx_word32_t prop_sum = 1;
   for (i=0;i<M;i++)
   {
      spx_word32_t tmp = 1;
      if (SPX_MDF_RVV_ON)
      {
         for (p=0;p<P;p++)
            tmp += spx_mdf_rvv_prop_sumsq_i16(&W[p*N*M + i*N], N, 18);
      }
      else
      {
         for (p=0;p<P;p++)
            for (j=0;j<N;j++)
               tmp += MULT16_16(EXTRACT16(SHR32(W[p*N*M + i*N+j],18)), EXTRACT16(SHR32(W[p*N*M + i*N+j],18)));
      }
      /* Just a security in case an overflow were to occur */
      tmp = MIN32(ABS32(tmp), 536870912);
      prop[i] = spx_sqrt(tmp);
      if (prop[i] > max_sum)
         max_sum = prop[i];
   }
   for (i=0;i<M;i++)
   {
      prop[i] += MULT16_16_Q15(QCONST16(.1f,15),max_sum);
      prop_sum += EXTEND32(prop[i]);
   }
   for (i=0;i<M;i++)
   {
      prop[i] = DIV32(MULT16_16(QCONST16(.99f,15),prop[i]),prop_sum);
   }
}

#else /* FLOATING_POINT && __riscv_float_abi_double */

void  spx_mdf_rvv_word2int_f32(spx_int16_t *out, const float *in, int len,
                               int stride);
void  spx_mdf_rvv_smul_accum_f32(const float *X, const float *Y, float *acc,
                                 int N, int M);
void  spx_mdf_rvv_wsmul_conj_f32(const float *w, const float *X, const float *Y,
                                 float *prod, int N, float p);
void  spx_mdf_rvv_power_spectrum_f32(const float *X, float *ps, int N);
void  spx_mdf_rvv_res_window_f32(float *y, const float *window,
                                 const float *last_y, int N);
void  spx_mdf_rvv_res_scale_f32(float *residual_echo, float leak2, int len);
void  spx_mdf_rvv_power_spectrum_accum_f32(const float *X, float *ps, int N);
float spx_mdf_rvv_inner_prod_f32(const float *x, const float *y, int len);
void  spx_mdf_rvv_weight_update_f32(float *w, const float *phi, int N);

/* Kernels behind the build-independent wrappers below. */
#define SPX_MDF_RVV_WEIGHT_UPDATE         spx_mdf_rvv_weight_update_f32
#define SPX_MDF_RVV_RES_WINDOW            spx_mdf_rvv_res_window_f32
#define SPX_MDF_RVV_RES_SCALE             spx_mdf_rvv_res_scale_f32
#define SPX_MDF_RVV_POWER_SPECTRUM        spx_mdf_rvv_power_spectrum_f32
#define SPX_MDF_RVV_POWER_SPECTRUM_ACCUM  spx_mdf_rvv_power_spectrum_accum_f32

/* Float only (no fixed-point twin): the fixed WORD2INT is a cheap branchy
 * clamp, but the float one calls floor() per sample. The de-emphasis
 * recurrence itself is serial (each tmp_out feeds the next through mem),
 * so it stays scalar; the recurrence runs in chunks through a small
 * buffer and the RVV kernel batch-converts each chunk. The chunked C
 * statements match the scalar fallback's exactly, and the kernel is
 * bit-exact vs WORD2INT (proven exhaustively over all floats), so RVV
 * output matches C bit for bit. */
#define OVERRIDE_MDF_DEEMPH_OUTPUT
static inline spx_word16_t mdf_deemph_output(spx_int16_t *out, const spx_word16_t *input, const spx_word16_t *e, spx_word16_t preemph, spx_word16_t mem, int len, int stride)
{
   int i;
   if (SPX_MDF_RVV_ON && len >= 8)
   {
      float tmp[256];
      while (len > 0)
      {
         int n = len < 256 ? len : 256;
         for (i=0;i<n;i++)
         {
            spx_word32_t tmp_out = SUB32(EXTEND32(input[i]), EXTEND32(e[i]));
            tmp_out = ADD32(tmp_out, EXTEND32(MULT16_16_P15(preemph, mem)));
            tmp[i] = tmp_out;
            mem = tmp_out;
         }
         spx_mdf_rvv_word2int_f32(out, tmp, n, stride);
         input += n;
         e += n;
         out += n*stride;
         len -= n;
      }
      return mem;
   }
   for (i=0;i<len;i++)
   {
      spx_word32_t tmp_out = SUB32(EXTEND32(input[i]), EXTEND32(e[i]));
      tmp_out = ADD32(tmp_out, EXTEND32(MULT16_16_P15(preemph, mem)));
      out[i*stride] = WORD2INT(tmp_out);
      mem = tmp_out;
   }
   return mem;
}

#define OVERRIDE_MDF_SPECTRAL_MUL_ACCUM
static inline void spectral_mul_accum(const spx_word16_t *X, const spx_word32_t *Y, spx_word16_t *acc, int N, int M)
{
   int i,j;
   if (SPX_MDF_RVV_ON && N > 2 && !(N & 1) && M > 0)
   {
      float a0 = 0, aN = 0;
      spx_mdf_rvv_smul_accum_f32(X, Y, acc, N, M);
      for (j=0;j<M;j++)
      {
         a0 += X[j*N]*Y[j*N];
         aN += X[(j+1)*N-1]*Y[(j+1)*N-1];
      }
      acc[0] = a0;
      acc[N-1] = aN;
      return;
   }
   for (i=0;i<N;i++)
      acc[i] = 0;
   for (j=0;j<M;j++)
   {
      acc[0] += X[0]*Y[0];
      for (i=1;i<N-1;i+=2)
      {
         acc[i] += (X[i]*Y[i] - X[i+1]*Y[i+1]);
         acc[i+1] += (X[i+1]*Y[i] + X[i]*Y[i+1]);
      }
      acc[i] += X[i]*Y[i];
      X += N;
      Y += N;
   }
}

#define OVERRIDE_MDF_WEIGHTED_SPECTRAL_MUL_CONJ
static inline void weighted_spectral_mul_conj(const spx_float_t *w, const spx_float_t p, const spx_word16_t *X, const spx_word16_t *Y, spx_word32_t *prod, int N)
{
   int i, j;
   spx_float_t W;
   if (SPX_MDF_RVV_ON && N > 2 && !(N & 1))
   {
      spx_mdf_rvv_wsmul_conj_f32(w, X, Y, prod, N, p);
      prod[0] = p*w[0]*(X[0]*Y[0]);
      prod[N-1] = p*w[N>>1]*(X[N-1]*Y[N-1]);
      return;
   }
   W = FLOAT_AMULT(p, w[0]);
   prod[0] = FLOAT_MUL32(W,MULT16_16(X[0],Y[0]));
   for (i=1,j=1;i<N-1;i+=2,j++)
   {
      W = FLOAT_AMULT(p, w[j]);
      prod[i] = FLOAT_MUL32(W,MAC16_16(MULT16_16(X[i],Y[i]), X[i+1],Y[i+1]));
      prod[i+1] = FLOAT_MUL32(W,MAC16_16(MULT16_16(-X[i+1],Y[i]), X[i],Y[i+1]));
   }
   W = FLOAT_AMULT(p, w[j]);
   prod[i] = FLOAT_MUL32(W,MULT16_16(X[i],Y[i]));
}

#define OVERRIDE_MDF_INNER_PROD
static inline spx_word32_t mdf_inner_prod(const spx_word16_t *x, const spx_word16_t *y, int len)
{
   spx_word32_t sum=0;
   /* below ~32 elements the vsetvli/reduction overhead beats the gain */
   if (SPX_MDF_RVV_ON && len >= 32)
      return spx_mdf_rvv_inner_prod_f32(x, y, len & ~1); /* C drops an odd tail */
   len >>= 1;
   while(len--)
   {
      spx_word32_t part=0;
      part = MAC16_16(part,*x++,*y++);
      part = MAC16_16(part,*x++,*y++);
      /* HINT: If you had a 40-bit accumulator, you could shift only at the end */
      sum = ADD32(sum,SHR32(part,6));
   }
   return sum;
}

#define OVERRIDE_MDF_ADJUST_PROP
static inline void mdf_adjust_prop(const spx_word32_t *W, int N, int M, int P, spx_word16_t *prop)
{
   int i, j, p;
   spx_word16_t max_sum = 1;
   spx_word32_t prop_sum = 1;
   for (i=0;i<M;i++)
   {
      spx_word32_t tmp = 1;
      if (SPX_MDF_RVV_ON)
      {
         for (p=0;p<P;p++)
            tmp += spx_mdf_rvv_inner_prod_f32(&W[p*N*M + i*N], &W[p*N*M + i*N], N);
      }
      else
      {
         for (p=0;p<P;p++)
            for (j=0;j<N;j++)
               tmp += MULT16_16(EXTRACT16(SHR32(W[p*N*M + i*N+j],18)), EXTRACT16(SHR32(W[p*N*M + i*N+j],18)));
      }
      prop[i] = spx_sqrt(tmp);
      if (prop[i] > max_sum)
         max_sum = prop[i];
   }
   for (i=0;i<M;i++)
   {
      prop[i] += MULT16_16_Q15(QCONST16(.1f,15),max_sum);
      prop_sum += EXTEND32(prop[i]);
   }
   for (i=0;i<M;i++)
   {
      prop[i] = DIV32(MULT16_16(QCONST16(.99f,15),prop[i]),prop_sum);
   }
}

#endif /* FIXED_POINT / float */

/* Wrappers whose C side is the same in both builds: the kernels differ only
 * in element type, selected through the SPX_MDF_RVV_* aliases above. */

#define OVERRIDE_MDF_WEIGHT_UPDATE
static inline void mdf_weight_update(spx_word32_t *w, const spx_word32_t *phi, int N)
{
   int i;
   if (SPX_MDF_RVV_ON && N >= 4)
   {
      SPX_MDF_RVV_WEIGHT_UPDATE(w, phi, N);
      return;
   }
   for (i=0;i<N;i++)
      w[i] += phi[i];
}

#define OVERRIDE_MDF_RESIDUAL_WINDOW
static inline void mdf_residual_window(spx_word16_t *y, const spx_word16_t *window, const spx_word16_t *last_y, int N)
{
   int i;
   if (SPX_MDF_RVV_ON && N >= 4)
   {
      SPX_MDF_RVV_RES_WINDOW(y, window, last_y, N);
      return;
   }
   for (i=0;i<N;i++)
      y[i] = MULT16_16_Q15(window[i], last_y[i]);
}

#define OVERRIDE_MDF_RESIDUAL_SCALE
static inline void mdf_residual_scale(spx_word32_t *residual_echo, spx_word16_t leak2, int len)
{
   int i;
   if (SPX_MDF_RVV_ON && len >= 4)
   {
      SPX_MDF_RVV_RES_SCALE(residual_echo, leak2, len);
      return;
   }
   for (i=0;i<len;i++)
      residual_echo[i] = (spx_int32_t)MULT16_32_Q15(leak2,residual_echo[i]);
}

#define OVERRIDE_MDF_POWER_SPECTRUM
static inline void power_spectrum(const spx_word16_t *X, spx_word32_t *ps, int N)
{
   int i, j;
   ps[0]=MULT16_16(X[0],X[0]);
   if (SPX_MDF_RVV_ON && N > 2 && !(N & 1))
   {
      SPX_MDF_RVV_POWER_SPECTRUM(X, ps, N);
      ps[N>>1]=MULT16_16(X[N-1],X[N-1]);
      return;
   }
   for (i=1,j=1;i<N-1;i+=2,j++)
   {
      ps[j] =  MULT16_16(X[i],X[i]) + MULT16_16(X[i+1],X[i+1]);
   }
   ps[j]=MULT16_16(X[i],X[i]);
}

#define OVERRIDE_MDF_POWER_SPECTRUM_ACCUM
static inline void power_spectrum_accum(const spx_word16_t *X, spx_word32_t *ps, int N)
{
   int i, j;
   ps[0]+=MULT16_16(X[0],X[0]);
   if (SPX_MDF_RVV_ON && N > 2 && !(N & 1))
   {
      SPX_MDF_RVV_POWER_SPECTRUM_ACCUM(X, ps, N);
      ps[N>>1]+=MULT16_16(X[N-1],X[N-1]);
      return;
   }
   for (i=1,j=1;i<N-1;i+=2,j++)
   {
      ps[j] +=  MULT16_16(X[i],X[i]) + MULT16_16(X[i+1],X[i+1]);
   }
   ps[j]+=MULT16_16(X[i],X[i]);
}

#endif /* FIXED_POINT || __riscv_float_abi_double */

#endif /* MDF_RVV_H */
