/* Copyright (C) 2003-2008 Jean-Marc Valin
 * Copyright (C) 2026 Tristan Matthews
 */
/**
   @file mdf_rvv.h
   @brief MDF echo canceller spectral kernels (RISC-V Vector extension,
          runtime-dispatched)
*/

/* Runtime-dispatched RVV kernels for mdf.c's O(M*N) spectral loops. The
 * vector code lives out-of-line in mdf_rvv_asm.S, so this header is plain
 * C and mdf.c stays base-ISA. mdf.c includes it right after its scalar
 * kernels: each <kernel>_dispatch wrapper below calls the asm when
 * SPX_MDF_RVV_ON (and the shape is worth vectorizing) and the scalar
 * kernel otherwise, and `#define <kernel> <kernel>_dispatch` then points
 * every later call site at the wrapper. The asm handles the interior
 * complex bins of the packed spectrum; the two real-only edge bins (0 and
 * N-1) stay in C here. Fixed point is bit-exact vs C; float pairs with a
 * checkasm tolerance (FMAs, reordered sums). checkasm defines
 * MDF_RVV_FORCE_ON to test the asm unconditionally. */

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
#define SPX_MDF_RVV_INNER_PROD(x, y, len) spx_mdf_rvv_inner_prod_i16(x, y, len, 6)
#define SPX_MDF_RVV_PROP_SUMSQ(w, len)    spx_mdf_rvv_prop_sumsq_i16(w, len, 18)

static inline void spectral_mul_accum_dispatch(const spx_word16_t *X, const spx_word32_t *Y, spx_word16_t *acc, int N, int M)
{
   if (SPX_MDF_RVV_ON && N > 2 && !(N & 1) && M > 0)
   {
      int j;
      spx_word32_t tmp1=0,tmp2=0;
      spx_mdf_rvv_smul_accum_i16(X, Y, acc, N, M, WEIGHT_SHIFT);
      for (j=0;j<M;j++)
      {
         tmp1 = MAC16_16(tmp1, X[j*N],TOP16(Y[j*N]));
         tmp2 = MAC16_16(tmp2, X[(j+1)*N-1],TOP16(Y[(j+1)*N-1]));
      }
      acc[0] = PSHR32(tmp1,WEIGHT_SHIFT);
      acc[N-1] = PSHR32(tmp2,WEIGHT_SHIFT);
   }
   else
      spectral_mul_accum(X, Y, acc, N, M);
}
#define spectral_mul_accum spectral_mul_accum_dispatch

static inline void spectral_mul_accum16_dispatch(const spx_word16_t *X, const spx_word16_t *Y, spx_word16_t *acc, int N, int M)
{
   if (SPX_MDF_RVV_ON && N > 2 && !(N & 1) && M > 0)
   {
      int j;
      spx_word32_t tmp1=0,tmp2=0;
      spx_mdf_rvv_smul_accum16_i16(X, Y, acc, N, M, WEIGHT_SHIFT);
      for (j=0;j<M;j++)
      {
         tmp1 = MAC16_16(tmp1, X[j*N],Y[j*N]);
         tmp2 = MAC16_16(tmp2, X[(j+1)*N-1],Y[(j+1)*N-1]);
      }
      acc[0] = PSHR32(tmp1,WEIGHT_SHIFT);
      acc[N-1] = PSHR32(tmp2,WEIGHT_SHIFT);
   }
   else
      spectral_mul_accum16(X, Y, acc, N, M);
}
#define spectral_mul_accum16 spectral_mul_accum16_dispatch

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
/* C drops an odd tail element */
#define SPX_MDF_RVV_INNER_PROD(x, y, len) spx_mdf_rvv_inner_prod_f32(x, y, (len) & ~1)
#define SPX_MDF_RVV_PROP_SUMSQ(w, len)    spx_mdf_rvv_inner_prod_f32(w, w, len)

/* Float only (no fixed-point twin): the fixed WORD2INT is a cheap branchy
 * clamp, but the float one calls floor() per sample. The de-emphasis
 * recurrence itself is serial (each tmp_out feeds the next through mem),
 * so it stays scalar; the recurrence runs in chunks through a small
 * buffer and the RVV kernel batch-converts each chunk. The chunked C
 * statements match the scalar kernel's exactly, and the kernel is
 * bit-exact vs WORD2INT (proven exhaustively over all floats), so RVV
 * output matches C bit for bit. */
static inline spx_word16_t mdf_deemph_output_dispatch(spx_int16_t *out, const spx_word16_t *input, const spx_word16_t *e, spx_word16_t preemph, spx_word16_t mem, int len, int stride)
{
   if (SPX_MDF_RVV_ON && len >= 8)
   {
      float tmp[256];
      while (len > 0)
      {
         int i, n = len < 256 ? len : 256;
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
   return mdf_deemph_output(out, input, e, preemph, mem, len, stride);
}
#define mdf_deemph_output mdf_deemph_output_dispatch

static inline void spectral_mul_accum_dispatch(const spx_word16_t *X, const spx_word32_t *Y, spx_word16_t *acc, int N, int M)
{
   if (SPX_MDF_RVV_ON && N > 2 && !(N & 1) && M > 0)
   {
      int j;
      float a0 = 0, aN = 0;
      spx_mdf_rvv_smul_accum_f32(X, Y, acc, N, M);
      for (j=0;j<M;j++)
      {
         a0 += X[j*N]*Y[j*N];
         aN += X[(j+1)*N-1]*Y[(j+1)*N-1];
      }
      acc[0] = a0;
      acc[N-1] = aN;
   }
   else
      spectral_mul_accum(X, Y, acc, N, M);
}
#define spectral_mul_accum spectral_mul_accum_dispatch

static inline void weighted_spectral_mul_conj_dispatch(const spx_float_t *w, const spx_float_t p, const spx_word16_t *X, const spx_word16_t *Y, spx_word32_t *prod, int N)
{
   if (SPX_MDF_RVV_ON && N > 2 && !(N & 1))
   {
      spx_mdf_rvv_wsmul_conj_f32(w, X, Y, prod, N, p);
      prod[0] = p*w[0]*(X[0]*Y[0]);
      prod[N-1] = p*w[N>>1]*(X[N-1]*Y[N-1]);
   }
   else
      weighted_spectral_mul_conj(w, p, X, Y, prod, N);
}
#define weighted_spectral_mul_conj weighted_spectral_mul_conj_dispatch

#endif /* FIXED_POINT / float */

/* Wrappers whose C side is the same in both builds: the kernels differ only
 * in element type, selected through the SPX_MDF_RVV_* aliases above. */

static inline spx_word32_t mdf_inner_prod_dispatch(const spx_word16_t *x, const spx_word16_t *y, int len)
{
   /* below ~32 elements the vsetvli/reduction overhead beats the gain */
   if (SPX_MDF_RVV_ON && len >= 32)
      return SPX_MDF_RVV_INNER_PROD(x, y, len);
   return mdf_inner_prod(x, y, len);
}
#define mdf_inner_prod mdf_inner_prod_dispatch

/* Only the per-block sum of squares is vectorized; the sqrt/normalisation
 * tail of the scalar kernel is repeated here since it follows the sums. */
static inline void mdf_adjust_prop_dispatch(const spx_word32_t *W, int N, int M, int P, spx_word16_t *prop)
{
   int i, p;
   spx_word16_t max_sum = 1;
   spx_word32_t prop_sum = 1;
   if (!SPX_MDF_RVV_ON)
   {
      mdf_adjust_prop(W, N, M, P, prop);
      return;
   }
   for (i=0;i<M;i++)
   {
      spx_word32_t tmp = 1;
      for (p=0;p<P;p++)
         tmp += SPX_MDF_RVV_PROP_SUMSQ(&W[p*N*M + i*N], N);
#ifdef FIXED_POINT
      /* Just a security in case an overflow were to occur */
      tmp = MIN32(ABS32(tmp), 536870912);
#endif
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
#define mdf_adjust_prop mdf_adjust_prop_dispatch

static inline void mdf_weight_update_dispatch(spx_word32_t *w, const spx_word32_t *phi, int N)
{
   if (SPX_MDF_RVV_ON && N >= 4)
      SPX_MDF_RVV_WEIGHT_UPDATE(w, phi, N);
   else
      mdf_weight_update(w, phi, N);
}
#define mdf_weight_update mdf_weight_update_dispatch

static inline void mdf_residual_window_dispatch(spx_word16_t *y, const spx_word16_t *window, const spx_word16_t *last_y, int N)
{
   if (SPX_MDF_RVV_ON && N >= 4)
      SPX_MDF_RVV_RES_WINDOW(y, window, last_y, N);
   else
      mdf_residual_window(y, window, last_y, N);
}
#define mdf_residual_window mdf_residual_window_dispatch

static inline void mdf_residual_scale_dispatch(spx_word32_t *residual_echo, spx_word16_t leak2, int len)
{
   if (SPX_MDF_RVV_ON && len >= 4)
      SPX_MDF_RVV_RES_SCALE(residual_echo, leak2, len);
   else
      mdf_residual_scale(residual_echo, leak2, len);
}
#define mdf_residual_scale mdf_residual_scale_dispatch

static inline void power_spectrum_dispatch(const spx_word16_t *X, spx_word32_t *ps, int N)
{
   if (SPX_MDF_RVV_ON && N > 2 && !(N & 1))
   {
      ps[0]=MULT16_16(X[0],X[0]);
      SPX_MDF_RVV_POWER_SPECTRUM(X, ps, N);
      ps[N>>1]=MULT16_16(X[N-1],X[N-1]);
   }
   else
      power_spectrum(X, ps, N);
}
#define power_spectrum power_spectrum_dispatch

static inline void power_spectrum_accum_dispatch(const spx_word16_t *X, spx_word32_t *ps, int N)
{
   if (SPX_MDF_RVV_ON && N > 2 && !(N & 1))
   {
      ps[0]+=MULT16_16(X[0],X[0]);
      SPX_MDF_RVV_POWER_SPECTRUM_ACCUM(X, ps, N);
      ps[N>>1]+=MULT16_16(X[N-1],X[N-1]);
   }
   else
      power_spectrum_accum(X, ps, N);
}
#define power_spectrum_accum power_spectrum_accum_dispatch

#endif /* FIXED_POINT || __riscv_float_abi_double */

#endif /* MDF_RVV_H */
