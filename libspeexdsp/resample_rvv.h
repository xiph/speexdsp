/* Copyright (C) 2007-2008 Jean-Marc Valin
 * Copyright (C) 2008 Thorvald Natvig
 * Copyright (C) 2026 Tyler Kim
 */
/**
   @file resample_rvv.h
   @brief Resampler inner products (RISC-V Vector extension, runtime-dispatched)
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

/* Runtime-dispatched RVV inner products. The vector kernels live out-of-line in
 * resample_rvv_asm.S, so this header is plain C and including it keeps
 * resample.c base-ISA. Each OVERRIDE_* leaf calls a kernel only when
 * spx_rvv_enabled, else runs the stock scalar loop -- so one build runs on both
 * V and non-V CPUs.
 *
 * The float kernels return in fa0 (lp64d/ilp32d hard-double), so they bind only
 * under __riscv_float_abi_double; other float ABIs get the scalar path. The
 * int16 kernels return integers and work under every ABI.
 *
 * checkasm defines RESAMPLE_RVV_FORCE_ON to test the asm unconditionally. */

#ifdef USE_RVV_ASM

#if defined(FIXED_POINT) || defined(__riscv_float_abi_double)
#  ifdef RESAMPLE_RVV_FORCE_ON
#    define SPX_RVV_ON 1
#  else
extern int spx_rvv_enabled;          /* defined in resample.c, detected at init */
unsigned int spx_resample_rvv_vlenb(void);
int spx_resample_rvv_compliant(void);
#    define SPX_RVV_ON spx_rvv_enabled
#    define RESAMPLE_RVV_RUNTIME 1    /* tells resample.c to define+detect the flag */
#  endif
#endif

#ifdef FIXED_POINT

spx_int32_t spx_resample_rvv_ip_i16(const spx_int16_t *a, const spx_int16_t *b,
                                    unsigned int len);
void spx_resample_rvv_interp4_i16(const spx_int16_t *a, const spx_int16_t *b,
                                  unsigned int len, unsigned int oversample,
                                  spx_int32_t *accum);

#define OVERRIDE_INNER_PRODUCT_SINGLE
static inline spx_word32_t inner_product_single(const spx_word16_t *a, const spx_word16_t *b, unsigned int len)
{
   spx_word32_t sum;
   if (SPX_RVV_ON)
      sum = spx_resample_rvv_ip_i16(a, b, len);
   else
   {
      unsigned int j;
      sum = 0;
      for (j = 0; j < len; j++)
         sum += MULT16_16(a[j], b[j]);
   }
   return SATURATE32PSHR(sum, 15, 32767);
}

#define OVERRIDE_INTERPOLATE_PRODUCT_SINGLE
static inline spx_word32_t interpolate_product_single(const spx_word16_t *a, const spx_word16_t *b, unsigned int len, const spx_uint32_t oversample, const spx_word16_t *frac)
{
   spx_int32_t accum[4];
   spx_word32_t sum;
   if (SPX_RVV_ON)
      spx_resample_rvv_interp4_i16(a, b, len, oversample, accum);
   else
   {
      unsigned int j;
      accum[0] = accum[1] = accum[2] = accum[3] = 0;
      for (j = 0; j < len; j++)
      {
         const spx_word16_t curr = a[j];
         const spx_word16_t *bi = b + j * oversample;
         accum[0] += MULT16_16(curr, bi[0]);
         accum[1] += MULT16_16(curr, bi[1]);
         accum[2] += MULT16_16(curr, bi[2]);
         accum[3] += MULT16_16(curr, bi[3]);
      }
   }
   sum = MULT16_32_Q15(frac[0], accum[0]) + MULT16_32_Q15(frac[1], accum[1])
       + MULT16_32_Q15(frac[2], accum[2]) + MULT16_32_Q15(frac[3], accum[3]);
   return SATURATE32PSHR(sum, 15, 32767);
}

#elif defined(FLOATING_POINT) && defined(__riscv_float_abi_double)

float  spx_resample_rvv_ip_f32(const float *a, const float *b, unsigned int len);
double spx_resample_rvv_ipd_f32(const float *a, const float *b, unsigned int len);
float  spx_resample_rvv_interp_f32(const float *a, const float *b, unsigned int len,
                                   unsigned int oversample, const float *frac);
double spx_resample_rvv_interpd_f32(const float *a, const float *b, unsigned int len,
                                    unsigned int oversample, const float *frac);

#define OVERRIDE_INNER_PRODUCT_SINGLE
static inline float inner_product_single(const float *a, const float *b, unsigned int len)
{
   if (SPX_RVV_ON)
      return spx_resample_rvv_ip_f32(a, b, len);
   {
      float sum = 0;
      unsigned int j;
      for (j = 0; j < len; j++)
         sum += a[j] * b[j];
      return sum;
   }
}

#define OVERRIDE_INNER_PRODUCT_DOUBLE
static inline double inner_product_double(const float *a, const float *b, unsigned int len)
{
   if (SPX_RVV_ON)
      return spx_resample_rvv_ipd_f32(a, b, len);
   {
      double accum[4] = {0, 0, 0, 0};
      unsigned int j;
      for (j = 0; j < len; j += 4)
      {
         accum[0] += a[j]   * b[j];
         accum[1] += a[j+1] * b[j+1];
         accum[2] += a[j+2] * b[j+2];
         accum[3] += a[j+3] * b[j+3];
      }
      return accum[0] + accum[1] + accum[2] + accum[3];
   }
}

#define OVERRIDE_INTERPOLATE_PRODUCT_SINGLE
static inline float interpolate_product_single(const float *a, const float *b, unsigned int len, const spx_uint32_t oversample, float *frac)
{
   if (SPX_RVV_ON)
      return spx_resample_rvv_interp_f32(a, b, len, oversample, frac);
   {
      float accum[4] = {0, 0, 0, 0};
      unsigned int j;
      for (j = 0; j < len; j++)
      {
         const float curr = a[j];
         const float *bi = b + j * oversample;
         accum[0] += curr * bi[0];
         accum[1] += curr * bi[1];
         accum[2] += curr * bi[2];
         accum[3] += curr * bi[3];
      }
      return frac[0] * accum[0] + frac[1] * accum[1]
           + frac[2] * accum[2] + frac[3] * accum[3];
   }
}

#define OVERRIDE_INTERPOLATE_PRODUCT_DOUBLE
static inline double interpolate_product_double(const float *a, const float *b, unsigned int len, const spx_uint32_t oversample, float *frac)
{
   if (SPX_RVV_ON)
      return spx_resample_rvv_interpd_f32(a, b, len, oversample, frac);
   {
      double accum[4] = {0, 0, 0, 0};
      unsigned int j;
      for (j = 0; j < len; j++)
      {
         const float curr = a[j];
         const float *bi = b + j * oversample;
         accum[0] += MULT16_16(curr, bi[0]);
         accum[1] += MULT16_16(curr, bi[1]);
         accum[2] += MULT16_16(curr, bi[2]);
         accum[3] += MULT16_16(curr, bi[3]);
      }
      return frac[0] * accum[0] + frac[1] * accum[1]
           + frac[2] * accum[2] + frac[3] * accum[3];
   }
}

#endif /* FIXED_POINT / FLOATING_POINT+lp64d */

#endif /* USE_RVV_ASM */
