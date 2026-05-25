/* Copyright (C) 2007-2008 Jean-Marc Valin
 * Copyright (C) 2008 Thorvald Natvig
 */
/**
   @file resample.h
   @brief Resampler functions (C version)
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

#ifndef RESAMPLE_H
#define RESAMPLE_H

#include "arch.h"

#ifndef OVERRIDE_INNER_PRODUCT_SINGLE
static inline spx_word32_t inner_product_single(const spx_word16_t *a, const spx_word16_t *b, unsigned int len)
{
   int j;
   spx_word32_t sum = 0;
   for(j=0;j<len;j++) sum += MULT16_16(a[j], b[j]);

/*    This code is slower on most DSPs which have only 2 accumulators.
      Plus this forces truncation to 32 bits and you lose the HW guard bits.
      I think we can trust the compiler and let it vectorize and/or unroll itself.
      spx_word32_t accum[4] = {0,0,0,0};
      for(j=0;j<N;j+=4) {
        accum[0] += MULT16_16(sinct[j], iptr[j]);
        accum[1] += MULT16_16(sinct[j+1], iptr[j+1]);
        accum[2] += MULT16_16(sinct[j+2], iptr[j+2]);
        accum[3] += MULT16_16(sinct[j+3], iptr[j+3]);
      }
      sum = accum[0] + accum[1] + accum[2] + accum[3];
*/
   return SATURATE32PSHR(sum, 15, 32767);
}
#endif

#ifndef OVERRIDE_INNER_PRODUCT_DOUBLE
static inline double inner_product_double(const spx_word16_t *a, const spx_word16_t *b, unsigned int len)
{
   int j;
   double accum[4] = {0,0,0,0};

   for(j=0;j<len;j+=4) {
      accum[0] += a[j]*b[j];
      accum[1] += a[j+1]*b[j+1];
      accum[2] += a[j+2]*b[j+2];
      accum[3] += a[j+3]*b[j+3];
   }
   return accum[0] + accum[1] + accum[2] + accum[3];
}
#endif

#ifndef OVERRIDE_INTERPOLATE_PRODUCT_SINGLE
static inline spx_word32_t interpolate_product_single(const spx_word16_t *a, const spx_word16_t *b, unsigned int len, const spx_uint32_t oversample, const spx_word16_t *interp)
{
   spx_word32_t sum;
   int j;
   spx_word32_t accum[4] = {0,0,0,0};

   for(j=0;j<len;j++) {
      const spx_word16_t curr_in=a[j];
      accum[0] += MULT16_16(curr_in, b[j*oversample]);
      accum[1] += MULT16_16(curr_in, b[j*oversample+1]);
      accum[2] += MULT16_16(curr_in, b[j*oversample+2]);
      accum[3] += MULT16_16(curr_in, b[j*oversample+3]);
   }

   sum = MULT16_32_Q15(interp[0],accum[0]) + MULT16_32_Q15(interp[1],accum[1]) + MULT16_32_Q15(interp[2],accum[2]) + MULT16_32_Q15(interp[3],accum[3]);
   return SATURATE32PSHR(sum, 15, 32767);
}
#endif

#ifndef OVERRIDE_INTERPOLATE_PRODUCT_DOUBLE
static inline double interpolate_product_double(const spx_word16_t *a, const spx_word16_t *b, unsigned int len, const spx_uint32_t oversample, const spx_word16_t *interp)
{
   int j;
   double sum;
   double accum[4] = {0,0,0,0};
   
   for(j=0;j<len;j++) {
      const double curr_in = a[j];
      accum[0] += curr_in * b[j*oversample];
      accum[1] += curr_in * b[j*oversample+1];
      accum[2] += curr_in * b[j*oversample+2];
      accum[3] += curr_in * b[j*oversample+3];
   }

   sum = interp[0]*accum[0] + interp[1]*accum[1] + interp[2]*accum[2] + interp[3]*accum[3];
   return sum;
}
#endif

#endif
