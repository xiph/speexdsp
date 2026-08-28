/* Copyright (C) 2006 Jean-Marc Valin
 * Copyright (C) 2026 Tristan Matthews
 */
/**
   @file filterbank_rvv.h
   @brief Filterbank psd16 kernel (RISC-V Vector extension,
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

/* Runtime-dispatched RVV kernel for filterbank_compute_psd16's per-bin
 * gather + weighted sum. The vector code lives out-of-line in
 * filterbank_rvv_asm.S, so this header is plain C and filterbank.c stays
 * base-ISA; the OVERRIDE_FBANK_PSD16 wrapper calls the kernel only when
 * SPX_FBANK_RVV_ON and falls back to the scalar loop otherwise.
 *
 * Float only: matches the preprocess RVV set (psd16's only caller) and
 * the profiled float pipeline. The fixed-point twin would be a
 * vwmul/vwmacc + rounding-narrow port -- straightforward if it ever
 * shows up in a fixed-point profile.
 *
 * The kernel fuses the second multiply-add (vfmacc), so results can
 * differ from the scalar two-rounding sum in the last ulp; checkasm
 * compares with the standard elementwise relative tolerance.
 *
 * checkasm defines FBANK_RVV_FORCE_ON to test the asm unconditionally. */

#ifndef FBANK_RVV_H
#define FBANK_RVV_H

#include "arch.h"

#if !defined(FIXED_POINT) && defined(__riscv_float_abi_double)

#ifdef FBANK_RVV_FORCE_ON
#  define SPX_FBANK_RVV_ON 1
#else
extern int spx_fbank_rvv_enabled;  /* defined in filterbank.c, detected at init */
unsigned int spx_fbank_rvv_vlenb(void);
int spx_fbank_rvv_compliant(void);
#  define SPX_FBANK_RVV_ON spx_fbank_rvv_enabled
#  define FBANK_RVV_RUNTIME 1        /* tells filterbank.c to define+detect the flag */
#endif

void spx_fbank_rvv_psd16_f32(const int *bank_left, const int *bank_right,
                             const float *filter_left, const float *filter_right,
                             const float *mel, float *ps, int len);

#define OVERRIDE_FBANK_PSD16
static inline void fbank_psd16(const int *bank_left, const int *bank_right,
                               const spx_word16_t *filter_left,
                               const spx_word16_t *filter_right,
                               const spx_word16_t *mel, spx_word16_t *ps, int len)
{
   int i;
   /* len >= 8: the strip-mine + gather setup breaks even with the scalar
    * loop at 8 elements (K1); real spectra are >= 64 bins. */
   if (SPX_FBANK_RVV_ON && len >= 8)
   {
      spx_fbank_rvv_psd16_f32(bank_left, bank_right, filter_left, filter_right,
                              mel, ps, len);
      return;
   }
   for (i=0;i<len;i++)
   {
      spx_word32_t tmp;
      tmp = MULT16_16(mel[bank_left[i]],filter_left[i]);
      tmp += MULT16_16(mel[bank_right[i]],filter_right[i]);
      ps[i] = EXTRACT16(PSHR32(tmp,15));
   }
}

#endif /* !FIXED_POINT && __riscv_float_abi_double */

#endif /* FBANK_RVV_H */
