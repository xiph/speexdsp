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
   modification, are permitted provided that the following conditions are
   met:

   1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
   IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
   OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
   HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
   STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

/* Runtime-dispatched RVV kernel for filterbank_compute_psd16's per-bin
 * gather + weighted sum. The vector code lives out-of-line in
 * filterbank_rvv_asm.S, so this header is plain C and filterbank.c stays
 * base-ISA. filterbank.c includes it right after its scalar fbank_psd16:
 * the dispatch wrapper below calls the kernel when the bank's arch allows and
 * the scalar loop otherwise, and the #define points the call site at it.
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

#include "cpu_support.h"

/* Passed to the dispatch wrapper by the redirect macro below: the calling
 * function's bank->arch, or SPX_ARCH_RVV under checkasm. */
#ifdef FBANK_RVV_FORCE_ON
#  define SPX_FBANK_ARCH SPX_ARCH_RVV
#else
#  define SPX_FBANK_ARCH (bank->arch)
#endif

void spx_fbank_rvv_psd16_f32(const int *bank_left, const int *bank_right,
                             const float *filter_left, const float *filter_right,
                             const float *mel, float *ps, int len);

static inline void fbank_psd16_dispatch(int arch, const int *bank_left, const int *bank_right,
                                        const spx_word16_t *filter_left,
                                        const spx_word16_t *filter_right,
                                        const spx_word16_t *mel, spx_word16_t *ps, int len)
{
   /* len >= 8: the strip-mine + gather setup breaks even with the scalar
    * loop at 8 elements (K1); real spectra are >= 64 bins. */
   if (arch >= SPX_ARCH_RVV && len >= 8)
      spx_fbank_rvv_psd16_f32(bank_left, bank_right, filter_left, filter_right,
                              mel, ps, len);
   else
      fbank_psd16(bank_left, bank_right, filter_left, filter_right, mel, ps, len);
}
#define fbank_psd16(bank_left, bank_right, filter_left, filter_right, mel, ps, len) \
   fbank_psd16_dispatch(SPX_FBANK_ARCH, bank_left, bank_right, filter_left, filter_right, mel, ps, len)

#endif /* !FIXED_POINT && __riscv_float_abi_double */

#endif /* FBANK_RVV_H */
