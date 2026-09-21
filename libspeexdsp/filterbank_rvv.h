/* Copyright (C) 2006 Jean-Marc Valin
 * Copyright (C) 2026 Tristan Matthews
 */
/**
   @file filterbank_rvv.h
   @brief Filterbank psd16 kernel (RISC-V Vector extension,
          runtime-dispatched)
*/

/* Runtime-dispatched RVV kernel for filterbank_compute_psd16's per-bin
 * gather + weighted sum. The vector code lives out-of-line in
 * filterbank_rvv_asm.S, so this header is plain C and filterbank.c stays
 * base-ISA. filterbank.c includes it right after its scalar fbank_psd16:
 * the dispatch wrapper below calls the kernel when SPX_FBANK_RVV_ON and
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

#ifdef FBANK_RVV_FORCE_ON
#  define SPX_FBANK_RVV_ON 1
#else
#include "rvv_cpu.h"
extern int spx_fbank_rvv_enabled;   /* defined in filterbank.c: -1 until its first init probes, then 0/1 */
#  define SPX_FBANK_RVV_ON (spx_fbank_rvv_enabled > 0)
#  define FBANK_RVV_RUNTIME 1        /* tells filterbank.c to define+detect the flag */
#endif

void spx_fbank_rvv_psd16_f32(const int *bank_left, const int *bank_right,
                             const float *filter_left, const float *filter_right,
                             const float *mel, float *ps, int len);

static inline void fbank_psd16_dispatch(const int *bank_left, const int *bank_right,
                                        const spx_word16_t *filter_left,
                                        const spx_word16_t *filter_right,
                                        const spx_word16_t *mel, spx_word16_t *ps, int len)
{
   /* len >= 8: the strip-mine + gather setup breaks even with the scalar
    * loop at 8 elements (K1); real spectra are >= 64 bins. */
   if (SPX_FBANK_RVV_ON && len >= 8)
      spx_fbank_rvv_psd16_f32(bank_left, bank_right, filter_left, filter_right,
                              mel, ps, len);
   else
      fbank_psd16(bank_left, bank_right, filter_left, filter_right, mel, ps, len);
}
#define fbank_psd16 fbank_psd16_dispatch

#endif /* !FIXED_POINT && __riscv_float_abi_double */

#endif /* FBANK_RVV_H */
