/* Copyright (C) 2003-2004 Mark Borgerding
 * Copyright (C) 2005-2007 Jean-Marc Valin
 * Copyright (C) 2026 Tristan Matthews
 */
/**
   @file kiss_fft_rvv.h
   @brief KISS FFT butterflies (RISC-V Vector extension, runtime-dispatched)
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

/* Runtime-dispatched RVV butterflies for kf_bfly2/3/4/5. The vector kernels
 * live out-of-line in kiss_fft_rvv_asm.S, so this header is plain C and
 * kiss_fft.c stays base-ISA; each kf_bfly* falls back to its scalar loop
 * unless SPX_KF_RVV_ON. Kernels take only pointer/integer args (no float-ABI
 * dependency). Fixed point is bit-exact vs C; float pairs with a checkasm
 * tolerance. checkasm defines KISS_FFT_RVV_FORCE_ON to test the asm
 * unconditionally. */

#ifndef KISS_FFT_RVV_H
#define KISS_FFT_RVV_H

#include "kiss_fft.h"

#ifdef KISS_FFT_RVV_FORCE_ON
#  define SPX_KF_RVV_ON 1
#else
extern int spx_kf_rvv_enabled;       /* defined in kiss_fft.c, detected at alloc */
unsigned int spx_kf_rvv_vlenb(void);
int spx_kf_rvv_compliant(void);
#  define SPX_KF_RVV_ON spx_kf_rvv_enabled
#  define KISS_FFT_RVV_RUNTIME 1     /* tells kiss_fft.c to define+detect the flag */
#endif

/* bfly2/bfly4 (and float bfly5) carry kf_work's batching parameters (N
 * sub-FFTs, mm apart); bfly3 and fixed-point bfly5 are called once per
 * sub-FFT like their C counterparts. SPX_KF_RVV_BFLY5_BATCH tells
 * kiss_fft.c which bfly5 calling convention the kernel uses; the batched
 * kernel also handles m==1 (across sub-FFTs), the per-sub-FFT one does not. */
#ifdef FIXED_POINT

void spx_kf_rvv_bfly2_i16(kiss_fft_cpx *Fout, const kiss_fft_cpx *tw,
                          size_t fstride, int m, int N, int mm, int inverse);
void spx_kf_rvv_bfly3_i16(kiss_fft_cpx *Fout, const kiss_fft_cpx *tw,
                          size_t fstride, int m, int inverse);
void spx_kf_rvv_bfly4_i16(kiss_fft_cpx *Fout, const kiss_fft_cpx *tw,
                          size_t fstride, int m, int N, int mm, int inverse);
void spx_kf_rvv_bfly5_i16(kiss_fft_cpx *Fout, const kiss_fft_cpx *tw,
                          size_t fstride, int m, int inverse);

#define SPX_KF_RVV_BFLY2 spx_kf_rvv_bfly2_i16
#define SPX_KF_RVV_BFLY3 spx_kf_rvv_bfly3_i16
#define SPX_KF_RVV_BFLY4 spx_kf_rvv_bfly4_i16
#define SPX_KF_RVV_BFLY5 spx_kf_rvv_bfly5_i16
#define SPX_KF_RVV_BFLY5_BATCH 0

#elif !defined(USE_SIMD) /* float kiss_fft_scalar */

void spx_kf_rvv_bfly2_f32(kiss_fft_cpx *Fout, const kiss_fft_cpx *tw,
                          size_t fstride, int m, int N, int mm, int inverse);
void spx_kf_rvv_bfly3_f32(kiss_fft_cpx *Fout, const kiss_fft_cpx *tw,
                          size_t fstride, int m, int inverse);
void spx_kf_rvv_bfly4_f32(kiss_fft_cpx *Fout, const kiss_fft_cpx *tw,
                          size_t fstride, int m, int N, int mm, int inverse);
void spx_kf_rvv_bfly5_f32(kiss_fft_cpx *Fout, const kiss_fft_cpx *tw,
                          size_t fstride, int m, int N, int mm, int inverse);

#define SPX_KF_RVV_BFLY2 spx_kf_rvv_bfly2_f32
#define SPX_KF_RVV_BFLY3 spx_kf_rvv_bfly3_f32
#define SPX_KF_RVV_BFLY4 spx_kf_rvv_bfly4_f32
#define SPX_KF_RVV_BFLY5 spx_kf_rvv_bfly5_f32
#define SPX_KF_RVV_BFLY5_BATCH 1

#endif /* FIXED_POINT / float */

#endif /* KISS_FFT_RVV_H */
