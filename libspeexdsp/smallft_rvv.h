/* Copyright (C) 1994-1996 Paul N. Swarztrauber (FFTPACK)
 * Copyright (C) The OggSquish source code, Xiph.Org
 * Copyright (C) 2026 Tristan Matthews
 */
/**
   @file smallft_rvv.h
   @brief smallft (FFTPACK real FFT) radix stages (RISC-V Vector extension,
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

/* Runtime-dispatched RVV kernels for the dradf4/dradb4 (and l1==1
 * dradf2/dradb2) stages that power-of-two smallft transforms spend all
 * their time in. The vector code lives out-of-line in smallft_rvv_asm.S,
 * so this header is plain C and smallft.c stays base-ISA. The dradN_rvv
 * wrappers dispatch per call and fall back to smallft.c's scalar stages
 * (forward-declared here, defined later in the same TU); drftf1/drftb1
 * call them through the SPX_DRAD* macros. The kernels reorder float sums
 * and pair with a checkasm tolerance. checkasm defines
 * SMALLFT_RVV_FORCE_ON to test the asm unconditionally. */

#ifndef SMALLFT_RVV_H
#define SMALLFT_RVV_H

#ifdef SMALLFT_RVV_FORCE_ON
#  define SPX_DRFT_RVV_ON 1
#else
extern int spx_drft_rvv_enabled;     /* defined in smallft.c, detected at init */
unsigned int spx_drft_rvv_vlenb(void);
int spx_drft_rvv_compliant(void);
#  define SPX_DRFT_RVV_ON spx_drft_rvv_enabled
#  define SMALLFT_RVV_RUNTIME 1      /* tells smallft.c to define+detect the flag */
#endif

void spx_drft_rvv_radf4_f32(int ido, int l1, const float *cc, float *ch,
                            const float *wa1, const float *wa2,
                            const float *wa3);
void spx_drft_rvv_radb4_f32(int ido, int l1, const float *cc, float *ch,
                            const float *wa1, const float *wa2,
                            const float *wa3);
void spx_drft_rvv_radf2_f32(int ido, const float *cc, float *ch,
                            const float *wa1);
void spx_drft_rvv_radb2_f32(int ido, const float *cc, float *ch,
                            const float *wa1);

/* smallft.c's scalar stages, defined after this header is included. */
static void dradf2(int ido, int l1, float *cc, float *ch, float *wa1);
static void dradf4(int ido, int l1, float *cc, float *ch, float *wa1,
                   float *wa2, float *wa3);
static void dradb2(int ido, int l1, float *cc, float *ch, float *wa1);
static void dradb4(int ido, int l1, float *cc, float *ch, float *wa1,
                   float *wa2, float *wa3);

/* Small twiddled shapes (4*ido*l1 < 256, i.e. n <= 128 transforms' interior
 * stages) lose to call/strided overhead on real hardware; the ido==1 fast
 * path always wins. */
static inline void dradf4_rvv(int ido, int l1, float *cc, float *ch,
                              float *wa1, float *wa2, float *wa3)
{
   if (SPX_DRFT_RVV_ON && (ido == 1 || 4*ido*l1 >= 256))
      spx_drft_rvv_radf4_f32(ido, l1, cc, ch, wa1, wa2, wa3);
   else
      dradf4(ido, l1, cc, ch, wa1, wa2, wa3);
}

static inline void dradb4_rvv(int ido, int l1, float *cc, float *ch,
                              float *wa1, float *wa2, float *wa3)
{
   if (SPX_DRFT_RVV_ON && (ido == 1 || 4*ido*l1 >= 256))
      spx_drft_rvv_radb4_f32(ido, l1, cc, ch, wa1, wa2, wa3);
   else
      dradb4(ido, l1, cc, ch, wa1, wa2, wa3);
}

/* drfti1 moves the (single) factor 2 to the front of ifac, so the radix-2
 * stage always runs with l1==1; the kernel only implements that shape. */
static inline void dradf2_rvv(int ido, int l1, float *cc, float *ch,
                              float *wa1)
{
   if (SPX_DRFT_RVV_ON && l1 == 1 && ido >= 32)
      spx_drft_rvv_radf2_f32(ido, cc, ch, wa1);
   else
      dradf2(ido, l1, cc, ch, wa1);
}

static inline void dradb2_rvv(int ido, int l1, float *cc, float *ch,
                              float *wa1)
{
   if (SPX_DRFT_RVV_ON && l1 == 1 && ido >= 32)
      spx_drft_rvv_radb2_f32(ido, cc, ch, wa1);
   else
      dradb2(ido, l1, cc, ch, wa1);
}

#define SPX_DRADF2 dradf2_rvv
#define SPX_DRADF4 dradf4_rvv
#define SPX_DRADB2 dradb2_rvv
#define SPX_DRADB4 dradb4_rvv

#endif /* SMALLFT_RVV_H */
