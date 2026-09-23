/* Copyright (C) 2026 Tristan Matthews */
/**
   @file fftwrap_rvv.h
   @brief fftwrap's forward-FFT scale loop (RISC-V Vector extension,
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

/* fftwrap.c's USE_SMALLFT spx_fft pre-scales its input by 1/n before
 * spx_drft_forward; that loop is per-sample scalar work on every forward
 * FFT. The kernel lives in smallft_rvv_asm.S (always linked alongside the
 * smallft RVV stages) and dispatches on the drft lookup's arch level.
 * Bit-exact vs the scalar loop (one multiply per element). */

#ifndef FFTWRAP_RVV_H
#define FFTWRAP_RVV_H

#include "cpu_support.h"

/* out[i] = *scale * in[i]; out == in is allowed. scale is passed by
 * reference to keep the asm float-ABI-independent. */
void spx_drft_rvv_scale_f32(float *out, const float *in,
                            const float *scale, int n);

#endif /* FFTWRAP_RVV_H */
