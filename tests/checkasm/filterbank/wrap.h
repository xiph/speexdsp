/* Copyright (C) 2026 Tristan Matthews */
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

#ifndef SPEEXDSP_TESTS_CHECKASM_FILTERBANK_WRAP_H
#define SPEEXDSP_TESTS_CHECKASM_FILTERBANK_WRAP_H

#include "config.h"
#include "arch.h"
#include "filterbank.h"
#include "../compare.h"

/* Tests for filterbank_compute_psd16's per-bin two-index gather +
 * weighted sum, C vs SIMD. Like the preprocess tests, each variant
 * #includes filterbank.c with its public API renamed
 * (wrap_fbank_rename.h): wrap_fbank_c.c forces the scalar loop,
 * wrap_fbank_rvv.c pins the RVV dispatch on. The bank tables come from
 * the real filterbank_new (renamed, C TU), so the index patterns are
 * exactly the ones preprocess.c produces. */

/* ------------- Per-ISA availability gates -------------
 * Mirrors filterbank_rvv.h's gate: float only, on an FP ABI. */
#if defined(USE_RVV) && !defined(FIXED_POINT) && defined(__riscv_float_abi_double)
#  define HAVE_RVV_FBANK 1
#endif

/* ------------- Functions under test -------------
 * Thin shims around filterbank.c's fbank_psd16 inner loop; one per
 * variant TU. Only meaningful for float builds (the RVV kernel is
 * float-only), so the shims use plain float signatures. The bank
 * constructor/destructor shims expose the (renamed) real table builder
 * to the test. */
#ifndef FIXED_POINT

FilterBank *fbank_new_c(int banks, float sampling, int len, int type);
void fbank_destroy_c(FilterBank *bank);
void fbank_psd16_c(const FilterBank *bank, const float *mel, float *ps);

#ifdef HAVE_RVV_FBANK
void fbank_psd16_rvv(const FilterBank *bank, const float *mel, float *ps);
#endif

/* ------------- Test-input fill ------------- */
static inline void fbank_fill_gain(float *buf, int n) /* [0, 1) -- band gains */
{
    checkasm_randomize_rangef(buf, n, 1.0f);
}

/* ------------- Output comparison -------------
 * The RVV kernel fuses the second multiply-add (vfmacc), so results can
 * differ from the scalar two-rounding sum in the last ulp; compare with
 * checkasm_f32_within_tol relative to the buffer peak. */
#define FBANK_PSD16_F32_REL_TOL 1e-6

#endif /* !FIXED_POINT */

#endif /* SPEEXDSP_TESTS_CHECKASM_FILTERBANK_WRAP_H */
