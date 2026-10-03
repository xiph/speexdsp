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

#include <stdio.h>
#include <string.h>
#include "../internal.h"
#include "../smallft/wrap.h"

#ifdef HAVE_RVV_SMALLFT
#include "fftwrap_rvv.h"
#endif

/* fftwrap.c's USE_SMALLFT spx_fft pre-scale (out[i] = scale*in[i]),
 * C vs RVV. Both do one multiply per element, so outputs are bit-exact
 * and compared with memcmp. The n set covers the AEC shapes (2*frame =
 * 256/512), power-of-two sweeps, ragged vsetvli tails, and sub-VLMAX
 * sizes; the in-place (out == in) call spx_fft warns about is checked
 * too. The C baseline lives in wrap_smallft_c.c (no-autovec TU), the
 * kernel in smallft_rvv_asm.S (linked via the smallft RVV wrap). */

enum { MAX_N = 2048 };

static const int configs[] = { 8, 33, 64, 128, 256, 480, 512, 1024, 2048 };

void checkasm_check_fftwrap_scale(void)
{
#ifdef HAVE_RVV_SMALLFT
    checkasm_declare(void, float *, const float *, const float *, int);

    CHECKASM_ALIGN(float in[MAX_N]);
    CHECKASM_ALIGN(float out_ref[MAX_N]);
    CHECKASM_ALIGN(float out_new[MAX_N]);

    for (size_t c = 0; c < sizeof(configs) / sizeof(configs[0]); c++) {
        const int n = configs[c];
        const float scale = 1. / n;   /* as spx_fft computes it */

        smallft_fill_input(in, n);

        if (checkasm_check_func(fftwrap_scale_c, "fft_scale_n%d", n))
            checkasm_bench_new(out_new, in, &scale, n);

        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(spx_drft_rvv_scale_f32, "fft_scale_n%d", n)) {
                checkasm_call_ref(out_ref, in, &scale, n);
                checkasm_call_new(out_new, in, &scale, n);
                if (memcmp(out_ref, out_new, n * sizeof *out_ref))
                    checkasm_fail();

                memcpy(out_ref, in, n * sizeof *out_ref);
                memcpy(out_new, in, n * sizeof *out_new);
                checkasm_call_ref(out_ref, out_ref, &scale, n);
                checkasm_call_new(out_new, out_new, &scale, n);
                if (memcmp(out_ref, out_new, n * sizeof *out_ref))
                    checkasm_fail();

                checkasm_bench_new(out_new, in, &scale, n);
            }
        }
    }

    checkasm_report("fftwrap_scale");
#endif /* HAVE_RVV_SMALLFT */
}
