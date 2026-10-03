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
#include "../internal.h"
#include "wrap.h"

/* kiss_fft: whole-transform integration test + benchmark, C vs SIMD, both
 * directions -- the number spx_fft users (mdf, preprocess) actually see.
 * Bit-exact in fixed point; peak-relative tolerance in float. */

#define FUNC_NAME "kiss_fft"

enum { MAX_NFFT = 3600 };

static const int configs[] = { 60, 240, 250, 256, 480, 512, 1024, 3600 };

void checkasm_check_fft_transform(void)
{
#ifdef HAVE_RVV_KF_BFLY
    checkasm_declare(void, kiss_fft_cfg, const kiss_fft_cpx *, kiss_fft_cpx *);

    CHECKASM_ALIGN(kiss_fft_cpx in[MAX_NFFT]);
    CHECKASM_ALIGN(kiss_fft_cpx out_ref[MAX_NFFT]);
    CHECKASM_ALIGN(kiss_fft_cpx out_new[MAX_NFFT]);

    for (size_t c = 0; c < sizeof(configs) / sizeof(configs[0]); c++) {
        const int nfft = configs[c];

        for (int inverse = 0; inverse <= 1; inverse++) {
            kiss_fft_cfg cfg = fft_make_cfg(nfft, inverse);
            if (!cfg) { fprintf(stderr, FUNC_NAME ": alloc failed n%d\n", nfft); continue; }
            const char *dir = inverse ? "inv" : "fwd";

            fft_fill_input(in, nfft);

            if (checkasm_check_func(fft_c, FUNC_NAME "_%s_n%d", dir, nfft))
                checkasm_bench_new(cfg, in, out_new);

            if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
                if (checkasm_check_func(fft_rvv, FUNC_NAME "_%s_n%d", dir, nfft)) {
                    checkasm_call_ref(cfg, in, out_ref);
                    checkasm_call_new(cfg, in, out_new);
                    if (!fft_buf_matches(out_ref, out_new, nfft, KF_FFT_F32_REL_TOL))
                        checkasm_fail();
                    checkasm_bench_new(cfg, in, out_new);
                }
            }

            fft_free_cfg(cfg);
        }
    }

    checkasm_report(FUNC_NAME);
#endif /* HAVE_RVV_KF_BFLY */
}
