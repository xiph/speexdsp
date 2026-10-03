/* Copyright (C) 2026 Tyler Kim
 * Copyright (C) 2026 Tristan Matthews
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

#include <stdio.h>
#include "../internal.h"
#include "wrap.h"

/* resampler_basic_interpolate_single: interpolated sinc table, single-precision
 * accumulator. Inner kernel is interpolate_product_single -> SIMD via SSE or
 * NEON (both float only; no fixed-point override) or RVV (both modes; bit-exact
 * in fixed point).
 *
 * Parameterized over several coprime-ratio (large den_rate) conversions that all
 * select the interpolate-single kernel, with different qualities/filter lengths
 * so the kernel is exercised over a range of inner-loop lengths. Each config
 * asserts it landed in the expected kind and is skipped (with a note) otherwise. */

#define FUNC_NAME "resampler_interpolate_product_single"

enum { OUT_LEN = 256, IN_LEN = 2048, IN_BUF = 4096 };

void test_resampler_basic_interpolate_single(void)
{
#if (defined(USE_SSE) && !defined(FIXED_POINT)) || defined(HAVE_NEON_INTERPOLATE_SINGLE) || defined(HAVE_RVV_INTERPOLATE_SINGLE)
    static const struct { unsigned in_rate, out_rate; int quality; } configs[] = {
        { 44100, 48000, 5 },   /* 147:160 up   */
        { 48000, 44100, 5 },   /* 160:147 down */
        { 44100, 32000, 5 },   /* 441:320 down */
        { 16000, 44100, 5 },   /* 160:441 up   */
        { 44100, 48000, 3 },   /* lower quality -> shorter filter */
    };

    checkasm_declare(int, SpeexResamplerState *, spx_uint32_t, const spx_word16_t *,
                     spx_uint32_t *, spx_word16_t *, spx_uint32_t *);

    CHECKASM_ALIGN(spx_word16_t in[IN_BUF]);
    CHECKASM_ALIGN(spx_word16_t out_ref[OUT_LEN]);
    CHECKASM_ALIGN(spx_word16_t out_new[OUT_LEN]);
    resample_fill_input(in, IN_BUF);

    spx_uint32_t inl, outl;

    for (size_t i = 0; i < sizeof(configs) / sizeof(configs[0]); i++) {
        const unsigned ir = configs[i].in_rate, orr = configs[i].out_rate;
        const int q = configs[i].quality;

        SpeexResamplerState *st = resample_make_state_c(ir, orr, q);
        if (!st) { fprintf(stderr, "%s: init failed %u->%u q%d\n", FUNC_NAME, ir, orr, q); continue; }
        if (resample_kind(st) != RESAMPLE_KIND_INTERPOLATE_SINGLE) {
            fprintf(stderr, "%s: %u->%u q%d unexpected kind %d (filt_len %u), skipping\n",
                    FUNC_NAME, ir, orr, q, resample_kind(st), resample_filt_len(st));
            resample_destroy_state(st);
            continue;
        }

        /* C baseline (also the reference the SIMD variants are paired against). */
        if (checkasm_check_func(resampler_basic_interpolate_single_c, FUNC_NAME "_%u_%u_q%d", ir, orr, q)) {
            inl = IN_LEN; outl = OUT_LEN;
            checkasm_bench_new(st, 0, in, &inl, out_new, &outl);
        }

#if defined(USE_SSE) && !defined(FIXED_POINT)
        if (active_flags & SPEEXDSP_CPU_FLAG_SSE) {
            if (checkasm_check_func(resampler_basic_interpolate_single_sse, FUNC_NAME "_%u_%u_q%d", ir, orr, q)) {
                inl = IN_LEN; outl = OUT_LEN;
                int rc = checkasm_call_ref(st, 0, in, &inl, out_ref, &outl);
                inl = IN_LEN; outl = OUT_LEN;
                int rn = checkasm_call_new(st, 0, in, &inl, out_new, &outl);
                if (rc != rn || !resample_buffer_within_tol(out_ref, out_new, (unsigned) rc, RESAMPLE_FLOAT_REL_TOL))
                    checkasm_fail();
                inl = IN_LEN; outl = OUT_LEN;
                checkasm_bench_new(st, 0, in, &inl, out_new, &outl);
            }
        }
#endif

#ifdef HAVE_NEON_INTERPOLATE_SINGLE
        if (active_flags & SPEEXDSP_CPU_FLAG_NEON) {
            if (checkasm_check_func(resampler_basic_interpolate_single_neon, FUNC_NAME "_%u_%u_q%d", ir, orr, q)) {
                inl = IN_LEN; outl = OUT_LEN;
                int rc = checkasm_call_ref(st, 0, in, &inl, out_ref, &outl);
                inl = IN_LEN; outl = OUT_LEN;
                int rn = checkasm_call_new(st, 0, in, &inl, out_new, &outl);
                if (rc != rn || !resample_buffer_within_tol(out_ref, out_new, (unsigned) rc, RESAMPLE_FLOAT_REL_TOL))
                    checkasm_fail();
                inl = IN_LEN; outl = OUT_LEN;
                checkasm_bench_new(st, 0, in, &inl, out_new, &outl);
            }
        }
#endif

#ifdef HAVE_RVV_INTERPOLATE_SINGLE
        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(resampler_basic_interpolate_single_rvv, FUNC_NAME "_%u_%u_q%d", ir, orr, q)) {
                inl = IN_LEN; outl = OUT_LEN;
                int rc = checkasm_call_ref(st, 0, in, &inl, out_ref, &outl);
                inl = IN_LEN; outl = OUT_LEN;
                int rn = checkasm_call_new(st, 0, in, &inl, out_new, &outl);
                /* fixed-point RVV is bit-exact vs C (see wrap.h) */
#ifdef FIXED_POINT
                if (rc != rn || !resample_buffer_bitexact(out_ref, out_new, (unsigned) rc))
#else
                if (rc != rn || !resample_buffer_within_tol(out_ref, out_new, (unsigned) rc, RESAMPLE_FLOAT_REL_TOL))
#endif
                    checkasm_fail();
                inl = IN_LEN; outl = OUT_LEN;
                checkasm_bench_new(st, 0, in, &inl, out_new, &outl);
            }
        }
#endif

        resample_destroy_state(st);
    }

    checkasm_report(FUNC_NAME);
#endif /* any SIMD available */
}
