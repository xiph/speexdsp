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

/* smallft radix-2/4 stages: per-stage unit test + micro-benchmark, and a
 * whole-transform integration test, C vs RVV. The n set covers the ido==1
 * fast path, both part-2 vectorization axes (across pairs and across
 * rows), the l1==1 radix-2 stage of 2*4^k sizes, and odd-ido radix-4
 * stages (n = 4*odd). Stages read cc and write every element of ch, so
 * ref/new are pre-filled with identical garbage: a lane the RVV kernel
 * misses shows up as a mismatch. */

enum { MAX_N = 2048, MAX_STAGES = 16 };

static const int configs[] = { 12, 20, 64, 128, 256, 260, 512, 1024, 2048 };

void checkasm_check_smallft_stage(void)
{
#ifdef HAVE_RVV_SMALLFT
    checkasm_declare(void, struct drft_lookup *, const struct drft_stage *,
                     int, float *, float *);

    CHECKASM_ALIGN(float in[MAX_N]);
    CHECKASM_ALIGN(float garbage[MAX_N]);
    CHECKASM_ALIGN(float out_ref[MAX_N]);
    CHECKASM_ALIGN(float out_new[MAX_N]);
    struct drft_stage st[MAX_STAGES];

    for (size_t c = 0; c < sizeof(configs) / sizeof(configs[0]); c++) {
        const int n = configs[c];
        struct drft_lookup *l = smallft_make_lookup(n);
        if (!l) { fprintf(stderr, "smallft_stage: alloc failed n%d\n", n); continue; }

        for (int backward = 0; backward <= 1; backward++) {
            const int nst = smallft_stages(l, backward, st, MAX_STAGES);
            const char *dir = backward ? "radb" : "radf";

            for (int s = 0; s < nst; s++) {
                smallft_fill_input(in, n);
                smallft_fill_input(garbage, n);

                /* C baseline, also the reference the RVV variant pairs with. */
                if (checkasm_check_func(smallft_stage_c, "drft_%s%d_n%d_ido%d_l%d",
                                        dir, st[s].ip, n, st[s].ido, st[s].l1)) {
                    checkasm_bench_new(l, &st[s], backward, in, out_new);
                }

                if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
                    if (checkasm_check_func(smallft_stage_rvv, "drft_%s%d_n%d_ido%d_l%d",
                                            dir, st[s].ip, n, st[s].ido, st[s].l1)) {
                        memcpy(out_ref, garbage, n * sizeof *out_ref);
                        memcpy(out_new, garbage, n * sizeof *out_new);
                        checkasm_call_ref(l, &st[s], backward, in, out_ref);
                        checkasm_call_new(l, &st[s], backward, in, out_new);
                        if (!checkasm_f32_within_tol(out_ref, out_new, n, DRFT_STAGE_REL_TOL))
                            checkasm_fail();
                        checkasm_bench_new(l, &st[s], backward, in, out_new);
                    }
                }
            }
        }

        smallft_free_lookup(l);
    }

    checkasm_report("smallft_stage");
#endif /* HAVE_RVV_SMALLFT */
}

void checkasm_check_smallft_transform(void)
{
#ifdef HAVE_RVV_SMALLFT
    checkasm_declare(void, struct drft_lookup *, const struct drft_stage *,
                     int, float *, float *);

    CHECKASM_ALIGN(float in[MAX_N]);
    CHECKASM_ALIGN(float buf_ref[MAX_N]);
    CHECKASM_ALIGN(float buf_new[MAX_N]);

/* The transform runs in place, so repeated benching grows the data by ~n
 * per pass and would overflow. Periodically restore it from in (tidx is
 * checkasm's bench-iteration index). */
#define BENCH_RESTORED(buf) ((tidx & 31) ? (buf) : memcpy((buf), in, (size_t) n * sizeof(float)))

    for (size_t c = 0; c < sizeof(configs) / sizeof(configs[0]); c++) {
        const int n = configs[c];
        struct drft_lookup *l = smallft_make_lookup(n);
        if (!l) { fprintf(stderr, "smallft_transform: alloc failed n%d\n", n); continue; }

        for (int backward = 0; backward <= 1; backward++) {
            const char *dir = backward ? "backward" : "forward";

            smallft_fill_input(in, n);

            if (checkasm_check_func(smallft_forward_c, "drft_%s_n%d", dir, n)) {
                memcpy(buf_new, in, n * sizeof *buf_new);
                checkasm_bench_new(l, NULL, backward, BENCH_RESTORED(buf_new), NULL);
            }

            if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
                if (checkasm_check_func(smallft_forward_rvv, "drft_%s_n%d", dir, n)) {
                    memcpy(buf_ref, in, n * sizeof *buf_ref);
                    memcpy(buf_new, in, n * sizeof *buf_new);
                    checkasm_call_ref(l, NULL, backward, buf_ref, NULL);
                    checkasm_call_new(l, NULL, backward, buf_new, NULL);
                    if (!checkasm_f32_within_tol(buf_ref, buf_new, n, DRFT_FFT_REL_TOL))
                        checkasm_fail();
                    memcpy(buf_new, in, n * sizeof *buf_new);
                    checkasm_bench_new(l, NULL, backward, BENCH_RESTORED(buf_new), NULL);
                }
            }
        }

        smallft_free_lookup(l);
    }

    checkasm_report("smallft_transform");
#undef BENCH_RESTORED
#endif /* HAVE_RVV_SMALLFT */
}
