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

#ifndef SPEEXDSP_TESTS_CHECKASM_COMPARE_H
#define SPEEXDSP_TESTS_CHECKASM_COMPARE_H

#include <stdio.h>
#include <math.h>
#include "config.h"
#include "arch.h"
#include <checkasm/utils.h>

/* Input-fill and output-comparison helpers shared by the per-module tests. */

/* Symmetric [-1, 1) floats: raw random bits would be NaN/Inf soup. */
static inline void checkasm_fill_symmetric_f32(float *buf, int n)
{
    checkasm_randomize_rangef(buf, n, 2.0f);
    for (int i = 0; i < n; i++)
        buf[i] -= 1.0f;
}

/* Peak-relative (L-infinity) comparison: |ref-res| must stay within rel_tol
 * of the buffer's peak |ref|. Normalizing by the peak rather than per-sample
 * |ref| avoids the cancellation flakiness a near-zero output would otherwise
 * cause (Higham, backward error of a dot product). */
static inline int checkasm_f32_within_tol(const float *ref, const float *res,
                                          int n, double rel_tol)
{
    double peak = 0.0;
    for (int i = 0; i < n; i++) {
        double v = fabs((double) ref[i]);
        if (v > peak) peak = v;
    }
    for (int i = 0; i < n; i++) {
        double diff = fabs((double) ref[i] - (double) res[i]);
        double rel  = peak > 0.0 ? diff / peak : diff;
        /* NaN-safe: rel > rel_tol would be false for NaN and wrongly pass. */
        if (!(rel <= rel_tol)) {
            fprintf(stderr, "FAILED: [%d] ref=%g res=%g diff=%g peak=%g "
                    "rel=%.2e (tol %g)\n",
                    i, (double) ref[i], (double) res[i], diff, peak, rel, rel_tol);
            return 0;
        }
    }
    return 1;
}

static inline int checkasm_i16_bitexact(const spx_int16_t *ref, const spx_int16_t *res, int n)
{
    for (int i = 0; i < n; i++) {
        if (ref[i] != res[i]) {
            fprintf(stderr, "FAILED: [%d] ref=%d res=%d (bit-exact required)\n",
                    i, (int) ref[i], (int) res[i]);
            return 0;
        }
    }
    return 1;
}

static inline int checkasm_i32_bitexact(const spx_int32_t *ref, const spx_int32_t *res, int n)
{
    for (int i = 0; i < n; i++) {
        if (ref[i] != res[i]) {
            fprintf(stderr, "FAILED: [%d] ref=%ld res=%ld (bit-exact required)\n",
                    i, (long) ref[i], (long) res[i]);
            return 0;
        }
    }
    return 1;
}

#endif /* SPEEXDSP_TESTS_CHECKASM_COMPARE_H */
