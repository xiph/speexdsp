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
