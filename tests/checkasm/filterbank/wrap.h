#ifndef SPEEXDSP_TESTS_CHECKASM_FILTERBANK_WRAP_H
#define SPEEXDSP_TESTS_CHECKASM_FILTERBANK_WRAP_H

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "config.h"
#include "arch.h"
#include "filterbank.h"

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
#include <checkasm/utils.h>

static inline void fbank_fill_gain(float *buf, int n) /* [0, 1) -- band gains */
{
    checkasm_randomize_rangef(buf, n, 1.0f);
}

/* ------------- Output comparison -------------
 * The RVV kernel fuses the second multiply-add (vfmacc), so results can
 * differ from the scalar two-rounding sum in the last ulp; compare
 * relative to the buffer peak. */
#define FBANK_PSD16_F32_REL_TOL 1e-6

static inline int fbank_buf_within_tol(const float *ref, const float *res, int n, double rel_tol)
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

#endif /* !FIXED_POINT */

#endif /* SPEEXDSP_TESTS_CHECKASM_FILTERBANK_WRAP_H */
