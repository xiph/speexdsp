#ifndef SPEEXDSP_TESTS_CHECKASM_SMALLFT_WRAP_H
#define SPEEXDSP_TESTS_CHECKASM_SMALLFT_WRAP_H

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "config.h"
#include "smallft.h"

/* Tests for the smallft (FFTPACK real FFT) dradf4/dradb4/dradf2/dradb2
 * stages (per stage and whole transform), C vs RVV. Like the fft tests,
 * each variant #includes smallft.c with its public API renamed
 * (wrap_smallft_rename.h): wrap_smallft_c.c forces the scalar stages,
 * wrap_smallft_rvv.c pins the RVV dispatch on. The kernels are plain
 * float code with no fixed/float-mode dependency, so the tests run in
 * either build mode. */

/* ------------- Per-ISA availability gates ------------- */
#ifdef USE_RVV
#  define HAVE_RVV_SMALLFT 1
#endif

/* ------------- Stage enumeration -------------
 * One entry per drftf1/drftb1 radix-2/4 stage with the arguments the
 * driver passes (dradfg/dradbg stages are skipped: no kernel). Owned by
 * wrap_smallft_c.c. */
struct drft_stage {
    int ip, ido, l1, iw;   /* iw is the 1-based FFTPACK twiddle offset */
};
int smallft_stages(struct drft_lookup *l, int backward,
                   struct drft_stage *out, int max);

/* ------------- State helpers (owned by wrap_smallft_c.c) -------------
 * trigcache/splitcache contents are identical in every TU, so one lookup
 * is shared between the C and RVV calls. */
struct drft_lookup *smallft_make_lookup(int n);
void smallft_free_lookup(struct drft_lookup *l);

/* ------------- Functions under test -------------
 * One shared signature; the shim picks the radix from st->ip and derives
 * the wa pointers from the lookup exactly as drftf1/drftb1 do. */
void smallft_stage_c(struct drft_lookup *l, const struct drft_stage *st,
                     int backward, float *cc, float *ch);
void smallft_forward_c(struct drft_lookup *l, const struct drft_stage *st,
                       int backward, float *data, float *unused);
#ifdef HAVE_RVV_SMALLFT
void smallft_stage_rvv(struct drft_lookup *l, const struct drft_stage *st,
                       int backward, float *cc, float *ch);
void smallft_forward_rvv(struct drft_lookup *l, const struct drft_stage *st,
                         int backward, float *data, float *unused);
#endif

/* ------------- Test-input fill ------------- */
#include <checkasm/utils.h>

static inline void smallft_fill_input(float *buf, int n)
{
    /* Symmetric [-1, 1) floats (raw random bits would be NaN/Inf soup). */
    checkasm_randomize_rangef(buf, n, 2.0f);
    for (int i = 0; i < n; i++)
        buf[i] -= 1.0f;
}

/* ------------- Output comparison -------------
 * The kernels use FMAs and reordered sums -> compare relative to the
 * buffer peak, tighter for a single stage than for the log(n)-deep
 * transform. */
#define DRFT_STAGE_REL_TOL 1e-5
#define DRFT_FFT_REL_TOL   1e-4

static inline int smallft_buf_within_tol(const float *ref, const float *res,
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

#endif /* SPEEXDSP_TESTS_CHECKASM_SMALLFT_WRAP_H */
