#ifndef SPEEXDSP_TESTS_CHECKASM_SMALLFT_WRAP_H
#define SPEEXDSP_TESTS_CHECKASM_SMALLFT_WRAP_H

#include <string.h>
#include "config.h"
#include "smallft.h"
#include "../compare.h"

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

/* ------------- fftwrap forward-FFT pre-scale -------------
 * out[i] = *scale * in[i]; C baseline owned by wrap_smallft_c.c (the
 * no-autovec TU), tested against smallft_rvv_asm.S's
 * spx_drft_rvv_scale_f32 by checkasm/fftwrap/fftwrap_scale.c. */
void fftwrap_scale_c(float *out, const float *in, const float *scale, int n);

/* ------------- Test-input fill ------------- */
#define smallft_fill_input checkasm_fill_symmetric_f32

/* ------------- Output comparison -------------
 * The kernels use FMAs and reordered sums -> compare with
 * checkasm_f32_within_tol relative to the buffer peak, tighter for a
 * single stage than for the log(n)-deep transform. */
#define DRFT_STAGE_REL_TOL 1e-5
#define DRFT_FFT_REL_TOL   1e-4

#endif /* SPEEXDSP_TESTS_CHECKASM_SMALLFT_WRAP_H */
