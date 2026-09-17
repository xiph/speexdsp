#ifndef SPEEXDSP_TESTS_CHECKASM_MDF_WRAP_H
#define SPEEXDSP_TESTS_CHECKASM_MDF_WRAP_H

#include <string.h>
#include "config.h"
#include "arch.h"
#include "pseudofloat.h"
#include "../compare.h"

/* Tests for mdf.c's spectral kernels (spectral_mul_accum(16),
 * weighted_spectral_mul_conj, power_spectrum(_accum), mdf_inner_prod,
 * mdf_adjust_prop), C vs SIMD. Like the fft tests, each variant #includes
 * mdf.c with its public API renamed (wrap_mdf_rename.h): wrap_mdf_c.c
 * forces the scalar loops, wrap_mdf_rvv.c pins the RVV dispatch on. */

/* ------------- Per-ISA availability gates -------------
 * Mirrors mdf_rvv.h's gate: fixed point always, float only on an FP ABI. */
#if defined(USE_RVV) && (defined(FIXED_POINT) || defined(__riscv_float_abi_double))
#  define HAVE_RVV_MDF 1
#endif

/* ------------- Functions under test -------------
 * Thin shims around mdf.c's static inline kernels; one per variant TU.
 * weighted_spectral_mul_conj's p scalar is passed by pointer so the shim
 * signature stays uniform across fixed/float spx_float_t. */
void mdf_smul_accum_c(const spx_word16_t *X, const spx_word32_t *Y,
                      spx_word16_t *acc, int N, int M);
void mdf_smul_accum16_c(const spx_word16_t *X, const spx_word16_t *Y,
                        spx_word16_t *acc, int N, int M);
void mdf_wsmul_conj_c(const spx_float_t *w, const spx_float_t *p,
                      const spx_word16_t *X, const spx_word16_t *Y,
                      spx_word32_t *prod, int N);
void mdf_power_spectrum_c(const spx_word16_t *X, spx_word32_t *ps, int N);
void mdf_power_spectrum_accum_c(const spx_word16_t *X, spx_word32_t *ps, int N);
spx_word32_t mdf_inner_prod_c(const spx_word16_t *x, const spx_word16_t *y, int len);
void mdf_adjust_prop_c(const spx_word32_t *W, int N, int M, int P,
                       spx_word16_t *prop);
void mdf_res_window_c(spx_word16_t *y, const spx_word16_t *window,
                      const spx_word16_t *last_y, int N);
/* leak2 by pointer, same convention as wsmul's p */
void mdf_res_scale_c(spx_word32_t *residual_echo, const spx_word16_t *leak2, int len);
void mdf_weight_update_c(spx_word32_t *w, const spx_word32_t *phi, int N);
/* preemph/mem by pointer, same convention; returns the final mem */
spx_word16_t mdf_deemph_output_c(spx_int16_t *out, const spx_word16_t *input,
                                 const spx_word16_t *e, const spx_word16_t *preemph,
                                 const spx_word16_t *mem, int len, int stride);

#ifdef HAVE_RVV_MDF
void mdf_smul_accum_rvv(const spx_word16_t *X, const spx_word32_t *Y,
                        spx_word16_t *acc, int N, int M);
void mdf_smul_accum16_rvv(const spx_word16_t *X, const spx_word16_t *Y,
                          spx_word16_t *acc, int N, int M);
void mdf_wsmul_conj_rvv(const spx_float_t *w, const spx_float_t *p,
                        const spx_word16_t *X, const spx_word16_t *Y,
                        spx_word32_t *prod, int N);
void mdf_power_spectrum_rvv(const spx_word16_t *X, spx_word32_t *ps, int N);
void mdf_power_spectrum_accum_rvv(const spx_word16_t *X, spx_word32_t *ps, int N);
spx_word32_t mdf_inner_prod_rvv(const spx_word16_t *x, const spx_word16_t *y, int len);
void mdf_adjust_prop_rvv(const spx_word32_t *W, int N, int M, int P,
                         spx_word16_t *prop);
void mdf_res_window_rvv(spx_word16_t *y, const spx_word16_t *window,
                        const spx_word16_t *last_y, int N);
void mdf_res_scale_rvv(spx_word32_t *residual_echo, const spx_word16_t *leak2, int len);
void mdf_weight_update_rvv(spx_word32_t *w, const spx_word32_t *phi, int N);
spx_word16_t mdf_deemph_output_rvv(spx_int16_t *out, const spx_word16_t *input,
                                   const spx_word16_t *e, const spx_word16_t *preemph,
                                   const spx_word16_t *mem, int len, int stride);
#endif

/* ------------- Test-input fill ------------- */
static inline void mdf_fill_w16(spx_word16_t *buf, int n)
{
#ifdef FIXED_POINT
    /* Full-range random int16: bit-exactness needs no headroom. */
    checkasm_init(buf, (size_t) n * sizeof *buf);
#else
    checkasm_fill_symmetric_f32(buf, n);
#endif
}

static inline void mdf_fill_w32(spx_word32_t *buf, int n)
{
#ifdef FIXED_POINT
    checkasm_init(buf, (size_t) n * sizeof *buf);
#else
    mdf_fill_w16(buf, n);
#endif
}

/* ------------- Output comparison -------------
 * Fixed point: the RVV kernels replicate the C arithmetic exactly (wrapping
 * adds commute), so any difference is a real bug -> bit-exact. Float: the
 * kernels use FMAs and reordered sums -> compare relative to the buffer
 * peak. The int16 output of mdf_deemph_output is bit-exact in both builds
 * (the float WORD2INT kernel replicates floor(.5+x) exactly). */
#define MDF_SMUL_F32_REL_TOL  1e-5
#define MDF_WSMUL_F32_REL_TOL 1e-5
#define MDF_PS_F32_REL_TOL    1e-6
#define MDF_IP_F32_REL_TOL    1e-5
#define MDF_PROP_F32_REL_TOL  1e-4
#define MDF_ACC_F32_TOL       1e-6  /* of sum|terms|, not of the result */

#ifdef FIXED_POINT
#define mdf_buf16_matches(ref, res, n, tol) checkasm_i16_bitexact(ref, res, n)
#define mdf_buf32_matches(ref, res, n, tol) checkasm_i32_bitexact(ref, res, n)
static inline int mdf_scalar_matches(spx_word32_t ref, spx_word32_t res,
                                     double scale, double tol)
{
    (void) scale;
    (void) tol;
    if (ref != res) {
        fprintf(stderr, "FAILED: scalar ref=%ld res=%ld (bit-exact required)\n",
                (long) ref, (long) res);
        return 0;
    }
    return 1;
}
#else
#define mdf_buf16_matches(ref, res, n, tol) checkasm_f32_within_tol(ref, res, n, tol)
#define mdf_buf32_matches(ref, res, n, tol) checkasm_f32_within_tol(ref, res, n, tol)
/* rel_tol of the result plus MDF_ACC_F32_TOL of `scale`, the magnitude it was
 * accumulated from (0 if the caller has none): a reduction whose terms cancel
 * leaves a result too small to judge reassociated rounding against. */
static inline int mdf_scalar_matches(spx_word32_t ref, spx_word32_t res,
                                     double scale, double rel_tol)
{
    double diff  = fabs((double) ref - (double) res);
    double limit = rel_tol * fabs((double) ref) + MDF_ACC_F32_TOL * scale;
    if (!(diff <= limit)) {
        fprintf(stderr, "FAILED: scalar ref=%g res=%g diff=%.2e (limit %.2e)\n",
                (double) ref, (double) res, diff, limit);
        return 0;
    }
    return 1;
}
#endif

#endif /* SPEEXDSP_TESTS_CHECKASM_MDF_WRAP_H */
