#ifndef SPEEXDSP_TESTS_CHECKASM_MDF_WRAP_H
#define SPEEXDSP_TESTS_CHECKASM_MDF_WRAP_H

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "config.h"
#include "arch.h"
#include "pseudofloat.h"

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
#endif

/* ------------- Test-input fill ------------- */
#include <checkasm/utils.h>

static inline void mdf_fill_w16(spx_word16_t *buf, int n)
{
#ifdef FIXED_POINT
    /* Full-range random int16: bit-exactness needs no headroom. */
    checkasm_init(buf, (size_t) n * sizeof *buf);
#else
    /* Symmetric [-1, 1) floats (raw random bits would be NaN/Inf soup). */
    checkasm_randomize_rangef(buf, n, 2.0f);
    for (int i = 0; i < n; i++)
        buf[i] -= 1.0f;
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
 * peak. */
#define MDF_SMUL_F32_REL_TOL  1e-5
#define MDF_WSMUL_F32_REL_TOL 1e-5
#define MDF_PS_F32_REL_TOL    1e-6
#define MDF_IP_F32_REL_TOL    1e-5
#define MDF_PROP_F32_REL_TOL  1e-4

#ifdef FIXED_POINT
#define mdf_buf16_matches(ref, res, n, tol) mdf_buf_bitexact_w16(ref, res, n)
#define mdf_buf32_matches(ref, res, n, tol) mdf_buf_bitexact_w32(ref, res, n)
static inline int mdf_buf_bitexact_w16(const spx_word16_t *ref, const spx_word16_t *res, int n)
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
static inline int mdf_buf_bitexact_w32(const spx_word32_t *ref, const spx_word32_t *res, int n)
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
static inline int mdf_scalar_matches(spx_word32_t ref, spx_word32_t res, double tol)
{
    (void) tol;
    if (ref != res) {
        fprintf(stderr, "FAILED: scalar ref=%ld res=%ld (bit-exact required)\n",
                (long) ref, (long) res);
        return 0;
    }
    return 1;
}
#else
static inline int mdf_buf_within_tol(const float *ref, const float *res, int n, double rel_tol)
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
#define mdf_buf16_matches(ref, res, n, tol) mdf_buf_within_tol(ref, res, n, tol)
#define mdf_buf32_matches(ref, res, n, tol) mdf_buf_within_tol(ref, res, n, tol)
static inline int mdf_scalar_matches(spx_word32_t ref, spx_word32_t res, double rel_tol)
{
    double diff = fabs((double) ref - (double) res);
    double rel  = fabs((double) ref) > 0.0 ? diff / fabs((double) ref) : diff;
    if (!(rel <= rel_tol)) {
        fprintf(stderr, "FAILED: scalar ref=%g res=%g rel=%.2e (tol %g)\n",
                (double) ref, (double) res, rel, rel_tol);
        return 0;
    }
    return 1;
}
#endif

#endif /* SPEEXDSP_TESTS_CHECKASM_MDF_WRAP_H */
