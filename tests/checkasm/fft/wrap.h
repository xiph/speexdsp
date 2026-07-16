#ifndef SPEEXDSP_TESTS_CHECKASM_FFT_WRAP_H
#define SPEEXDSP_TESTS_CHECKASM_FFT_WRAP_H

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include "config.h"
#include "kiss_fft.h"

/* Tests for the kf_bfly2/3/4/5 butterflies (per stage and whole transform),
 * C vs SIMD. Like the resampler tests, each variant #includes kiss_fft.c
 * with its public API renamed (wrap_fft_rename.h): wrap_fft_c.c forces the
 * scalar butterflies, wrap_fft_rvv.c pins the RVV dispatch on. */

/* ------------- Per-ISA availability gates ------------- */
#ifdef USE_RVV
#  define HAVE_RVV_KF_BFLY 1
#endif

/* ------------- Stage enumeration -------------
 * One entry per kf_work butterfly stage, with the arguments kf_work passes.
 * Owned by wrap_fft_c.c (kf_factor is in scope there). */
struct fft_stage {
    int p, m, fstride, N, mm;
};
int fft_stages(int nfft, struct fft_stage *out, int max);

/* ------------- State helpers (owned by wrap_fft_c.c) -------------
 * The cfg layout and contents are identical in every TU, so one cfg is
 * shared between the C and SIMD calls. */
kiss_fft_cfg fft_make_cfg(int nfft, int inverse);
void fft_free_cfg(kiss_fft_cfg cfg);

/* ------------- Functions under test -------------
 * One shared signature; the radix-3/5 wrappers replicate kf_work's
 * per-sub-FFT loop (their C butterflies take no N/mm). */
void kf_bfly2_c(kiss_fft_cfg cfg, kiss_fft_cpx *Fout, int fstride, int m, int N, int mm);
void kf_bfly3_c(kiss_fft_cfg cfg, kiss_fft_cpx *Fout, int fstride, int m, int N, int mm);
void kf_bfly4_c(kiss_fft_cfg cfg, kiss_fft_cpx *Fout, int fstride, int m, int N, int mm);
void kf_bfly5_c(kiss_fft_cfg cfg, kiss_fft_cpx *Fout, int fstride, int m, int N, int mm);
void fft_c(kiss_fft_cfg cfg, const kiss_fft_cpx *fin, kiss_fft_cpx *fout);

#ifdef HAVE_RVV_KF_BFLY
void kf_bfly2_rvv(kiss_fft_cfg cfg, kiss_fft_cpx *Fout, int fstride, int m, int N, int mm);
void kf_bfly3_rvv(kiss_fft_cfg cfg, kiss_fft_cpx *Fout, int fstride, int m, int N, int mm);
void kf_bfly4_rvv(kiss_fft_cfg cfg, kiss_fft_cpx *Fout, int fstride, int m, int N, int mm);
void kf_bfly5_rvv(kiss_fft_cfg cfg, kiss_fft_cpx *Fout, int fstride, int m, int N, int mm);
void fft_rvv(kiss_fft_cfg cfg, const kiss_fft_cpx *fin, kiss_fft_cpx *fout);
#endif

/* ------------- Test-input fill ------------- */
#include <checkasm/utils.h>

static inline void fft_fill_input(kiss_fft_cpx *buf, int nfft)
{
#ifdef FIXED_POINT
    /* Full-range random int16: bit-exactness needs no headroom. */
    checkasm_init(buf, (size_t) nfft * sizeof *buf);
#else
    /* Symmetric [-1, 1) floats (raw random bits would be NaN/Inf soup). */
    float *f = (float *) buf;
    checkasm_randomize_rangef(f, 2 * nfft, 2.0f);
    for (int i = 0; i < 2 * nfft; i++)
        f[i] -= 1.0f;
#endif
}

/* ------------- Output comparison -------------
 * Fixed point: the RVV butterflies replicate the C arithmetic exactly, so
 * any difference is a real bug -> bit-exact. Float: the kernels use FMAs
 * and reordered sums -> compare relative to the buffer peak, tighter for a
 * single stage than for the log(nfft)-deep transform. */
#define KF_BFLY_F32_REL_TOL 1e-5
#define KF_FFT_F32_REL_TOL  1e-4

#ifdef FIXED_POINT
static inline int fft_buf_bitexact(const kiss_fft_cpx *ref, const kiss_fft_cpx *res,
        int nfft)
{
    const spx_int16_t *r = (const spx_int16_t *) ref, *s = (const spx_int16_t *) res;
    for (int i = 0; i < 2 * nfft; i++) {
        if (r[i] != s[i]) {
            fprintf(stderr, "FAILED: %s[%d] ref=%d res=%d (bit-exact required)\n",
                    (i & 1) ? "im" : "re", i / 2, (int) r[i], (int) s[i]);
            return 0;
        }
    }
    return 1;
}
#else
static inline int fft_buf_within_tol(const kiss_fft_cpx *ref, const kiss_fft_cpx *res,
        int nfft, double rel_tol)
{
    const float *r = (const float *) ref, *s = (const float *) res;
    double peak = 0.0;
    for (int i = 0; i < 2 * nfft; i++) {
        double v = fabs((double) r[i]);
        if (v > peak) peak = v;
    }
    for (int i = 0; i < 2 * nfft; i++) {
        double diff = fabs((double) r[i] - (double) s[i]);
        double rel  = peak > 0.0 ? diff / peak : diff;
        /* NaN-safe: rel > rel_tol would be false for NaN and wrongly pass. */
        if (!(rel <= rel_tol)) {
            fprintf(stderr, "FAILED: %s[%d] ref=%g res=%g diff=%g peak=%g "
                    "rel=%.2e (tol %g)\n",
                    (i & 1) ? "im" : "re", i / 2, (double) r[i], (double) s[i],
                    diff, peak, rel, rel_tol);
            return 0;
        }
    }
    return 1;
}
#endif

/* One comparator name for both modes so the tests stay #ifdef-free. */
static inline int fft_buf_matches(const kiss_fft_cpx *ref, const kiss_fft_cpx *res,
        int nfft, double rel_tol)
{
#ifdef FIXED_POINT
    (void) rel_tol;
    return fft_buf_bitexact(ref, res, nfft);
#else
    return fft_buf_within_tol(ref, res, nfft, rel_tol);
#endif
}

#endif /* SPEEXDSP_TESTS_CHECKASM_FFT_WRAP_H */
