#ifndef SPEEXDSP_TESTS_CHECKASM_FFT_WRAP_H
#define SPEEXDSP_TESTS_CHECKASM_FFT_WRAP_H

#include <string.h>
#include "config.h"
#include "kiss_fft.h"
#include "../compare.h"

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
static inline void fft_fill_input(kiss_fft_cpx *buf, int nfft)
{
#ifdef FIXED_POINT
    /* Full-range random int16: bit-exactness needs no headroom. */
    checkasm_init(buf, (size_t) nfft * sizeof *buf);
#else
    checkasm_fill_symmetric_f32((float *) buf, 2 * nfft);
#endif
}

/* ------------- Output comparison -------------
 * Fixed point: the RVV butterflies replicate the C arithmetic exactly, so
 * any difference is a real bug -> bit-exact. Float: the kernels use FMAs
 * and reordered sums -> compare relative to the buffer peak, tighter for a
 * single stage than for the log(nfft)-deep transform. One comparator name
 * for both modes so the tests stay #ifdef-free. */
#define KF_BFLY_F32_REL_TOL 1e-5
#define KF_FFT_F32_REL_TOL  1e-4

static inline int fft_buf_matches(const kiss_fft_cpx *ref, const kiss_fft_cpx *res,
        int nfft, double rel_tol)
{
#ifdef FIXED_POINT
    (void) rel_tol;
    return checkasm_i16_bitexact((const spx_int16_t *) ref, (const spx_int16_t *) res, 2 * nfft);
#else
    return checkasm_f32_within_tol((const float *) ref, (const float *) res, 2 * nfft, rel_tol);
#endif
}

#endif /* SPEEXDSP_TESTS_CHECKASM_FFT_WRAP_H */
