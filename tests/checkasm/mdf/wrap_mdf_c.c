/* Pure-C reference build of the mdf.c spectral kernels. #undef USE_RVV
 * keeps the loops scalar; this TU is the benchmark baseline, so it is
 * built with the no-autovec flags (checkasm_c_ref_args in
 * tests/meson.build). */
#define CKA_PREFIX ckamc_
#include "wrap_mdf_rename.h"
#include "wrap.h"

#undef USE_RVV

#include "wrap_mdf_impl.h"

void mdf_smul_accum_c(const spx_word16_t *X, const spx_word32_t *Y,
                      spx_word16_t *acc, int N, int M)
{
    spectral_mul_accum(X, Y, acc, N, M);
}

void mdf_smul_accum16_c(const spx_word16_t *X, const spx_word16_t *Y,
                        spx_word16_t *acc, int N, int M)
{
    spectral_mul_accum16(X, Y, acc, N, M);
}

void mdf_wsmul_conj_c(const spx_float_t *w, const spx_float_t *p,
                      const spx_word16_t *X, const spx_word16_t *Y,
                      spx_word32_t *prod, int N)
{
    weighted_spectral_mul_conj(w, *p, X, Y, prod, N);
}

void mdf_power_spectrum_c(const spx_word16_t *X, spx_word32_t *ps, int N)
{
    power_spectrum(X, ps, N);
}

void mdf_power_spectrum_accum_c(const spx_word16_t *X, spx_word32_t *ps, int N)
{
    power_spectrum_accum(X, ps, N);
}

spx_word32_t mdf_inner_prod_c(const spx_word16_t *x, const spx_word16_t *y, int len)
{
    return mdf_inner_prod(x, y, len);
}

void mdf_adjust_prop_c(const spx_word32_t *W, int N, int M, int P,
                       spx_word16_t *prop)
{
    mdf_adjust_prop(W, N, M, P, prop);
}
