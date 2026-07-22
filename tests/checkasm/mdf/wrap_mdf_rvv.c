/* RVV build of the mdf.c spectral kernels for checkasm: base-ISA C, the V
 * instructions live in mdf_rvv_asm.S (linked alongside). MDF_RVV_FORCE_ON
 * pins the dispatch on so the asm is tested deterministically, not per the
 * build host's getauxval. */
#define CKA_PREFIX ckamrvv_
#include "wrap_mdf_rename.h"
#include "wrap.h"

#ifdef USE_RVV
#  define MDF_RVV_FORCE_ON
#endif

#include "wrap_mdf_impl.h"

#ifdef HAVE_RVV_MDF

void mdf_smul_accum_rvv(const spx_word16_t *X, const spx_word32_t *Y,
                        spx_word16_t *acc, int N, int M)
{
    spectral_mul_accum(X, Y, acc, N, M);
}

void mdf_smul_accum16_rvv(const spx_word16_t *X, const spx_word16_t *Y,
                          spx_word16_t *acc, int N, int M)
{
    spectral_mul_accum16(X, Y, acc, N, M);
}

void mdf_wsmul_conj_rvv(const spx_float_t *w, const spx_float_t *p,
                        const spx_word16_t *X, const spx_word16_t *Y,
                        spx_word32_t *prod, int N)
{
    weighted_spectral_mul_conj(w, *p, X, Y, prod, N);
}

void mdf_power_spectrum_rvv(const spx_word16_t *X, spx_word32_t *ps, int N)
{
    power_spectrum(X, ps, N);
}

void mdf_power_spectrum_accum_rvv(const spx_word16_t *X, spx_word32_t *ps, int N)
{
    power_spectrum_accum(X, ps, N);
}

spx_word32_t mdf_inner_prod_rvv(const spx_word16_t *x, const spx_word16_t *y, int len)
{
    return mdf_inner_prod(x, y, len);
}

void mdf_adjust_prop_rvv(const spx_word32_t *W, int N, int M, int P,
                         spx_word16_t *prop)
{
    mdf_adjust_prop(W, N, M, P, prop);
}

#endif /* HAVE_RVV_MDF */
