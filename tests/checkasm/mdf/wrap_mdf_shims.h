/* External-linkage shims around mdf.c's static inline kernels, one per
 * kernel under test. Included by wrap_mdf_c.c and wrap_mdf_rvv.c after the
 * module, with CKA_VARIANT naming the TU (foo -> foo_c / foo_rvv); wrap.h
 * carries the matching prototypes. Scalars (wsmul's p, leak2, preemph,
 * mem) are passed by pointer so the checked-call trampoline never has to
 * marshal float register args. */

void CKA_SHIM(mdf_smul_accum)(const spx_word16_t *X, const spx_word32_t *Y,
                              spx_word16_t *acc, int N, int M)
{
    spectral_mul_accum(X, Y, acc, N, M);
}

void CKA_SHIM(mdf_smul_accum16)(const spx_word16_t *X, const spx_word16_t *Y,
                                spx_word16_t *acc, int N, int M)
{
    spectral_mul_accum16(X, Y, acc, N, M);
}

void CKA_SHIM(mdf_wsmul_conj)(const spx_float_t *w, const spx_float_t *p,
                              const spx_word16_t *X, const spx_word16_t *Y,
                              spx_word32_t *prod, int N)
{
    weighted_spectral_mul_conj(w, *p, X, Y, prod, N);
}

void CKA_SHIM(mdf_power_spectrum)(const spx_word16_t *X, spx_word32_t *ps, int N)
{
    power_spectrum(X, ps, N);
}

void CKA_SHIM(mdf_power_spectrum_accum)(const spx_word16_t *X, spx_word32_t *ps, int N)
{
    power_spectrum_accum(X, ps, N);
}

spx_word32_t CKA_SHIM(mdf_inner_prod)(const spx_word16_t *x, const spx_word16_t *y, int len)
{
    return mdf_inner_prod(x, y, len);
}

void CKA_SHIM(mdf_adjust_prop)(const spx_word32_t *W, int N, int M, int P,
                               spx_word16_t *prop)
{
    mdf_adjust_prop(W, N, M, P, prop);
}

void CKA_SHIM(mdf_res_window)(spx_word16_t *y, const spx_word16_t *window,
                              const spx_word16_t *last_y, int N)
{
    mdf_residual_window(y, window, last_y, N);
}

void CKA_SHIM(mdf_res_scale)(spx_word32_t *residual_echo, const spx_word16_t *leak2, int len)
{
    mdf_residual_scale(residual_echo, *leak2, len);
}

void CKA_SHIM(mdf_weight_update)(spx_word32_t *w, const spx_word32_t *phi, int N)
{
    mdf_weight_update(w, phi, N);
}

spx_word16_t CKA_SHIM(mdf_deemph_output)(spx_int16_t *out, const spx_word16_t *input,
                                         const spx_word16_t *e, const spx_word16_t *preemph,
                                         const spx_word16_t *mem, int len, int stride)
{
    return mdf_deemph_output(out, input, e, *preemph, *mem, len, stride);
}
