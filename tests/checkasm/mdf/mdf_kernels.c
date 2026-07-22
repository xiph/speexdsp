#include <stdio.h>
#include "../internal.h"
#include "wrap.h"

/* mdf.c spectral kernels: unit test + micro-benchmark, C vs RVV. The
 * kernels operate on half-complex (packed) spectra of even length N; the
 * N set exercises interior-bin counts both divisible and not divisible by
 * any vector length ((N-2)/2 = 31, 124, 127, 255), and M covers single-
 * block up to long-tail filters. Fixed point is bit-exact; float pairs
 * with small relative tolerances. */

#ifdef HAVE_RVV_MDF

static const struct { int N, M; } smul_configs[] = {
    { 64, 1 }, { 64, 8 }, { 250, 8 }, { 256, 1 }, { 256, 8 }, { 256, 16 },
    { 512, 8 },
};

enum { MAX_N = 512, MAX_M = 16, MAX_P = 2 };

static void test_mdf_smul_accum(void)
{
    checkasm_declare(void, const spx_word16_t *, const spx_word32_t *,
                     spx_word16_t *, int, int);

    CHECKASM_ALIGN(spx_word16_t X[MAX_N * MAX_M]);
    CHECKASM_ALIGN(spx_word32_t Y[MAX_N * MAX_M]);
    CHECKASM_ALIGN(spx_word16_t acc_ref[MAX_N]);
    CHECKASM_ALIGN(spx_word16_t acc_new[MAX_N]);

    for (size_t c = 0; c < sizeof(smul_configs) / sizeof(smul_configs[0]); c++) {
        const int N = smul_configs[c].N, M = smul_configs[c].M;

        mdf_fill_w16(X, N * M);
        mdf_fill_w32(Y, N * M);

        /* C baseline, also the reference the RVV variant pairs with. */
        if (checkasm_check_func(mdf_smul_accum_c, "mdf_smul_accum_N%d_M%d", N, M))
            checkasm_bench_new(X, Y, acc_new, N, M);

        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(mdf_smul_accum_rvv, "mdf_smul_accum_N%d_M%d", N, M)) {
                checkasm_call_ref(X, Y, acc_ref, N, M);
                checkasm_call_new(X, Y, acc_new, N, M);
                if (!mdf_buf16_matches(acc_ref, acc_new, N, MDF_SMUL_F32_REL_TOL))
                    checkasm_fail();
                checkasm_bench_new(X, Y, acc_new, N, M);
            }
        }
    }

    checkasm_report("mdf_smul_accum");
}

#ifdef FIXED_POINT
/* Float spectral_mul_accum16 is an alias of spectral_mul_accum; only the
 * fixed-point 16-bit-Y variant is a distinct kernel. */
static void test_mdf_smul_accum16(void)
{
    checkasm_declare(void, const spx_word16_t *, const spx_word16_t *,
                     spx_word16_t *, int, int);

    CHECKASM_ALIGN(spx_word16_t X[MAX_N * MAX_M]);
    CHECKASM_ALIGN(spx_word16_t Y[MAX_N * MAX_M]);
    CHECKASM_ALIGN(spx_word16_t acc_ref[MAX_N]);
    CHECKASM_ALIGN(spx_word16_t acc_new[MAX_N]);

    for (size_t c = 0; c < sizeof(smul_configs) / sizeof(smul_configs[0]); c++) {
        const int N = smul_configs[c].N, M = smul_configs[c].M;

        mdf_fill_w16(X, N * M);
        mdf_fill_w16(Y, N * M);

        if (checkasm_check_func(mdf_smul_accum16_c, "mdf_smul_accum16_N%d_M%d", N, M))
            checkasm_bench_new(X, Y, acc_new, N, M);

        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(mdf_smul_accum16_rvv, "mdf_smul_accum16_N%d_M%d", N, M)) {
                checkasm_call_ref(X, Y, acc_ref, N, M);
                checkasm_call_new(X, Y, acc_new, N, M);
                if (!mdf_buf16_matches(acc_ref, acc_new, N, MDF_SMUL_F32_REL_TOL))
                    checkasm_fail();
                checkasm_bench_new(X, Y, acc_new, N, M);
            }
        }
    }

    checkasm_report("mdf_smul_accum16");
}
#else
/* weighted_spectral_mul_conj's fixed-point pseudofloat arithmetic has no
 * RVV kernel; only the float variant is under test. */
static void test_mdf_wsmul_conj(void)
{
    checkasm_declare(void, const spx_float_t *, const spx_float_t *,
                     const spx_word16_t *, const spx_word16_t *,
                     spx_word32_t *, int);

    CHECKASM_ALIGN(spx_float_t w[MAX_N / 2 + 1]);
    CHECKASM_ALIGN(spx_word16_t X[MAX_N]);
    CHECKASM_ALIGN(spx_word16_t Y[MAX_N]);
    CHECKASM_ALIGN(spx_word32_t prod_ref[MAX_N]);
    CHECKASM_ALIGN(spx_word32_t prod_new[MAX_N]);
    spx_float_t p;

    for (size_t c = 0; c < sizeof(smul_configs) / sizeof(smul_configs[0]); c++) {
        const int N = smul_configs[c].N;
        if (c > 0 && N == smul_configs[c - 1].N)
            continue; /* M is irrelevant here */

        mdf_fill_w16(w, N / 2 + 1);
        mdf_fill_w16(X, N);
        mdf_fill_w16(Y, N);
        p = 0.5f;

        if (checkasm_check_func(mdf_wsmul_conj_c, "mdf_wsmul_conj_N%d", N))
            checkasm_bench_new(w, &p, X, Y, prod_new, N);

        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(mdf_wsmul_conj_rvv, "mdf_wsmul_conj_N%d", N)) {
                checkasm_call_ref(w, &p, X, Y, prod_ref, N);
                checkasm_call_new(w, &p, X, Y, prod_new, N);
                if (!mdf_buf32_matches(prod_ref, prod_new, N, MDF_WSMUL_F32_REL_TOL))
                    checkasm_fail();
                checkasm_bench_new(w, &p, X, Y, prod_new, N);
            }
        }
    }

    checkasm_report("mdf_wsmul_conj");
}
#endif /* FIXED_POINT / float */

static void test_mdf_power_spectrum(void)
{
    checkasm_declare(void, const spx_word16_t *, spx_word32_t *, int);

    CHECKASM_ALIGN(spx_word16_t X[MAX_N]);
    CHECKASM_ALIGN(spx_word32_t ps_base[MAX_N / 2 + 1]);
    CHECKASM_ALIGN(spx_word32_t ps_ref[MAX_N / 2 + 1]);
    CHECKASM_ALIGN(spx_word32_t ps_new[MAX_N / 2 + 1]);

    for (int accum = 0; accum <= 1; accum++) {
        const char *name = accum ? "mdf_power_spectrum_accum" : "mdf_power_spectrum";

        for (size_t c = 0; c < sizeof(smul_configs) / sizeof(smul_configs[0]); c++) {
            const int N = smul_configs[c].N;
            const int nps = N / 2 + 1;
            if (c > 0 && N == smul_configs[c - 1].N)
                continue;

            mdf_fill_w16(X, N);
            mdf_fill_w32(ps_base, nps); /* live accumulator contents */

            if (checkasm_check_func(accum ? mdf_power_spectrum_accum_c : mdf_power_spectrum_c,
                                    "%s_N%d", name, N)) {
                memcpy(ps_new, ps_base, nps * sizeof *ps_new);
                checkasm_bench_new(X, ps_new, N);
            }

            if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
                if (checkasm_check_func(accum ? mdf_power_spectrum_accum_rvv : mdf_power_spectrum_rvv,
                                        "%s_N%d", name, N)) {
                    memcpy(ps_ref, ps_base, nps * sizeof *ps_ref);
                    memcpy(ps_new, ps_base, nps * sizeof *ps_new);
                    checkasm_call_ref(X, ps_ref, N);
                    checkasm_call_new(X, ps_new, N);
                    if (!mdf_buf32_matches(ps_ref, ps_new, nps, MDF_PS_F32_REL_TOL))
                        checkasm_fail();
                    checkasm_bench_new(X, ps_new, N);
                }
            }
        }
    }

    checkasm_report("mdf_power_spectrum");
}

static void test_mdf_inner_prod(void)
{
    checkasm_declare(spx_word32_t, const spx_word16_t *, const spx_word16_t *, int);

    static const int lens[] = { 2, 8, 128, 129, 256, 1024 };
    CHECKASM_ALIGN(spx_word16_t x[1024]);
    CHECKASM_ALIGN(spx_word16_t y[1024]);

    mdf_fill_w16(x, 1024);
    mdf_fill_w16(y, 1024);

    for (size_t c = 0; c < sizeof(lens) / sizeof(lens[0]); c++) {
        const int len = lens[c];

        if (checkasm_check_func(mdf_inner_prod_c, "mdf_inner_prod_%d", len))
            checkasm_bench_new(x, y, len);

        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(mdf_inner_prod_rvv, "mdf_inner_prod_%d", len)) {
                spx_word32_t ref = checkasm_call_ref(x, y, len);
                spx_word32_t res = checkasm_call_new(x, y, len);
                if (!mdf_scalar_matches(ref, res, MDF_IP_F32_REL_TOL))
                    checkasm_fail();
                checkasm_bench_new(x, y, len);
            }
        }
    }

    checkasm_report("mdf_inner_prod");
}

static void test_mdf_adjust_prop(void)
{
    checkasm_declare(void, const spx_word32_t *, int, int, int, spx_word16_t *);

    static const struct { int N, M, P; } configs[] = {
        { 64, 4, 1 }, { 256, 8, 1 }, { 256, 8, 2 }, { 512, 8, 1 },
    };
    CHECKASM_ALIGN(spx_word32_t W[MAX_N * 8 * MAX_P]);
    CHECKASM_ALIGN(spx_word16_t prop_ref[8]);
    CHECKASM_ALIGN(spx_word16_t prop_new[8]);

    for (size_t c = 0; c < sizeof(configs) / sizeof(configs[0]); c++) {
        const int N = configs[c].N, M = configs[c].M, P = configs[c].P;

        mdf_fill_w32(W, N * M * P);

        if (checkasm_check_func(mdf_adjust_prop_c, "mdf_adjust_prop_N%d_M%d_P%d", N, M, P))
            checkasm_bench_new(W, N, M, P, prop_new);

        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(mdf_adjust_prop_rvv, "mdf_adjust_prop_N%d_M%d_P%d", N, M, P)) {
                checkasm_call_ref(W, N, M, P, prop_ref);
                checkasm_call_new(W, N, M, P, prop_new);
                if (!mdf_buf16_matches(prop_ref, prop_new, M, MDF_PROP_F32_REL_TOL))
                    checkasm_fail();
                checkasm_bench_new(W, N, M, P, prop_new);
            }
        }
    }

    checkasm_report("mdf_adjust_prop");
}

#endif /* HAVE_RVV_MDF */

void checkasm_check_mdf_kernels(void)
{
#ifdef HAVE_RVV_MDF
    test_mdf_smul_accum();
#ifdef FIXED_POINT
    test_mdf_smul_accum16();
#else
    test_mdf_wsmul_conj();
#endif
    test_mdf_power_spectrum();
    test_mdf_inner_prod();
    test_mdf_adjust_prop();
#endif /* HAVE_RVV_MDF */
}
