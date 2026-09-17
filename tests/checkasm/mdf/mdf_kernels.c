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
                double scale = 0.0;
                for (int i = 0; i < len; i++)
                    scale += fabs((double) x[i] * (double) y[i]);
                if (!mdf_scalar_matches(ref, res, scale, MDF_IP_F32_REL_TOL))
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

/* speex_echo_get_residual's two elementwise loops. One multiply per
 * element in both builds, so C and RVV must match bit for bit (float
 * included: a single vfmul rounds identically to the scalar multiply). */
static void test_mdf_res_window(void)
{
    checkasm_declare(void, spx_word16_t *, const spx_word16_t *,
                     const spx_word16_t *, int);

    static const int sizes[] = { 129, 256, 320, 512, 960 };
    CHECKASM_ALIGN(spx_word16_t window[960]);
    CHECKASM_ALIGN(spx_word16_t last_y[960]);
    CHECKASM_ALIGN(spx_word16_t y_ref[960]);
    CHECKASM_ALIGN(spx_word16_t y_new[960]);

    for (size_t c = 0; c < sizeof(sizes) / sizeof(sizes[0]); c++) {
        const int N = sizes[c];

        mdf_fill_w16(window, N);
        mdf_fill_w16(last_y, N);
        /* different garbage in the outputs so a skipped lane is caught */
        mdf_fill_w16(y_ref, N);
        mdf_fill_w16(y_new, N);

        if (checkasm_check_func(mdf_res_window_c, "mdf_res_window_%d", N))
            checkasm_bench_new(y_new, window, last_y, N);

        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(mdf_res_window_rvv, "mdf_res_window_%d", N)) {
                checkasm_call_ref(y_ref, window, last_y, N);
                checkasm_call_new(y_new, window, last_y, N);
                if (!mdf_buf16_matches(y_ref, y_new, N, 0.0))
                    checkasm_fail();
                checkasm_bench_new(y_new, window, last_y, N);
            }
        }
    }

    checkasm_report("mdf_res_window");
}

static void test_mdf_res_scale(void)
{
    checkasm_declare(void, spx_word32_t *, const spx_word16_t *, int);

    static const int lens[] = { 65, 129, 241, 257 };
    CHECKASM_ALIGN(spx_word32_t res_base[257]);
    CHECKASM_ALIGN(spx_word32_t res_ref[257]);
    CHECKASM_ALIGN(spx_word32_t res_new[257]);
    spx_word16_t leak2;

    for (size_t c = 0; c < sizeof(lens) / sizeof(lens[0]); c++) {
        const int len = lens[c];

        mdf_fill_w32(res_base, len);
#ifdef FIXED_POINT
        leak2 = checkasm_rand() & 0x7fff;
#else
        /* residual_echo holds power-spectrum magnitudes, and the C loop
         * truncates through an int32 cast -- scale the fill so the
         * truncation is exercised on non-zero values */
        for (int i = 0; i < len; i++)
            res_base[i] *= 1e6f;
        leak2 = 0.7f;
#endif

        if (checkasm_check_func(mdf_res_scale_c, "mdf_res_scale_%d", len)) {
            memcpy(res_new, res_base, len * sizeof *res_new);
            checkasm_bench_new(res_new, &leak2, len);
        }

        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(mdf_res_scale_rvv, "mdf_res_scale_%d", len)) {
                memcpy(res_ref, res_base, len * sizeof *res_ref);
                memcpy(res_new, res_base, len * sizeof *res_new);
                checkasm_call_ref(res_ref, &leak2, len);
                checkasm_call_new(res_new, &leak2, len);
                if (!mdf_buf32_matches(res_ref, res_new, len, 0.0))
                    checkasm_fail();
                checkasm_bench_new(res_new, &leak2, len);
            }
        }
    }

    checkasm_report("mdf_res_scale");
}

/* The weight update W[block] += PHI: one add per element in both builds
 * (wrapping int32 in fixed, in-order vfadd in float), so C and RVV must
 * match bit for bit. */
static void test_mdf_weight_update(void)
{
    checkasm_declare(void, spx_word32_t *, const spx_word32_t *, int);

    static const int lens[] = { 4, 65, 129, 256, 512 };
    CHECKASM_ALIGN(spx_word32_t w_base[512]);
    CHECKASM_ALIGN(spx_word32_t w_ref[512]);
    CHECKASM_ALIGN(spx_word32_t w_new[512]);
    CHECKASM_ALIGN(spx_word32_t phi[512]);

    for (size_t c = 0; c < sizeof(lens) / sizeof(lens[0]); c++) {
        const int N = lens[c];

        mdf_fill_w32(w_base, N);
        mdf_fill_w32(phi, N);

        if (checkasm_check_func(mdf_weight_update_c, "mdf_weight_update_%d", N)) {
            memcpy(w_new, w_base, N * sizeof *w_new);
            checkasm_bench_new(w_new, phi, N);
        }

        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(mdf_weight_update_rvv, "mdf_weight_update_%d", N)) {
                memcpy(w_ref, w_base, N * sizeof *w_ref);
                memcpy(w_new, w_base, N * sizeof *w_new);
                checkasm_call_ref(w_ref, phi, N);
                checkasm_call_new(w_new, phi, N);
                if (!mdf_buf32_matches(w_ref, w_new, N, 0.0))
                    checkasm_fail();
                checkasm_bench_new(w_new, phi, N);
            }
        }
    }

    checkasm_report("mdf_weight_update");
}

#ifndef FIXED_POINT
/* The de-emphasized int16 output: the recurrence is scalar either way,
 * but the RVV path batch-converts with the frm=RDN WORD2INT kernel,
 * which must replicate floor(.5+x) bit for bit -- clamps, interior and
 * ties included. The "ties" configs feed exact multiples of 0.5 through
 * a zeroed recurrence (preemph = 0, e = 0) so every floor tie is hit. */
static void test_mdf_deemph_output(void)
{
    checkasm_declare(spx_word16_t, spx_int16_t *, const spx_word16_t *,
                     const spx_word16_t *, const spx_word16_t *,
                     const spx_word16_t *, int, int);

    static const struct { int len, stride, ties; } configs[] = {
        { 8, 1, 0 }, { 128, 1, 0 }, { 250, 1, 0 }, { 512, 1, 0 },
        { 128, 2, 0 }, { 256, 1, 1 }, { 512, 1, 1 },
    };
    CHECKASM_ALIGN(spx_word16_t input[512]);
    CHECKASM_ALIGN(spx_word16_t e[512]);
    CHECKASM_ALIGN(spx_int16_t out_ref[1024]);
    CHECKASM_ALIGN(spx_int16_t out_new[1024]);
    spx_word16_t preemph, mem0;

    for (size_t c = 0; c < sizeof(configs) / sizeof(configs[0]); c++) {
        const int len = configs[c].len, stride = configs[c].stride;
        int i;

        if (configs[c].ties) {
            /* exact half-integers pass through untouched */
            for (i = 0; i < len; i++) {
                input[i] = (float)((int)(checkasm_rand() % 131072) - 65536) * 0.5f;
                e[i] = 0.0f;
            }
            preemph = 0.0f;
            mem0 = 0.0f;
        } else {
            /* +-8000 inputs; the preemph=.9 recurrence amplifies mem to
             * ~+-50000 peaks, exercising clamps and interior together */
            mdf_fill_w16(input, len);
            mdf_fill_w16(e, len);
            for (i = 0; i < len; i++) {
                input[i] *= 8000.0f;
                e[i] *= 8000.0f;
            }
            preemph = 0.9f;
            mem0 = 123.25f;
        }
        /* same garbage in both outputs: untouched gap lanes (stride > 1)
         * must stay identical, a clobbered or skipped lane is caught */
        checkasm_init(out_ref, sizeof out_ref);
        memcpy(out_new, out_ref, sizeof out_new);

        if (checkasm_check_func(mdf_deemph_output_c, "mdf_deemph_output_%d_s%d%s",
                                len, stride, configs[c].ties ? "_ties" : ""))
            checkasm_bench_new(out_new, input, e, &preemph, &mem0, len, stride);

        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(mdf_deemph_output_rvv, "mdf_deemph_output_%d_s%d%s",
                                    len, stride, configs[c].ties ? "_ties" : "")) {
                spx_word16_t mem_ref = checkasm_call_ref(out_ref, input, e, &preemph, &mem0, len, stride);
                spx_word16_t mem_new = checkasm_call_new(out_new, input, e, &preemph, &mem0, len, stride);
                if (!checkasm_i16_bitexact(out_ref, out_new, len * stride) ||
                    !mdf_scalar_matches(mem_ref, mem_new, 0.0, 0.0))
                    checkasm_fail();
                checkasm_bench_new(out_new, input, e, &preemph, &mem0, len, stride);
            }
        }
    }

    checkasm_report("mdf_deemph_output");
}
#endif /* !FIXED_POINT */

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
    test_mdf_res_window();
    test_mdf_res_scale();
    test_mdf_weight_update();
#ifndef FIXED_POINT
    test_mdf_deemph_output();
#endif
#endif /* HAVE_RVV_MDF */
}
