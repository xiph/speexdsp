#include <stdio.h>
#include "../internal.h"
#include "wrap.h"

/* preprocess.c per-bin kernels: unit test + micro-benchmark, C vs RVV
 * (float builds only -- the RVV preprocess kernels have no fixed-point
 * twins). The N set covers a sub-vector-length frame (5), interior-bin
 * counts both divisible and not divisible by any vector length, the
 * typical narrowband/wideband sizes (160/256) and an odd size (257).
 * In/out buffers are re-seeded from a base copy before each variant so C
 * and RVV see identical inputs. */

#if !defined(FIXED_POINT) && defined(HAVE_RVV_PREPROC)

static const int test_sizes[] = { 5, 32, 160, 250, 256, 257, 512 };
#define NUM_SIZES (sizeof(test_sizes) / sizeof(test_sizes[0]))

enum { MAX_N = 512, NBANDS = 24 };

static void test_preproc_window(void)
{
    checkasm_declare(void, float *, const float *, int);

    CHECKASM_ALIGN(float frame_base[2 * MAX_N]);
    CHECKASM_ALIGN(float frame_ref[2 * MAX_N]);
    CHECKASM_ALIGN(float frame_new[2 * MAX_N]);
    CHECKASM_ALIGN(float window[2 * MAX_N]);

    for (size_t c = 0; c < NUM_SIZES; c++) {
        const int len = 2 * test_sizes[c];

        preproc_fill_signal(frame_base, len);
        preproc_fill_gain(window, len);

        if (checkasm_check_func(preproc_window_c, "preproc_window_%d", len)) {
            memcpy(frame_new, frame_base, len * sizeof *frame_new);
            checkasm_bench_new(frame_new, window, len);
        }

        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(preproc_window_rvv, "preproc_window_%d", len)) {
                memcpy(frame_ref, frame_base, len * sizeof *frame_ref);
                memcpy(frame_new, frame_base, len * sizeof *frame_new);
                checkasm_call_ref(frame_ref, window, len);
                checkasm_call_new(frame_new, window, len);
                if (!checkasm_f32_within_tol(frame_ref, frame_new, len, PREPROC_ELTWISE_F32_REL_TOL))
                    checkasm_fail();
                checkasm_bench_new(frame_new, window, len);
            }
        }
    }

    checkasm_report("preproc_window");
}

static void test_preproc_power_spectrum(void)
{
    checkasm_declare(void, const float *, float *, int);

    CHECKASM_ALIGN(float ft[2 * MAX_N]);
    CHECKASM_ALIGN(float ps_ref[MAX_N]);
    CHECKASM_ALIGN(float ps_new[MAX_N]);

    for (size_t c = 0; c < NUM_SIZES; c++) {
        const int N = test_sizes[c];

        preproc_fill_signal(ft, 2 * N);
        /* pre-fill the outputs with different garbage so a skipped bin is
         * caught */
        preproc_fill_power(ps_ref, N);
        preproc_fill_power(ps_new, N);

        if (checkasm_check_func(preproc_power_spectrum_c, "preproc_power_spectrum_%d", N))
            checkasm_bench_new(ft, ps_new, N);

        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(preproc_power_spectrum_rvv, "preproc_power_spectrum_%d", N)) {
                checkasm_call_ref(ft, ps_ref, N);
                checkasm_call_new(ft, ps_new, N);
                if (!checkasm_f32_within_tol(ps_ref, ps_new, N, PREPROC_ELTWISE_F32_REL_TOL))
                    checkasm_fail();
                checkasm_bench_new(ft, ps_new, N);
            }
        }
    }

    checkasm_report("preproc_power_spectrum");
}

static void test_preproc_smooth_spectrum(void)
{
    checkasm_declare(void, float *, const float *, int);

    CHECKASM_ALIGN(float S_base[MAX_N]);
    CHECKASM_ALIGN(float S_ref[MAX_N]);
    CHECKASM_ALIGN(float S_new[MAX_N]);
    CHECKASM_ALIGN(float ps[MAX_N]);

    for (size_t c = 0; c < NUM_SIZES; c++) {
        const int N = test_sizes[c];

        preproc_fill_power(S_base, N);
        preproc_fill_power(ps, N);

        if (checkasm_check_func(preproc_smooth_spectrum_c, "preproc_smooth_spectrum_%d", N)) {
            memcpy(S_new, S_base, N * sizeof *S_new);
            checkasm_bench_new(S_new, ps, N);
        }

        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(preproc_smooth_spectrum_rvv, "preproc_smooth_spectrum_%d", N)) {
                memcpy(S_ref, S_base, N * sizeof *S_ref);
                memcpy(S_new, S_base, N * sizeof *S_new);
                checkasm_call_ref(S_ref, ps, N);
                checkasm_call_new(S_new, ps, N);
                if (!checkasm_f32_within_tol(S_ref, S_new, N, PREPROC_ELTWISE_F32_REL_TOL))
                    checkasm_fail();
                checkasm_bench_new(S_new, ps, N);
            }
        }
    }

    checkasm_report("preproc_smooth_spectrum");
}

static void test_preproc_min_track(void)
{
    checkasm_declare(void, float *, float *, const float *, int);

    CHECKASM_ALIGN(float Smin_base[MAX_N]);
    CHECKASM_ALIGN(float Stmp_base[MAX_N]);
    CHECKASM_ALIGN(float Smin_ref[MAX_N]);
    CHECKASM_ALIGN(float Stmp_ref[MAX_N]);
    CHECKASM_ALIGN(float Smin_new[MAX_N]);
    CHECKASM_ALIGN(float Stmp_new[MAX_N]);
    CHECKASM_ALIGN(float S[MAX_N]);

    for (int swap = 0; swap <= 1; swap++) {
        const char *name = swap ? "preproc_min_track_swap" : "preproc_min_track";

        for (size_t c = 0; c < NUM_SIZES; c++) {
            const int N = test_sizes[c];

            preproc_fill_power(Smin_base, N);
            preproc_fill_power(Stmp_base, N);
            preproc_fill_power(S, N);

            if (checkasm_check_func(swap ? preproc_min_track_swap_c : preproc_min_track_c,
                                    "%s_%d", name, N)) {
                memcpy(Smin_new, Smin_base, N * sizeof *Smin_new);
                memcpy(Stmp_new, Stmp_base, N * sizeof *Stmp_new);
                checkasm_bench_new(Smin_new, Stmp_new, S, N);
            }

            if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
                if (checkasm_check_func(swap ? preproc_min_track_swap_rvv : preproc_min_track_rvv,
                                        "%s_%d", name, N)) {
                    memcpy(Smin_ref, Smin_base, N * sizeof *Smin_ref);
                    memcpy(Stmp_ref, Stmp_base, N * sizeof *Stmp_ref);
                    memcpy(Smin_new, Smin_base, N * sizeof *Smin_new);
                    memcpy(Stmp_new, Stmp_base, N * sizeof *Stmp_new);
                    checkasm_call_ref(Smin_ref, Stmp_ref, S, N);
                    checkasm_call_new(Smin_new, Stmp_new, S, N);
                    /* pure min/copy: bit-exact, but reuse the tol helper */
                    if (!checkasm_f32_within_tol(Smin_ref, Smin_new, N, 0.0) ||
                        !checkasm_f32_within_tol(Stmp_ref, Stmp_new, N, 0.0))
                        checkasm_fail();
                    checkasm_bench_new(Smin_new, Stmp_new, S, N);
                }
            }
        }
    }

    checkasm_report("preproc_min_track");
}

static void test_preproc_update_prob(void)
{
    checkasm_declare(void, const float *, const float *, int *, int);

    CHECKASM_ALIGN(float S[MAX_N]);
    CHECKASM_ALIGN(float Smin[MAX_N]);
    CHECKASM_ALIGN(int up_ref[MAX_N]);
    CHECKASM_ALIGN(int up_new[MAX_N]);

    for (size_t c = 0; c < NUM_SIZES; c++) {
        const int N = test_sizes[c];

        preproc_fill_power(S, N);
        /* keep Smin in the same range so both branches get exercised */
        preproc_fill_power(Smin, N);
        for (int i = 0; i < N; i++)
            Smin[i] *= 0.4f;

        if (checkasm_check_func(preproc_update_prob_c, "preproc_update_prob_%d", N))
            checkasm_bench_new(S, Smin, up_new, N);

        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(preproc_update_prob_rvv, "preproc_update_prob_%d", N)) {
                checkasm_call_ref(S, Smin, up_ref, N);
                checkasm_call_new(S, Smin, up_new, N);
                /* the threshold multiply rounds identically in both, so
                 * the flags must match exactly */
                if (!checkasm_i32_bitexact((const spx_int32_t *) up_ref, (const spx_int32_t *) up_new, N))
                    checkasm_fail();
                checkasm_bench_new(S, Smin, up_new, N);
            }
        }
    }

    checkasm_report("preproc_update_prob");
}

static void test_preproc_noise_update(void)
{
    checkasm_declare(void, const int *, const float *, float *, const float *, int);

    CHECKASM_ALIGN(float noise_base[MAX_N]);
    CHECKASM_ALIGN(float noise_ref[MAX_N]);
    CHECKASM_ALIGN(float noise_new[MAX_N]);
    CHECKASM_ALIGN(float ps[MAX_N]);
    CHECKASM_ALIGN(int up[MAX_N]);
    static const float beta[2] = { 0.03f, 0.97f };

    for (size_t c = 0; c < NUM_SIZES; c++) {
        const int N = test_sizes[c];

        preproc_fill_power(noise_base, N);
        preproc_fill_power(ps, N);
        preproc_fill_prob(up, N);

        if (checkasm_check_func(preproc_noise_update_c, "preproc_noise_update_%d", N)) {
            memcpy(noise_new, noise_base, N * sizeof *noise_new);
            checkasm_bench_new(up, ps, noise_new, beta, N);
        }

        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(preproc_noise_update_rvv, "preproc_noise_update_%d", N)) {
                memcpy(noise_ref, noise_base, N * sizeof *noise_ref);
                memcpy(noise_new, noise_base, N * sizeof *noise_new);
                checkasm_call_ref(up, ps, noise_ref, beta, N);
                checkasm_call_new(up, ps, noise_new, beta, N);
                if (!checkasm_f32_within_tol(noise_ref, noise_new, N, PREPROC_SNR_F32_REL_TOL))
                    checkasm_fail();
                checkasm_bench_new(up, ps, noise_new, beta, N);
            }
        }
    }

    checkasm_report("preproc_noise_update");
}

static void test_preproc_snr_update(void)
{
    checkasm_declare(void, const float *, const float *, const float *,
                     const float *, const float *, float *, float *, int);

    CHECKASM_ALIGN(float ps[MAX_N + NBANDS]);
    CHECKASM_ALIGN(float noise[MAX_N + NBANDS]);
    CHECKASM_ALIGN(float echo[MAX_N + NBANDS]);
    CHECKASM_ALIGN(float reverb[MAX_N + NBANDS]);
    CHECKASM_ALIGN(float old_ps[MAX_N + NBANDS]);
    CHECKASM_ALIGN(float post_ref[MAX_N + NBANDS]);
    CHECKASM_ALIGN(float prior_ref[MAX_N + NBANDS]);
    CHECKASM_ALIGN(float post_new[MAX_N + NBANDS]);
    CHECKASM_ALIGN(float prior_new[MAX_N + NBANDS]);

    for (size_t c = 0; c < NUM_SIZES; c++) {
        const int len = test_sizes[c] + NBANDS; /* N+M bins in preprocess.c */

        preproc_fill_power(ps, len);
        preproc_fill_power(noise, len);
        preproc_fill_power(echo, len);
        preproc_fill_power(reverb, len);
        preproc_fill_power(old_ps, len);

        if (checkasm_check_func(preproc_snr_update_c, "preproc_snr_update_%d", len))
            checkasm_bench_new(ps, noise, echo, reverb, old_ps, post_new, prior_new, len);

        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(preproc_snr_update_rvv, "preproc_snr_update_%d", len)) {
                checkasm_call_ref(ps, noise, echo, reverb, old_ps, post_ref, prior_ref, len);
                checkasm_call_new(ps, noise, echo, reverb, old_ps, post_new, prior_new, len);
                if (!checkasm_f32_within_tol(post_ref, post_new, len, PREPROC_SNR_F32_REL_TOL) ||
                    !checkasm_f32_within_tol(prior_ref, prior_new, len, PREPROC_SNR_F32_REL_TOL))
                    checkasm_fail();
                checkasm_bench_new(ps, noise, echo, reverb, old_ps, post_new, prior_new, len);
            }
        }
    }

    checkasm_report("preproc_snr_update");
}

static void test_preproc_zeta_smooth(void)
{
    checkasm_declare(void, float *, const float *, int, int);

    CHECKASM_ALIGN(float zeta_base[MAX_N + NBANDS]);
    CHECKASM_ALIGN(float zeta_ref[MAX_N + NBANDS]);
    CHECKASM_ALIGN(float zeta_new[MAX_N + NBANDS]);
    CHECKASM_ALIGN(float prior[MAX_N + NBANDS]);

    for (size_t c = 0; c < NUM_SIZES; c++) {
        const int N = test_sizes[c];
        const int len = N + NBANDS;

        preproc_fill_snr(zeta_base, len);
        preproc_fill_snr(prior, len);

        if (checkasm_check_func(preproc_zeta_smooth_c, "preproc_zeta_smooth_%d", N)) {
            memcpy(zeta_new, zeta_base, len * sizeof *zeta_new);
            checkasm_bench_new(zeta_new, prior, N, NBANDS);
        }

        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(preproc_zeta_smooth_rvv, "preproc_zeta_smooth_%d", N)) {
                memcpy(zeta_ref, zeta_base, len * sizeof *zeta_ref);
                memcpy(zeta_new, zeta_base, len * sizeof *zeta_new);
                checkasm_call_ref(zeta_ref, prior, N, NBANDS);
                checkasm_call_new(zeta_new, prior, N, NBANDS);
                if (!checkasm_f32_within_tol(zeta_ref, zeta_new, len, PREPROC_ELTWISE_F32_REL_TOL))
                    checkasm_fail();
                checkasm_bench_new(zeta_new, prior, N, NBANDS);
            }
        }
    }

    checkasm_report("preproc_zeta_smooth");
}

static void test_preproc_em_gain(void)
{
    checkasm_declare(void, const float *, const float *, const float *,
                     const float *, float *, float *, float *, int);

    CHECKASM_ALIGN(float prior[MAX_N]);
    CHECKASM_ALIGN(float post[MAX_N]);
    CHECKASM_ALIGN(float ps[MAX_N]);
    CHECKASM_ALIGN(float gain_floor[MAX_N]);
    CHECKASM_ALIGN(float gain_base[MAX_N]);
    CHECKASM_ALIGN(float gain2_base[MAX_N]);
    CHECKASM_ALIGN(float old_ps_base[MAX_N]);
    CHECKASM_ALIGN(float gain_ref[MAX_N]);
    CHECKASM_ALIGN(float gain2_ref[MAX_N]);
    CHECKASM_ALIGN(float old_ps_ref[MAX_N]);
    CHECKASM_ALIGN(float gain_new[MAX_N]);
    CHECKASM_ALIGN(float gain2_new[MAX_N]);
    CHECKASM_ALIGN(float old_ps_new[MAX_N]);

    for (size_t c = 0; c < NUM_SIZES; c++) {
        const int N = test_sizes[c];

        preproc_fill_snr(prior, N);
        /* post ranges over [-1, 100] in preprocess.c */
        preproc_fill_snr(post, N);
        for (int i = 0; i < N; i++)
            post[i] -= 1.0f;
        preproc_fill_power(ps, N);
        preproc_fill_gain(gain_floor, N);
        preproc_fill_gain(gain_base, N);
        preproc_fill_gain(gain2_base, N);
        preproc_fill_power(old_ps_base, N);

        /* Force a few bins onto hypergeom_gain's tails: theta > 10 needs
         * prior >> post's reach, theta < 0 cannot occur (kept clamped),
         * but small theta (ind == 0) needs tiny prior. */
        if (N >= 8) {
            prior[1] = 0.0f;                 /* theta = 0, ind = 0 */
            prior[2] = 100.0f; post[2] = 100.0f; /* theta ~ 100, ind > 19 */
        }

        if (checkasm_check_func(preproc_em_gain_c, "preproc_em_gain_%d", N)) {
            memcpy(gain_new, gain_base, N * sizeof *gain_new);
            memcpy(gain2_new, gain2_base, N * sizeof *gain2_new);
            memcpy(old_ps_new, old_ps_base, N * sizeof *old_ps_new);
            checkasm_bench_new(prior, post, ps, gain_floor, gain_new, gain2_new, old_ps_new, N);
        }

        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(preproc_em_gain_rvv, "preproc_em_gain_%d", N)) {
                memcpy(gain_ref, gain_base, N * sizeof *gain_ref);
                memcpy(gain2_ref, gain2_base, N * sizeof *gain2_ref);
                memcpy(old_ps_ref, old_ps_base, N * sizeof *old_ps_ref);
                memcpy(gain_new, gain_base, N * sizeof *gain_new);
                memcpy(gain2_new, gain2_base, N * sizeof *gain2_new);
                memcpy(old_ps_new, old_ps_base, N * sizeof *old_ps_new);
                checkasm_call_ref(prior, post, ps, gain_floor, gain_ref, gain2_ref, old_ps_ref, N);
                checkasm_call_new(prior, post, ps, gain_floor, gain_new, gain2_new, old_ps_new, N);
                if (!checkasm_f32_within_tol(gain_ref, gain_new, N, PREPROC_EM_F32_REL_TOL) ||
                    !checkasm_f32_within_tol(gain2_ref, gain2_new, N, PREPROC_EM_F32_REL_TOL) ||
                    !checkasm_f32_within_tol(old_ps_ref, old_ps_new, N, PREPROC_EM_F32_REL_TOL))
                    checkasm_fail();
                checkasm_bench_new(prior, post, ps, gain_floor, gain_new, gain2_new, old_ps_new, N);
            }
        }
    }

    checkasm_report("preproc_em_gain");
}

static void test_preproc_apply_gain(void)
{
    checkasm_declare(void, const float *, float *, int);

    CHECKASM_ALIGN(float gain2[MAX_N]);
    CHECKASM_ALIGN(float ft_base[2 * MAX_N]);
    CHECKASM_ALIGN(float ft_ref[2 * MAX_N]);
    CHECKASM_ALIGN(float ft_new[2 * MAX_N]);

    for (size_t c = 0; c < NUM_SIZES; c++) {
        const int N = test_sizes[c];

        preproc_fill_gain(gain2, N);
        preproc_fill_signal(ft_base, 2 * N);

        if (checkasm_check_func(preproc_apply_gain_c, "preproc_apply_gain_%d", N)) {
            memcpy(ft_new, ft_base, 2 * N * sizeof *ft_new);
            checkasm_bench_new(gain2, ft_new, N);
        }

        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(preproc_apply_gain_rvv, "preproc_apply_gain_%d", N)) {
                memcpy(ft_ref, ft_base, 2 * N * sizeof *ft_ref);
                memcpy(ft_new, ft_base, 2 * N * sizeof *ft_new);
                checkasm_call_ref(gain2, ft_ref, N);
                checkasm_call_new(gain2, ft_new, N);
                if (!checkasm_f32_within_tol(ft_ref, ft_new, 2 * N, PREPROC_ELTWISE_F32_REL_TOL))
                    checkasm_fail();
                checkasm_bench_new(gain2, ft_new, N);
            }
        }
    }

    checkasm_report("preproc_apply_gain");
}

/* The overlap-add int16 output: the RVV kernel adds under the ambient
 * rounding mode and replicates WORD2INT's floor(.5+x) with an frm=RDN
 * convert, so it must match the C loop bit for bit -- clamps, interior
 * and ties included. The "ties" pass feeds exact multiples of 0.5 so
 * every floor tie is hit. */
static void test_preproc_overlap_output(void)
{
    checkasm_declare(void, spx_int16_t *, const float *, const float *, int);

    CHECKASM_ALIGN(float outbuf[MAX_N]);
    CHECKASM_ALIGN(float frame[MAX_N]);
    CHECKASM_ALIGN(spx_int16_t x_ref[MAX_N]);
    CHECKASM_ALIGN(spx_int16_t x_new[MAX_N]);

    for (int ties = 0; ties <= 1; ties++) {
        for (size_t c = 0; c < NUM_SIZES; c++) {
            const int len = test_sizes[c];

            if (ties) {
                /* both terms exact multiples of 0.5: the sums stay exact
                 * (well under 2^24), so every half-integer tie is hit */
                for (int i = 0; i < len; i++) {
                    outbuf[i] = (float)((int)(checkasm_rand() % 131072) - 65536) * 0.5f;
                    frame[i]  = (float)((int)(checkasm_rand() % 65536) - 32768) * 0.5f;
                }
            } else {
                /* +-40000 sums: clamps and interior together */
                preproc_fill_signal(outbuf, len);
                preproc_fill_signal(frame, len);
                for (int i = 0; i < len; i++) {
                    outbuf[i] *= 30000.0f;
                    frame[i]  *= 10000.0f;
                }
            }
            /* different garbage per output so a skipped lane is caught */
            checkasm_init(x_ref, len * sizeof *x_ref);
            checkasm_init(x_new, len * sizeof *x_new);

            if (checkasm_check_func(preproc_overlap_output_c, "preproc_overlap_output_%d%s",
                                    len, ties ? "_ties" : ""))
                checkasm_bench_new(x_new, outbuf, frame, len);

            if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
                if (checkasm_check_func(preproc_overlap_output_rvv, "preproc_overlap_output_%d%s",
                                        len, ties ? "_ties" : "")) {
                    checkasm_call_ref(x_ref, outbuf, frame, len);
                    checkasm_call_new(x_new, outbuf, frame, len);
                    if (!checkasm_i16_bitexact(x_ref, x_new, len))
                        checkasm_fail();
                    checkasm_bench_new(x_new, outbuf, frame, len);
                }
            }
        }
    }

    checkasm_report("preproc_overlap_output");
}

#endif /* !FIXED_POINT && HAVE_RVV_PREPROC */

void checkasm_check_preproc_kernels(void)
{
#if !defined(FIXED_POINT) && defined(HAVE_RVV_PREPROC)
    test_preproc_window();
    test_preproc_power_spectrum();
    test_preproc_smooth_spectrum();
    test_preproc_min_track();
    test_preproc_update_prob();
    test_preproc_noise_update();
    test_preproc_snr_update();
    test_preproc_zeta_smooth();
    test_preproc_em_gain();
    test_preproc_apply_gain();
    test_preproc_overlap_output();
#endif
}
