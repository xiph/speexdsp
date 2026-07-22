#include <stdio.h>
#include "../internal.h"
#include "wrap.h"

/* smallft radix-2/4 stages: per-stage unit test + micro-benchmark, and a
 * whole-transform integration test, C vs RVV. The n set covers the ido==1
 * fast path, both part-2 vectorization axes (across pairs and across
 * rows), the l1==1 radix-2 stage of 2*4^k sizes, and odd-ido radix-4
 * stages (n = 4*odd). Stages read cc and write every element of ch, so
 * ref/new are pre-filled with identical garbage: a lane the RVV kernel
 * misses shows up as a mismatch. */

enum { MAX_N = 2048, MAX_STAGES = 16 };

static const int configs[] = { 12, 20, 64, 128, 256, 512, 1024, 2048 };

void checkasm_check_smallft_stage(void)
{
#ifdef HAVE_RVV_SMALLFT
    checkasm_declare(void, struct drft_lookup *, const struct drft_stage *,
                     int, float *, float *);

    CHECKASM_ALIGN(float in[MAX_N]);
    CHECKASM_ALIGN(float garbage[MAX_N]);
    CHECKASM_ALIGN(float out_ref[MAX_N]);
    CHECKASM_ALIGN(float out_new[MAX_N]);
    struct drft_stage st[MAX_STAGES];

    for (size_t c = 0; c < sizeof(configs) / sizeof(configs[0]); c++) {
        const int n = configs[c];
        struct drft_lookup *l = smallft_make_lookup(n);
        if (!l) { fprintf(stderr, "smallft_stage: alloc failed n%d\n", n); continue; }

        for (int backward = 0; backward <= 1; backward++) {
            const int nst = smallft_stages(l, backward, st, MAX_STAGES);
            const char *dir = backward ? "radb" : "radf";

            for (int s = 0; s < nst; s++) {
                smallft_fill_input(in, n);
                smallft_fill_input(garbage, n);

                /* C baseline, also the reference the RVV variant pairs with. */
                if (checkasm_check_func(smallft_stage_c, "drft_%s%d_n%d_ido%d_l%d",
                                        dir, st[s].ip, n, st[s].ido, st[s].l1)) {
                    checkasm_bench_new(l, &st[s], backward, in, out_new);
                }

                if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
                    if (checkasm_check_func(smallft_stage_rvv, "drft_%s%d_n%d_ido%d_l%d",
                                            dir, st[s].ip, n, st[s].ido, st[s].l1)) {
                        memcpy(out_ref, garbage, n * sizeof *out_ref);
                        memcpy(out_new, garbage, n * sizeof *out_new);
                        checkasm_call_ref(l, &st[s], backward, in, out_ref);
                        checkasm_call_new(l, &st[s], backward, in, out_new);
                        if (!smallft_buf_within_tol(out_ref, out_new, n, DRFT_STAGE_REL_TOL))
                            checkasm_fail();
                        checkasm_bench_new(l, &st[s], backward, in, out_new);
                    }
                }
            }
        }

        smallft_free_lookup(l);
    }

    checkasm_report("smallft_stage");
#endif /* HAVE_RVV_SMALLFT */
}

void checkasm_check_smallft_transform(void)
{
#ifdef HAVE_RVV_SMALLFT
    checkasm_declare(void, struct drft_lookup *, const struct drft_stage *,
                     int, float *, float *);

    CHECKASM_ALIGN(float in[MAX_N]);
    CHECKASM_ALIGN(float buf_ref[MAX_N]);
    CHECKASM_ALIGN(float buf_new[MAX_N]);

/* The transform runs in place, so repeated benching grows the data by ~n
 * per pass and would overflow. Periodically restore it from in (tidx is
 * checkasm's bench-iteration index). */
#define BENCH_RESTORED(buf) ((tidx & 31) ? (buf) : memcpy((buf), in, (size_t) n * sizeof(float)))

    for (size_t c = 0; c < sizeof(configs) / sizeof(configs[0]); c++) {
        const int n = configs[c];
        struct drft_lookup *l = smallft_make_lookup(n);
        if (!l) { fprintf(stderr, "smallft_transform: alloc failed n%d\n", n); continue; }

        for (int backward = 0; backward <= 1; backward++) {
            const char *dir = backward ? "backward" : "forward";

            smallft_fill_input(in, n);

            if (checkasm_check_func(smallft_forward_c, "drft_%s_n%d", dir, n)) {
                memcpy(buf_new, in, n * sizeof *buf_new);
                checkasm_bench_new(l, NULL, backward, BENCH_RESTORED(buf_new), NULL);
            }

            if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
                if (checkasm_check_func(smallft_forward_rvv, "drft_%s_n%d", dir, n)) {
                    memcpy(buf_ref, in, n * sizeof *buf_ref);
                    memcpy(buf_new, in, n * sizeof *buf_new);
                    checkasm_call_ref(l, NULL, backward, buf_ref, NULL);
                    checkasm_call_new(l, NULL, backward, buf_new, NULL);
                    if (!smallft_buf_within_tol(buf_ref, buf_new, n, DRFT_FFT_REL_TOL))
                        checkasm_fail();
                    memcpy(buf_new, in, n * sizeof *buf_new);
                    checkasm_bench_new(l, NULL, backward, BENCH_RESTORED(buf_new), NULL);
                }
            }
        }

        smallft_free_lookup(l);
    }

    checkasm_report("smallft_transform");
#undef BENCH_RESTORED
#endif /* HAVE_RVV_SMALLFT */
}
