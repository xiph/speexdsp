#include <stdio.h>
#include "../internal.h"
#include "wrap.h"

/* kf_bfly2/3/4/5: per-stage unit test + micro-benchmark. Every radix-2..5
 * stage of each nfft is driven with kf_work's real arguments, forward and
 * inverse. The nfft set makes the radix-2/4 kernels hit both vectorization
 * shapes (m >= N and N > m, down to the m=1 tail stages) and includes strip
 * lengths that are not multiples of any vector length. */

#define FUNC_NAME "kf_bfly"

enum { MAX_NFFT = 3600 };

static const int configs[] = { 60, 240, 250, 256, 480, 512, 1024, 3600 };

typedef void (*bfly_fn)(kiss_fft_cfg, kiss_fft_cpx *, int, int, int, int);

void checkasm_check_fft_bfly(void)
{
#ifdef HAVE_RVV_KF_BFLY
    checkasm_declare(void, kiss_fft_cfg, kiss_fft_cpx *, int, int, int, int);

    CHECKASM_ALIGN(kiss_fft_cpx in[MAX_NFFT]);
    CHECKASM_ALIGN(kiss_fft_cpx out_ref[MAX_NFFT]);
    CHECKASM_ALIGN(kiss_fft_cpx out_new[MAX_NFFT]);
    struct fft_stage st[32];

/* The kernels run in place, so out_new degrades to Inf/NaN (or 0 for
 * fixed-point) after ~128 iterations. Periodically restore it from in
 * during the bench (tidx is checkasm's bench-iteration index). */
#define BENCH_RESTORED(buf) ((tidx & 31) ? (buf) : memcpy((buf), in, (size_t) nfft * sizeof *in))

    for (size_t c = 0; c < sizeof(configs) / sizeof(configs[0]); c++) {
        const int nfft = configs[c];
        const int nst = fft_stages(nfft, st, 32);

        for (int inverse = 0; inverse <= 1; inverse++) {
            kiss_fft_cfg cfg = fft_make_cfg(nfft, inverse);
            if (!cfg) { fprintf(stderr, FUNC_NAME ": alloc failed n%d\n", nfft); continue; }
            const char *dir = inverse ? "inv" : "fwd";

            for (int s = 0; s < nst; s++) {
                bfly_fn fn_c, fn_rvv;
                switch (st[s].p) {
                case 2: fn_c = kf_bfly2_c; fn_rvv = kf_bfly2_rvv; break;
                case 3: fn_c = kf_bfly3_c; fn_rvv = kf_bfly3_rvv; break;
                case 4: fn_c = kf_bfly4_c; fn_rvv = kf_bfly4_rvv; break;
                case 5: fn_c = kf_bfly5_c; fn_rvv = kf_bfly5_rvv; break;
                default: continue; /* kf_bfly_generic stage: no kernel */
                }
                const int fs = st[s].fstride, m = st[s].m, N = st[s].N, mm = st[s].mm;

                fft_fill_input(in, nfft);

                /* C baseline, also the reference the RVV variant pairs with. */
                if (checkasm_check_func(fn_c, FUNC_NAME "%d_%s_n%d_m%d", st[s].p, dir, nfft, m)) {
                    checkasm_bench_new(cfg, BENCH_RESTORED(out_new), fs, m, N, mm);
                }

                if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
                    if (checkasm_check_func(fn_rvv, FUNC_NAME "%d_%s_n%d_m%d", st[s].p, dir, nfft, m)) {
                        memcpy(out_ref, in, (size_t) nfft * sizeof *in);
                        checkasm_call_ref(cfg, out_ref, fs, m, N, mm);
                        memcpy(out_new, in, (size_t) nfft * sizeof *in);
                        checkasm_call_new(cfg, out_new, fs, m, N, mm);
                        if (!fft_buf_matches(out_ref, out_new, nfft, KF_BFLY_F32_REL_TOL))
                            checkasm_fail();
                        checkasm_bench_new(cfg, BENCH_RESTORED(out_new), fs, m, N, mm);
                    }
                }
            }

            fft_free_cfg(cfg);
        }
    }

    checkasm_report(FUNC_NAME);
#undef BENCH_RESTORED
#endif /* HAVE_RVV_KF_BFLY */
}
