#include <stdio.h>
#include "../internal.h"
#include "wrap.h"

/* kiss_fft: whole-transform integration test + benchmark, C vs SIMD, both
 * directions -- the number spx_fft users (mdf, preprocess) actually see.
 * Bit-exact in fixed point; peak-relative tolerance in float. */

#define FUNC_NAME "kiss_fft"

enum { MAX_NFFT = 3600 };

static const int configs[] = { 60, 240, 250, 256, 480, 512, 1024, 3600 };

void checkasm_check_fft_transform(void)
{
#ifdef HAVE_RVV_KF_BFLY
    checkasm_declare(void, kiss_fft_cfg, const kiss_fft_cpx *, kiss_fft_cpx *);

    CHECKASM_ALIGN(kiss_fft_cpx in[MAX_NFFT]);
    CHECKASM_ALIGN(kiss_fft_cpx out_ref[MAX_NFFT]);
    CHECKASM_ALIGN(kiss_fft_cpx out_new[MAX_NFFT]);

    for (size_t c = 0; c < sizeof(configs) / sizeof(configs[0]); c++) {
        const int nfft = configs[c];

        for (int inverse = 0; inverse <= 1; inverse++) {
            kiss_fft_cfg cfg = fft_make_cfg(nfft, inverse);
            if (!cfg) { fprintf(stderr, FUNC_NAME ": alloc failed n%d\n", nfft); continue; }
            const char *dir = inverse ? "inv" : "fwd";

            fft_fill_input(in, nfft);

            if (checkasm_check_func(fft_c, FUNC_NAME "_%s_n%d", dir, nfft))
                checkasm_bench_new(cfg, in, out_new);

            if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
                if (checkasm_check_func(fft_rvv, FUNC_NAME "_%s_n%d", dir, nfft)) {
                    checkasm_call_ref(cfg, in, out_ref);
                    checkasm_call_new(cfg, in, out_new);
                    if (!fft_buf_matches(out_ref, out_new, nfft, KF_FFT_F32_REL_TOL))
                        checkasm_fail();
                    checkasm_bench_new(cfg, in, out_new);
                }
            }

            fft_free_cfg(cfg);
        }
    }

    checkasm_report(FUNC_NAME);
#endif /* HAVE_RVV_KF_BFLY */
}
