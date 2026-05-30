#include <stdio.h>
#include "../internal.h"
#include "wrap.h"

/* resampler_basic_interpolate_single: interpolated sinc table, single-precision
 * accumulator. Inner kernel is interpolate_product_single -> SIMD via SSE (float
 * only; NEON has no override yet). Driven with a 44100->48000 (den_rate=160),
 * quality-5 state so update_filter selects the interpolated, single variant. */

enum { OUT_LEN = 256, IN_LEN = 1024, IN_BUF = 1024 };

void test_resampler_basic_interpolate_single(void)
{
#if defined(USE_SSE) && !defined(FIXED_POINT)
    SpeexResamplerState *st = resample_make_state(44100, 48000, 5);
    if (!st) { fprintf(stderr, "resampler_basic_interpolate_single: init failed\n"); return; }
    if (resample_kind(st) != RESAMPLE_KIND_INTERPOLATE_SINGLE) {
        fprintf(stderr, "resampler_basic_interpolate_single: unexpected kind %d\n", resample_kind(st));
        resample_destroy_state(st);
        return;
    }

    checkasm_declare(int, SpeexResamplerState *, spx_uint32_t, const spx_word16_t *,
                     spx_uint32_t *, spx_word16_t *, spx_uint32_t *);

    CHECKASM_ALIGN(spx_word16_t in[IN_BUF]);
    CHECKASM_ALIGN(spx_word16_t out_ref[OUT_LEN]);
    CHECKASM_ALIGN(spx_word16_t out_new[OUT_LEN]);
    resample_fill_input(in, IN_BUF);

    spx_uint32_t inl, outl;

    if (checkasm_check_func(resampler_basic_interpolate_single_c, "resampler_basic_interpolate_single")) {
        inl = IN_LEN; outl = OUT_LEN;
        checkasm_bench_new(st, 0, in, &inl, out_new, &outl);
    }

    if (active_flags & SPEEXDSP_CPU_FLAG_SSE) {
        if (checkasm_check_func(resampler_basic_interpolate_single_sse, "resampler_basic_interpolate_single")) {
            inl = IN_LEN; outl = OUT_LEN;
            int rc = checkasm_call_ref(st, 0, in, &inl, out_ref, &outl);
            inl = IN_LEN; outl = OUT_LEN;
            int rn = checkasm_call_new(st, 0, in, &inl, out_new, &outl);
            if (rc != rn || !resample_buffer_within_tol(out_ref, out_new, (unsigned) rc, RESAMPLE_FLOAT_REL_TOL))
                checkasm_fail();
            inl = IN_LEN; outl = OUT_LEN;
            checkasm_bench_new(st, 0, in, &inl, out_new, &outl);
        }
    }

    resample_destroy_state(st);
    checkasm_report("resampler_basic_interpolate_single");
#endif /* SSE && !FIXED_POINT */
}
