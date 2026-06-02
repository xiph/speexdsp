#include <stdio.h>
#include "../internal.h"
#include "wrap.h"

/* resampler_basic_direct_double: direct sinc table, double-precision accumulator.
 * Inner kernel is inner_product_double -> SIMD via SSE2 (floating point only).
 * Driven with a 24000->48000 (den_rate=2), quality-10 state (quality>8 selects
 * the double-precision variant). Both the C reference and the SSE2 path form
 * single-precision products, so they agree to ~double eps. */

#define FUNC_NAME "resampler_inner_product_double"

enum { OUT_LEN = 256, IN_LEN = 1024, IN_BUF = 1024 };

void test_resampler_basic_direct_double(void)
{
#if defined(USE_SSE2) && !defined(FIXED_POINT)
    SpeexResamplerState *st = resample_make_state(24000, 48000, 10);
    if (!st) { fprintf(stderr, "resampler_basic_direct_double: init failed\n"); return; }
    if (resample_kind(st) != RESAMPLE_KIND_DIRECT_DOUBLE) {
        fprintf(stderr, "resampler_basic_direct_double: unexpected kind %d\n", resample_kind(st));
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

    if (checkasm_check_func(resampler_basic_direct_double_c, FUNC_NAME)) {
        inl = IN_LEN; outl = OUT_LEN;
        checkasm_bench_new(st, 0, in, &inl, out_new, &outl);
    }

    if (active_flags & SPEEXDSP_CPU_FLAG_SSE2) {
        if (checkasm_check_func(resampler_basic_direct_double_sse2, FUNC_NAME)) {
            inl = IN_LEN; outl = OUT_LEN;
            int rc = checkasm_call_ref(st, 0, in, &inl, out_ref, &outl);
            inl = IN_LEN; outl = OUT_LEN;
            int rn = checkasm_call_new(st, 0, in, &inl, out_new, &outl);
            if (rc != rn || !resample_buffer_within_tol(out_ref, out_new, (unsigned) rc, RESAMPLE_DOUBLE_REL_TOL))
                checkasm_fail();
            inl = IN_LEN; outl = OUT_LEN;
            checkasm_bench_new(st, 0, in, &inl, out_new, &outl);
        }
    }

    resample_destroy_state(st);
    checkasm_report(FUNC_NAME);
#endif /* SSE2 && !FIXED_POINT */
}
