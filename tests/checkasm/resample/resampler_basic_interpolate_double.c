#include <stdio.h>
#include "../internal.h"
#include "wrap.h"

/* resampler_basic_interpolate_double: interpolated sinc table, double-precision
 * accumulator. Inner kernel is interpolate_product_double -> SIMD via SSE2
 * (floating point only). The C reference promotes inputs to double before
 * multiplying while the SSE2 path multiplies in single precision, so the looser
 * float-product tolerance applies.
 *
 * Parameterized over several coprime-ratio (large den_rate) conversions at
 * quality > 8, with a range of filter lengths. Each config asserts its kind and
 * is skipped (with a note) otherwise. */

#define FUNC_NAME "resampler_interpolate_product_double"

enum { OUT_LEN = 256, IN_LEN = 2048, IN_BUF = 4096 };

void test_resampler_basic_interpolate_double(void)
{
#if defined(USE_SSE2) && !defined(FIXED_POINT)
    static const struct { unsigned in_rate, out_rate; int quality; } configs[] = {
        { 44100, 48000, 10 },   /* 147:160 up   */
        { 48000, 44100,  9 },   /* 160:147 down */
        { 44100, 32000, 10 },   /* 441:320 down */
        { 16000, 44100,  9 },   /* 160:441 up   */
    };

    checkasm_declare(int, SpeexResamplerState *, spx_uint32_t, const spx_word16_t *,
                     spx_uint32_t *, spx_word16_t *, spx_uint32_t *);

    CHECKASM_ALIGN(spx_word16_t in[IN_BUF]);
    CHECKASM_ALIGN(spx_word16_t out_ref[OUT_LEN]);
    CHECKASM_ALIGN(spx_word16_t out_new[OUT_LEN]);
    resample_fill_input(in, IN_BUF);

    spx_uint32_t inl, outl;

    for (size_t i = 0; i < sizeof(configs) / sizeof(configs[0]); i++) {
        const unsigned ir = configs[i].in_rate, orr = configs[i].out_rate;
        const int q = configs[i].quality;

        SpeexResamplerState *st = resample_make_state(ir, orr, q);
        if (!st) { fprintf(stderr, "%s: init failed %u->%u q%d\n", FUNC_NAME, ir, orr, q); continue; }
        if (resample_kind(st) != RESAMPLE_KIND_INTERPOLATE_DOUBLE) {
            fprintf(stderr, "%s: %u->%u q%d unexpected kind %d (filt_len %u), skipping\n",
                    FUNC_NAME, ir, orr, q, resample_kind(st), resample_filt_len(st));
            resample_destroy_state(st);
            continue;
        }

        if (checkasm_check_func(resampler_basic_interpolate_double_c, FUNC_NAME "_%u_%u_q%d", ir, orr, q)) {
            inl = IN_LEN; outl = OUT_LEN;
            checkasm_bench_new(st, 0, in, &inl, out_new, &outl);
        }

        if (active_flags & SPEEXDSP_CPU_FLAG_SSE2) {
            if (checkasm_check_func(resampler_basic_interpolate_double_sse2, FUNC_NAME "_%u_%u_q%d", ir, orr, q)) {
                inl = IN_LEN; outl = OUT_LEN;
                int rc = checkasm_call_ref(st, 0, in, &inl, out_ref, &outl);
                inl = IN_LEN; outl = OUT_LEN;
                int rn = checkasm_call_new(st, 0, in, &inl, out_new, &outl);
                if (rc != rn || !resample_buffer_within_tol(out_ref, out_new, (unsigned) rc, RESAMPLE_DOUBLE_FLOATMUL_REL_TOL))
                    checkasm_fail();
                inl = IN_LEN; outl = OUT_LEN;
                checkasm_bench_new(st, 0, in, &inl, out_new, &outl);
            }
        }

        resample_destroy_state(st);
    }

    checkasm_report(FUNC_NAME);
#endif /* SSE2 && !FIXED_POINT */
}
