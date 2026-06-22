#include <stdio.h>
#include "../internal.h"
#include "wrap.h"

/* resampler_basic_direct_single: direct sinc table, single-precision (spx_word32_t)
 * accumulator. Inner kernel is inner_product_single -> SIMD via SSE (float),
 * NEON (both modes) or RVV (both modes; bit-exact in fixed point).
 *
 * Parameterized over several (rate, quality) configs that all select the
 * direct-single kernel but produce different filter lengths, so the SIMD inner
 * products are exercised at lengths both divisible and not divisible by the
 * vector width. In particular the aarch64 NEON kernel strides 16 elements with a
 * `len % 16` remainder loop, and filt_len is always a multiple of 8, so configs
 * with filt_len % 16 == 8 (e.g. the 3:2 downsamples below -> 72, 120) hit the
 * remainder path while the others (8, 48, 80, 160) do not. Each config asserts
 * it landed in the expected kind and is skipped (with a note) otherwise. */

#define FUNC_NAME "resampler_inner_product_single"

enum { OUT_LEN = 256, IN_LEN = 2048, IN_BUF = 4096 };

static const struct { unsigned in_rate, out_rate; int quality; } configs[] = {
    { 24000, 48000, 5 },   /* 1:2  up,   filt_len 80  (% 16 == 0) */
    { 48000, 32000, 5 },   /* 3:2  down, filt_len 120 (% 16 == 8, remainder) */
    {  8000, 16000, 3 },   /* 1:2  up,   filt_len 48  (% 16 == 0) */
    { 48000, 32000, 3 },   /* 3:2  down, filt_len 72  (% 16 == 8, remainder) */
    { 44100, 22050, 5 },   /* 2:1  down, filt_len 160 (% 16 == 0) */
    { 24000, 48000, 0 },   /* 1:2  up,   filt_len 8   (% 16 == 8, all remainder) */
};

void test_resampler_basic_direct_single(void)
{
#if (defined(USE_SSE) && !defined(FIXED_POINT)) || defined(HAVE_NEON_DIRECT_SINGLE) || defined(HAVE_RVV_DIRECT_SINGLE)
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
        if (resample_kind(st) != RESAMPLE_KIND_DIRECT_SINGLE) {
            fprintf(stderr, "%s: %u->%u q%d unexpected kind %d (filt_len %u), skipping\n",
                    FUNC_NAME, ir, orr, q, resample_kind(st), resample_filt_len(st));
            resample_destroy_state(st);
            continue;
        }

        /* C baseline (also the reference the SIMD variants are paired against). */
        if (checkasm_check_func(resampler_basic_direct_single_c, FUNC_NAME "_%u_%u_q%d", ir, orr, q)) {
            inl = IN_LEN; outl = OUT_LEN;
            checkasm_bench_new(st, 0, in, &inl, out_new, &outl);
        }

#if defined(USE_SSE) && !defined(FIXED_POINT)
        if (active_flags & SPEEXDSP_CPU_FLAG_SSE) {
            if (checkasm_check_func(resampler_basic_direct_single_sse, FUNC_NAME "_%u_%u_q%d", ir, orr, q)) {
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
#endif

#ifdef HAVE_NEON_DIRECT_SINGLE
        if (active_flags & SPEEXDSP_CPU_FLAG_NEON) {
            if (checkasm_check_func(resampler_basic_direct_single_neon, FUNC_NAME "_%u_%u_q%d", ir, orr, q)) {
                inl = IN_LEN; outl = OUT_LEN;
                int rc = checkasm_call_ref(st, 0, in, &inl, out_ref, &outl);
                inl = IN_LEN; outl = OUT_LEN;
                int rn = checkasm_call_new(st, 0, in, &inl, out_new, &outl);
#ifdef FIXED_POINT
                if (rc != rn || !resample_buffer_within_tol(out_ref, out_new, (unsigned) rc, 0.0))
#else
                if (rc != rn || !resample_buffer_within_tol(out_ref, out_new, (unsigned) rc, RESAMPLE_FLOAT_REL_TOL))
#endif
                    checkasm_fail();
                inl = IN_LEN; outl = OUT_LEN;
                checkasm_bench_new(st, 0, in, &inl, out_new, &outl);
            }
        }
#endif

#ifdef HAVE_RVV_DIRECT_SINGLE
        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(resampler_basic_direct_single_rvv, FUNC_NAME "_%u_%u_q%d", ir, orr, q)) {
                inl = IN_LEN; outl = OUT_LEN;
                int rc = checkasm_call_ref(st, 0, in, &inl, out_ref, &outl);
                inl = IN_LEN; outl = OUT_LEN;
                int rn = checkasm_call_new(st, 0, in, &inl, out_new, &outl);
                /* fixed-point RVV is bit-exact vs C (see wrap.h) */
#ifdef FIXED_POINT
                if (rc != rn || !resample_buffer_bitexact(out_ref, out_new, (unsigned) rc))
#else
                if (rc != rn || !resample_buffer_within_tol(out_ref, out_new, (unsigned) rc, RESAMPLE_FLOAT_REL_TOL))
#endif
                    checkasm_fail();
                inl = IN_LEN; outl = OUT_LEN;
                checkasm_bench_new(st, 0, in, &inl, out_new, &outl);
            }
        }
#endif

        resample_destroy_state(st);
    }

    checkasm_report(FUNC_NAME);
#endif /* any SIMD available */
}
