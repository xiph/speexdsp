#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include "../internal.h"
#include "wrap.h"

#ifdef FIXED_POINT
typedef int32_t out_t;
typedef int16_t in_t;
#else
typedef float out_t;
typedef float in_t;
#endif

void test_interpolate_product_single(void)
{
#if !defined(HAVE_NEON_INTERPOLATE_PRODUCT_SINGLE) \
 && !(defined(USE_SSE) && !defined(FIXED_POINT))
    /* No SIMD impl for this routine in this build. */
    return;
#else
    checkasm_declare(out_t, const in_t *, const in_t *, unsigned int,
                     spx_uint32_t, const in_t *);

    /* a: filter input window;  b: sinc table (indexed b[j*oversample + 0..3]).
     * Buffers are sized for the longest length swept below; shorter lengths
     * read a strict subset. For len=256, oversample=8 the sinc table needs
     * >= 2052 elements. */
    enum {
        MAX_LEN    = 256,
        OVERSAMPLE = 8,
        BUF_A_SIZE = MAX_LEN,
        BUF_B_SIZE = MAX_LEN * OVERSAMPLE + 4,
    };
    CHECKASM_ALIGN(in_t a[BUF_A_SIZE]);
    CHECKASM_ALIGN(in_t b[BUF_B_SIZE]);
#ifdef FIXED_POINT
    resample_fill_word16(a, BUF_A_SIZE);
    resample_fill_word16(b, BUF_B_SIZE);
#else
    resample_fill_float(a, BUF_A_SIZE);
    resample_fill_float(b, BUF_B_SIZE);
#endif

    /* Cubic coefficients matching cubic_coef(0.5, interp) — see resample.c. */
#ifdef FIXED_POINT
    const in_t interp[4] = { -2048, 18432, 18432, -2048 };  /* Q15 of {-1/16, 9/16, 9/16, -1/16} */
#else
    const in_t interp[4] = { -0.0625f, 0.5625f, 0.5625f, -0.0625f };
#endif

    /* Sweep two filter lengths so the asm gets exercised at multiple sizes.
     * Both are even (the SSE loop steps i+=2). */
    static const unsigned lens[] = { 128, 256 };
    for (unsigned li = 0; li < sizeof(lens)/sizeof(lens[0]); li++) {
        const unsigned len = lens[li];

        if (checkasm_check_func(interpolate_product_single_c,
                                "interpolate_product_single_len%u", len))
            checkasm_bench_new(a, b, len, (spx_uint32_t)OVERSAMPLE, interp);

#if defined(USE_SSE) && !defined(FIXED_POINT)
        if (active_flags & SPEEXDSP_CPU_FLAG_SSE) {
            if (checkasm_check_func(interpolate_product_single_sse,
                                    "interpolate_product_single_len%u", len)) {
                out_t ref = checkasm_call_ref(a, b, len, (spx_uint32_t)OVERSAMPLE, interp);
                out_t res = checkasm_call_new(a, b, len, (spx_uint32_t)OVERSAMPLE, interp);
                double scale = resample_abs_interp_scale(a, b, len, OVERSAMPLE, interp);
                if (!is_resample_result_within_tolerance_float(ref, res, scale, len))
                    checkasm_fail();
                checkasm_bench_new(a, b, len, (spx_uint32_t)OVERSAMPLE, interp);
            }
        }
#endif

#ifdef HAVE_NEON_INTERPOLATE_PRODUCT_SINGLE
        if (active_flags & SPEEXDSP_CPU_FLAG_NEON) {
            if (checkasm_check_func(interpolate_product_single_neon,
                                    "interpolate_product_single_len%u", len)) {
                out_t ref = checkasm_call_ref(a, b, len, (spx_uint32_t)OVERSAMPLE, interp);
                out_t res = checkasm_call_new(a, b, len, (spx_uint32_t)OVERSAMPLE, interp);
#ifdef FIXED_POINT
                if (!is_resample_result_within_tolerance_word16(ref, res, len))
                    checkasm_fail();
#else
                double scale = resample_abs_interp_scale(a, b, len, OVERSAMPLE, interp);
                if (!is_resample_result_within_tolerance_float(ref, res, scale, len))
                    checkasm_fail();
#endif
                checkasm_bench_new(a, b, len, (spx_uint32_t)OVERSAMPLE, interp);
            }
        }
#endif
    }

    checkasm_report("interpolate_product_single");
#endif
}
