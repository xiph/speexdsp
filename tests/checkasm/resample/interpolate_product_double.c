#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../internal.h"
#include "wrap.h"

void test_interpolate_product_double(void)
{
#ifdef FIXED_POINT
    /* No double-precision routine in fixed-point mode. */
    return;
#else
#if !defined(USE_SSE2) && !defined(HAVE_NEON_INTERPOLATE_PRODUCT_DOUBLE)
    /* No SIMD impl for this routine in this build. */
    return;
#else
    typedef double out_t;
    typedef float in_t;
    checkasm_declare(out_t, const in_t *, const in_t *, unsigned int,
                     spx_uint32_t, const in_t *);

    enum {
        LEN        = 256,
        OVERSAMPLE = 8,
        BUF_A_SIZE = LEN,
        BUF_B_SIZE = LEN * OVERSAMPLE + 4,
    };
    CHECKASM_ALIGN(in_t a[BUF_A_SIZE]);
    CHECKASM_ALIGN(in_t b[BUF_B_SIZE]);
    resample_fill_float(a, BUF_A_SIZE);
    resample_fill_float(b, BUF_B_SIZE);

    const in_t interp[4] = { -0.0625f, 0.5625f, 0.5625f, -0.0625f };

    if (checkasm_check_func(interpolate_product_double_c,
                            "interpolate_product_double_len%u", LEN))
        checkasm_bench_new(a, b, LEN, (spx_uint32_t)OVERSAMPLE, interp);

#ifdef USE_SSE2
    if (active_flags & SPEEXDSP_CPU_FLAG_SSE2) {
        if (checkasm_check_func(interpolate_product_double_sse2,
                                "interpolate_product_double_len%u", LEN)) {
            out_t ref = checkasm_call_ref(a, b, LEN, (spx_uint32_t)OVERSAMPLE, interp);
            out_t res = checkasm_call_new(a, b, LEN, (spx_uint32_t)OVERSAMPLE, interp);
            double scale = resample_abs_interp_scale(a, b, LEN, OVERSAMPLE, interp);
            if (!is_resample_result_within_tolerance_double(ref, res, scale, LEN))
                checkasm_fail();
            checkasm_bench_new(a, b, LEN, (spx_uint32_t)OVERSAMPLE, interp);
        }
    }
#endif

#ifdef HAVE_NEON_INTERPOLATE_PRODUCT_DOUBLE
    if (active_flags & SPEEXDSP_CPU_FLAG_NEON) {
        if (checkasm_check_func(interpolate_product_double_neon,
                                "interpolate_product_double_len%u", LEN)) {
            out_t ref = checkasm_call_ref(a, b, LEN, (spx_uint32_t)OVERSAMPLE, interp);
            out_t res = checkasm_call_new(a, b, LEN, (spx_uint32_t)OVERSAMPLE, interp);
            double scale = resample_abs_interp_scale(a, b, LEN, OVERSAMPLE, interp);
            if (!is_resample_result_within_tolerance_double(ref, res, scale, LEN))
                checkasm_fail();
            checkasm_bench_new(a, b, LEN, (spx_uint32_t)OVERSAMPLE, interp);
        }
    }
#endif

    checkasm_report("interpolate_product_double");
#endif
#endif
}
