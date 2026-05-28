#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../internal.h"
#include "wrap.h"

void test_inner_product_double(void)
{
#ifdef FIXED_POINT
    /* No double-precision routine in fixed-point mode. */
    return;
#else
#if !defined(USE_SSE2) && !defined(HAVE_NEON_INNER_PRODUCT_DOUBLE)
    /* No SIMD impl for this routine in this build. */
    return;
#else
    typedef double out_t;
    typedef float in_t;
    checkasm_declare(out_t, const in_t *, const in_t *, unsigned int);

    enum { BUF_SIZE = 1024 };
    CHECKASM_ALIGN(in_t a[BUF_SIZE]);
    CHECKASM_ALIGN(in_t b[BUF_SIZE]);
    resample_fill_float(a, BUF_SIZE);
    resample_fill_float(b, BUF_SIZE);

    const unsigned len = 256;

    if (checkasm_check_func(inner_product_double_c,
                            "inner_product_double_len%u", len))
        checkasm_bench_new(a, b, len);

#ifdef USE_SSE2
    if (active_flags & SPEEXDSP_CPU_FLAG_SSE2) {
        if (checkasm_check_func(inner_product_double_sse2,
                                "inner_product_double_len%u", len)) {
            out_t ref = checkasm_call_ref(a, b, len);
            out_t res = checkasm_call_new(a, b, len);
            if (!is_resample_result_within_tolerance_double(ref, res, len))
                checkasm_fail();
            checkasm_bench_new(a, b, len);
        }
    }
#endif

#ifdef HAVE_NEON_INNER_PRODUCT_DOUBLE
    if (active_flags & SPEEXDSP_CPU_FLAG_NEON) {
        if (checkasm_check_func(inner_product_double_neon,
                                "inner_product_double_len%u", len)) {
            out_t ref = checkasm_call_ref(a, b, len);
            out_t res = checkasm_call_new(a, b, len);
            if (!is_resample_result_within_tolerance_double(ref, res, len))
                checkasm_fail();
            checkasm_bench_new(a, b, len);
        }
    }
#endif

    checkasm_report("inner_product_double");
#endif
#endif
}
