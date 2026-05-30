#include <stdio.h>
#include <stdlib.h>
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

void test_inner_product_single(void)
{
#if !defined(HAVE_NEON_INNER_PRODUCT_SINGLE) && !(defined(USE_SSE) && !defined(FIXED_POINT))
    /* No SIMD impl for this routine in this build. */
    return;
#else
    checkasm_declare(out_t, const in_t *, const in_t *, unsigned int);

    enum { BUF_SIZE = 1024 };
    CHECKASM_ALIGN(in_t a[BUF_SIZE]);
    CHECKASM_ALIGN(in_t b[BUF_SIZE]);

#ifdef FIXED_POINT
    resample_fill_word16(a, BUF_SIZE);
    resample_fill_word16(b, BUF_SIZE);
#else
    resample_fill_float(a, BUF_SIZE);
    resample_fill_float(b, BUF_SIZE);
#endif

    /* Filter Length: 8 (quality 0, fastest) to 256 (quality 10, highest
     * quality). Sweep two lengths so the asm's main-loop and remainder-loop
     * paths get exercised differently across runs. */
    static const unsigned lens[] = { 128, 256 };
    for (unsigned li = 0; li < sizeof(lens)/sizeof(lens[0]); li++) {
        const unsigned len = lens[li];

        /* C baseline */
        if (checkasm_check_func(inner_product_single_c,
                                "inner_product_single_len%u", len))
            checkasm_bench_new(a, b, len);

#ifdef HAVE_NEON_INNER_PRODUCT_SINGLE
        if (active_flags & SPEEXDSP_CPU_FLAG_NEON) {
            /* Name must match the C baseline so checkasm_call_ref dispatches
             * to the C function (the harness keys ref/new pairing by name). */
            if (checkasm_check_func(inner_product_single_neon,
                                    "inner_product_single_len%u", len)) {
                out_t ref = checkasm_call_ref(a, b, len);
                out_t res = checkasm_call_new(a, b, len);
#ifdef FIXED_POINT
                if (!is_resample_result_within_tolerance_word16(ref, res, len))
                    checkasm_fail();
#else
                double scale = resample_abs_dot_scale(a, b, len);
                if (!is_resample_result_within_tolerance_float(ref, res, scale, len))
                    checkasm_fail();
#endif
                checkasm_bench_new(a, b, len);
            }
        }
#endif

#if defined(USE_SSE) && !defined(FIXED_POINT)
        if (active_flags & SPEEXDSP_CPU_FLAG_SSE) {
            if (checkasm_check_func(inner_product_single_sse,
                                    "inner_product_single_len%u", len)) {
                out_t ref = checkasm_call_ref(a, b, len);
                out_t res = checkasm_call_new(a, b, len);
                double scale = resample_abs_dot_scale(a, b, len);
                if (!is_resample_result_within_tolerance_float(ref, res, scale, len))
                    checkasm_fail();
                checkasm_bench_new(a, b, len);
            }
        }
#endif
    }

    checkasm_report("inner_product_single");
#endif /* any SIMD available */
}
