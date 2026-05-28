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

#define BUF_SIZE 1024

void test_inner_product_single(void)
{
#if !defined(HAVE_NEON_INNER_PRODUCT_SINGLE) && !(defined(USE_SSE) && !defined(FIXED_POINT))
    /* No SIMD impl for this routine in this build. */
    return;
#else
    checkasm_declare(out_t, const in_t *, const in_t *, unsigned int);

    CHECKASM_ALIGN(in_t a[BUF_SIZE]);
    CHECKASM_ALIGN(in_t b[BUF_SIZE]);
    INITIALIZE_BUF(a);
    INITIALIZE_BUF(b);

#ifdef FIXED_POINT
    /* Bound input magnitude so the C reference's int32 accumulator cannot
     * overflow; NEON accumulates in int64 then saturates. Without this mask
     * the two paths diverge nondeterministically for any non-trivial len. */
    for (unsigned i = 0; i < BUF_SIZE; i++) {
        a[i] >>= 6;
        b[i] >>= 6;
    }
#endif

    /* Representative production filter length: Q7 upsample base / Q10 with
     * a 2x downsample. Real speexdsp filt_len values fall in 80..256. */
    const unsigned len = 256;

    /* C baseline. */
    if (checkasm_check_func(inner_product_single_c,
                            "inner_product_single_len%u", len))
        checkasm_bench_new(a, b, len);

#ifdef HAVE_NEON_INNER_PRODUCT_SINGLE
    if (active_flags & SPEEXDSP_CPU_FLAG_NEON) {
        if (checkasm_check_func(inner_product_single_neon,
                                "inner_product_single_len%u", len)) {
            out_t ref = checkasm_call_ref(a, b, len);
            out_t res = checkasm_call_new(a, b, len);
#ifdef FIXED_POINT
            /* C clamps symmetrically to +/-32767 via SATURATE32PSHR; NEON's
             * sqrshrn saturates to the full int16 range and rounds-half-away-
             * from-zero on negative ties. Tolerate a 1-LSB divergence. */
            int32_t diff = ref > res ? ref - res : res - ref;
            if (diff > 1) {
                fprintf(stderr,
                    "FAILED: len=%u ref=%" PRId32 " res=%" PRId32 " diff=%" PRId32 "\n",
                    len, ref, res, diff);
                checkasm_fail();
            }
#else
            if (fabsf((float)(ref - res)) > 1e-4f) {
                fprintf(stderr, "FAILED: len=%u ref=%f res=%f\n",
                    len, (float)ref, (float)res);
                checkasm_fail();
            }
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
            if (fabsf((float)(ref - res)) > 1e-4f) {
                fprintf(stderr, "FAILED: len=%u ref=%f res=%f\n",
                    len, (float)ref, (float)res);
                checkasm_fail();
            }
            checkasm_bench_new(a, b, len);
        }
    }
#endif

    checkasm_report("inner_product_single");
#endif /* any SIMD available */
}
