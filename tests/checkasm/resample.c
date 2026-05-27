#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "internal.h"
#include "../../libspeexdsp/resample.h"

#ifdef FIXED_POINT
typedef int32_t out_t;
typedef int16_t in_t;
extern out_t inner_product_single_c(const in_t *a, const in_t *b, unsigned int len);
extern out_t inner_product_single_neon(const in_t *a, const in_t *b, unsigned int len);
#else
typedef float out_t;
typedef float in_t;
extern out_t inner_product_single_c(const in_t *a, const in_t *b, unsigned int len);
extern out_t inner_product_single_neon(const in_t *a, const in_t *b, unsigned int len);
#endif

#define BUF_SIZE 1024

static void check_inner_product_single(void)
{
    checkasm_declare(out_t, const in_t *, const in_t *, unsigned int);

    // Input test vectors
    CHECKASM_ALIGN(in_t a[BUF_SIZE]);
    CHECKASM_ALIGN(in_t b[BUF_SIZE]);
    INITIALIZE_BUF(a);
    INITIALIZE_BUF(b);

    // Test a variety of lengths, including those satisfying NEON constraints (len % 4 == 0 and len >= 4)
    // and some that don't, to ensure robust behavior.
    const unsigned lengths[] = { 4, 8, 12, 16, 32, 64, 128, 256, 512, 1024 };

    for (int i = 0; i < sizeof(lengths) / sizeof(lengths[0]); i++) {
        unsigned len = lengths[i];

        // Benchmark C reference as baseline
        if (checkasm_check_func(inner_product_single_c, "inner_product_single")) {
            checkasm_bench_new(a, b, len);
        }

        // Test and Benchmark NEON
        if (active_flags & SPEEXDSP_CPU_FLAG_NEON) {
            if (checkasm_check_func(inner_product_single_neon, "inner_product_single")) {
                // NEON implementation constraint: len % 4 == 0 and len >= 4.
                out_t ref = checkasm_call_ref(a, b, len);
                out_t res = checkasm_call_new(a, b, len);
#ifdef FIXED_POINT
                if (ref != res) {
                    fprintf(stderr, "FAILED: a=%hd, b=%hd, len=%u, ref=%d, res=%d\n", *a, *b, len, ref, res);
                    checkasm_fail();
                }
#else
                if (fabsf((float)(ref - res)) > 1e-4f) {
                    fprintf(stderr, "FAILED: len=%u, ref=%f, res=%f\n", len, (float)ref, (float)res);
                    checkasm_fail();
                }
#endif
                checkasm_bench_new(a, b, len);
            }
        }
    }

    checkasm_report("inner_product_single");
}

void checkasm_check_resampler(void)
{
    check_inner_product_single();
}
