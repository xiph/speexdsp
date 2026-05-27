#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include "internal.h"
#include "../../libspeexdsp/resample.h"

#ifdef FIXED_POINT
typedef int32_t out_t;
typedef int16_t in_t;
#else
typedef float out_t;
typedef float in_t;
#endif

extern out_t inner_product_single_c(const in_t *a, const in_t *b, unsigned int len);
#ifdef USE_NEON
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

#ifdef FIXED_POINT
    // Bound input magnitude so the C reference's int32 accumulator cannot
    // overflow. NEON accumulates in int64 then saturates, so without this
    // mask the two paths diverge nondeterministically for any non-trivial
    // length. With values in [-2^9, 2^9), worst case sum is len * 2^18,
    // well inside int32 for len up to BUF_SIZE.
    for (unsigned i = 0; i < BUF_SIZE; i++) {
        a[i] >>= 6;
        b[i] >>= 6;
    }
#endif

    // Representative filter length for production use: Q7 upsample base, or
    // Q10 with a 2x downsample. Real speexdsp filt_len values typically fall
    // in 80..256; 256 sits at the high end of common quality settings.
    const unsigned len = 256;

    // Benchmark C reference as baseline
    if (checkasm_check_func(inner_product_single_c, "inner_product_single_len%u", len)) {
        checkasm_bench_new(a, b, len);
    }

#ifdef USE_NEON
    // Test and Benchmark NEON
    if (active_flags & SPEEXDSP_CPU_FLAG_NEON) {
        if (checkasm_check_func(inner_product_single_neon, "inner_product_single_len%u", len)) {
            // NEON implementation constraint: len % 4 == 0 and len >= 4.
            out_t ref = checkasm_call_ref(a, b, len);
            out_t res = checkasm_call_new(a, b, len);
#ifdef FIXED_POINT
            // C clamps symmetrically to +/-32767 via SATURATE32PSHR;
            // NEON's sqrshrn saturates to the full int16 range and
            // rounds-half-away-from-zero on negative ties. Tolerate
            // a 1-LSB divergence on those corner cases.
            int32_t diff = ref > res ? ref - res : res - ref;
            if (diff > 1) {
                fprintf(stderr,
                    "FAILED: len=%u ref=%" PRId32 " res=%" PRId32 " diff=%" PRId32 "\n",
                    len, ref, res, diff);
                checkasm_fail();
            }
            //fprintf(stdout, "PASSED: len=%u ref=%" PRId32 " res=%" PRId32 " diff=%" PRId32 "\n", len, ref, res, diff);
#else
            if (fabsf((float)(ref - res)) > 1e-4f) {
                fprintf(stderr, "FAILED: len=%u, ref=%f, res=%f\n", len, (float)ref, (float)res);
                checkasm_fail();
            }
#endif
            checkasm_bench_new(a, b, len);
        }
    }
#endif

    checkasm_report("inner_product_single");
}

void checkasm_check_resampler(void)
{
    check_inner_product_single();
}
