#include "internal.h"
#include "arch.h"

// Declarations of functions to benchmark
extern float test_routine_c(const float *a, const float *b, int len);
extern float test_routine_sse(const float *a, const float *b, int len);
extern float test_routine_neon(const float *a, const float *b, int len);

void checkasm_check_test_routine(void)
{
    // Always register the reference implementation
    checkasm_check_func(test_routine_c, "c");

    // Only register SIMD implementations if the runtime check passes
    if (active_flags & SPEEXDSP_CPU_FLAG_SSE)
        checkasm_check_func(test_routine_sse, "sse");

    if (active_flags & SPEEXDSP_CPU_FLAG_NEON)
        checkasm_check_func(test_routine_neon, "neon");

    // Trigger the benchmark for this group
    checkasm_report("test_routine");
}
