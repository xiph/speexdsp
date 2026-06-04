#include <stdio.h>
#include "internal.h"
#include "config.h"

static const CheckasmTest tests[] = {
    /* unit micro-benchmarks of the SIMD kernels */
    { "resample_kernels", checkasm_check_resample_kernels },

    /* full-pipeline integration benchmark across sample-rate conversions */
    { "resample_process", checkasm_check_resample_process },

    { 0 }
};

/* List of cpu flags to check */
static const CheckasmCpuInfo flags[] = {
#if defined(USE_SSE)
    { "SSE",  "sse",  SPEEXDSP_CPU_FLAG_SSE },
#endif
#if defined(USE_SSE2)
    { "SSE2", "sse2", SPEEXDSP_CPU_FLAG_SSE2 },
#endif
#if defined(USE_NEON)
    { "NEON", "neon", SPEEXDSP_CPU_FLAG_NEON },
#endif
    { 0 }
};

CheckasmCpu active_flags;

static void set_cpu_flags(CheckasmCpu flags)
{
    active_flags = flags;
}

int main(int argc, const char *argv[])
{
    CheckasmConfig cfg = {
        .cpu_flags      = flags,
        .tests          = tests,
        .set_cpu_flags  = set_cpu_flags,
    };

    cfg.cpu = detect_cpu_flags();

    return checkasm_main(&cfg, argc, argv);
}
