#include <stdio.h>
#include "internal.h"
#include "config.h"

static const CheckasmTest tests[] = {
    { "resample", checkasm_check_resampler },
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

    cfg.cpu = 0;
#if defined(USE_SSE)
    cfg.cpu = SPEEXDSP_CPU_FLAG_SSE;
#endif
#if defined(USE_SSE2)
    cfg.cpu = SPEEXDSP_CPU_FLAG_SSE2;
#endif
#if defined(USE_NEON)
    cfg.cpu = SPEEXDSP_CPU_FLAG_NEON;
#endif
    
    return checkasm_main(&cfg, argc, argv);
}
