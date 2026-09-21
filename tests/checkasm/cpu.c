#include "internal.h"

#if (defined(__arm__) || defined(__thumb2__)) && !defined(__aarch64__) && defined(__linux__)
#include <sys/auxv.h>
#include <asm/hwcap.h>
#endif

#if defined(__riscv) && defined(USE_RVV)
#include "rvv_cpu.h"   /* the library's own probe: HWCAP + RVV 0.7.1 rejection */
#endif

CheckasmCpu detect_cpu_flags(void)
{
    CheckasmCpu flags = 0;

#if defined(__aarch64__)
#if defined(USE_NEON)
    flags |= SPEEXDSP_CPU_FLAG_NEON;
#endif
#elif defined(__arm__) || defined(__thumb2__)
#if defined(USE_NEON) && defined(__linux__)
    if (getauxval(AT_HWCAP) & HWCAP_NEON)
        flags |= SPEEXDSP_CPU_FLAG_NEON;
#endif
#elif defined(__i386__) || defined(__x86_64__)
    __builtin_cpu_init();
#if defined(USE_SSE)
    if (__builtin_cpu_supports("sse"))
        flags |= SPEEXDSP_CPU_FLAG_SSE;
#endif
#if defined(USE_SSE2)
    if (__builtin_cpu_supports("sse2"))
        flags |= SPEEXDSP_CPU_FLAG_SSE2;
#endif
#elif defined(__riscv)
#if defined(USE_RVV)
    if (spx_rvv_detect())
        flags |= SPEEXDSP_CPU_FLAG_RVV;
#endif
#endif

    return flags;
}
