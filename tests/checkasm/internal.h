#ifndef SPEEXDSP_TESTS_CHECKASM_INTERNAL_H
#define SPEEXDSP_TESTS_CHECKASM_INTERNAL_H

#include "config.h"
#include <checkasm/checkasm.h>
#include <checkasm/test.h>
#include <checkasm/utils.h>

/* CPU Flags */
#define SPEEXDSP_CPU_FLAG_SSE   (1 << 0)
#define SPEEXDSP_CPU_FLAG_SSE2  (1 << 1)
#define SPEEXDSP_CPU_FLAG_NEON  (1 << 2)

/* Runtime CPU flags */
extern CheckasmCpu active_flags;

/* Runtime CPU feature detection */
CheckasmCpu detect_cpu_flags(void);

/* Declarations of test functions */
void checkasm_check_resampler(void);

#endif /* SPEEXDSP_TESTS_CHECKASM_INTERNAL_H */
