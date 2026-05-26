#ifndef SPEEXDSP_TESTS_CHECKASM_INTERNAL_H
#define SPEEXDSP_TESTS_CHECKASM_INTERNAL_H

#include "config.h"
#include <checkasm/checkasm.h>
#include <checkasm/test.h>

/* Random number generation */
#define rnd checkasm_rand

/* CPU Flags */
#define SPEEXDSP_CPU_FLAG_SSE   (1 << 0)
#define SPEEXDSP_CPU_FLAG_SSE2  (1 << 1)
#define SPEEXDSP_CPU_FLAG_NEON  (1 << 2)

/* Runtime CPU flags */
extern CheckasmCpu active_flags;

/* Declarations of test functions */
void checkasm_check_test_routine(void);

#endif /* SPEEXDSP_TESTS_CHECKASM_INTERNAL_H */
