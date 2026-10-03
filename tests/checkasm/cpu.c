/* Copyright (C) 2026 Tyler Kim
 * Copyright (C) 2026 Tristan Matthews
 */
/*
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

   - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   - Neither the name of the Xiph.org Foundation nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE FOUNDATION OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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
