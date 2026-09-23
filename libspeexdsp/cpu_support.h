/* Copyright (C) 2026 Tristan Matthews */
/**
   @file cpu_support.h
   @brief Runtime selection of the kernel set a state will use
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

#ifndef CPU_SUPPORT_H
#define CPU_SUPPORT_H

/* Per-state kernel selection: each state stores the level spx_select_arch()
 * picks at init, so there is no global to race on.
 * Levels are cumulative within a target; dispatch tests arch >= SPX_ARCH_X.
 * Only RISC-V dispatches at runtime; SSE/NEON stay compile-time. */

#define SPX_ARCH_C 0

#if defined(USE_RVV)

#include "rvv_cpu.h"

#define SPX_ARCH_RVV 1

static inline int spx_select_arch(void)
{
   return spx_rvv_detect() ? SPX_ARCH_RVV : SPX_ARCH_C;
}

#else

static inline int spx_select_arch(void)
{
   return SPX_ARCH_C;
}

#endif

#endif /* CPU_SUPPORT_H */
