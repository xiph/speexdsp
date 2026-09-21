/* Copyright (C) 2026 Tristan Matthews */
/**
   @file rvv_cpu.c
   @brief Runtime probe shared by the RISC-V Vector (RVV) kernels
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "rvv_cpu.h"

#if defined(__linux__)
#include <sys/auxv.h>

/* The two probes below are the only vector instructions outside the
 * *_rvv_asm.S files. `.option push / arch, +v / pop` scopes the extension
 * to each statement, the same mechanism the build probes for, so this file
 * still compiles with a base-ISA -march. Call them only once the 'V' HWCAP
 * bit says the vector unit is present and enabled. */

/* VLEN in bytes (vlmax for e8m1). */
static unsigned int rvv_vlenb(void)
{
   unsigned long vl;
   __asm__ volatile (".option push\n\t"
                     ".option arch, +v\n\t"
                     "vsetvli %0, zero, e8, m1, ta, ma\n\t"
                     ".option pop"
                     : "=r" (vl));
   return (unsigned int)vl;
}

/* Reject draft RVV 0.7.1 hardware, which also reports 'V' through HWCAP.
 * The tail/mask-agnostic (ta, ma) vtype config was only introduced in RVV
 * 0.9, so 0.7.1 hardware flags it illegal by setting the VILL bit -- vtype's
 * sign bit (XLEN-1) -- while compliant 1.0 hardware leaves it clear, giving
 * a positive vtype. A signed > 0 test thus returns 1 on 1.0 and 0 on 0.7.1.
 * Mirrors dav1d's has_compliant_rvv. */
static int rvv_compliant(void)
{
   long vtype;
   __asm__ volatile (".option push\n\t"
                     ".option arch, +v\n\t"
                     "vsetvli t0, zero, e8, m1, ta, ma\n\t"
                     "csrr %0, vtype\n\t"
                     ".option pop"
                     : "=r" (vtype) : : "t0");
   return vtype > 0;
}
#endif /* __linux__ */

int spx_rvv_detect(void)
{
   static int detected = -1;
   if (detected < 0)
   {
      detected = 0;
#if defined(__linux__)
      /* 'V' HWCAP bit, then reject draft RVV 0.7.1 hardware (which also
         sets it) via the vtype/VILL probe, and require VLEN >= 128 (the
         kernels' precondition; V mandates Zvl128b, but verify it directly). */
      if (getauxval(AT_HWCAP) & (1UL << ('V' - 'A')))
         detected = rvv_compliant() && rvv_vlenb() >= 16;
#endif
   }
   return detected;
}
