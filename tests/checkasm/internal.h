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
#define SPEEXDSP_CPU_FLAG_RVV   (1 << 3)

/* Runtime CPU flags */
extern CheckasmCpu active_flags;

/* Runtime CPU feature detection */
CheckasmCpu detect_cpu_flags(void);

/* Declarations of test functions */
void checkasm_check_resample_kernels(void);
void checkasm_check_resample_process(void);
void checkasm_check_fft_bfly(void);
void checkasm_check_fft_transform(void);
void checkasm_check_mdf_kernels(void);
void checkasm_check_preproc_kernels(void);
void checkasm_check_smallft_stage(void);
void checkasm_check_smallft_transform(void);
void checkasm_check_fftwrap_scale(void);
void checkasm_check_fbank_psd16(void);

#endif /* SPEEXDSP_TESTS_CHECKASM_INTERNAL_H */
