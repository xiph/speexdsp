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

#include <stdio.h>
#include "internal.h"
#include "config.h"

static const CheckasmTest tests[] = {
    /* unit micro-benchmarks of the SIMD kernels */
    { "resample_kernels", checkasm_check_resample_kernels },

    /* full-pipeline integration benchmark across sample-rate conversions */
    { "resample_process", checkasm_check_resample_process },

    /* per-stage micro-benchmarks of the KISS FFT butterflies */
    { "fft_bfly", checkasm_check_fft_bfly },

    /* whole-transform integration benchmark */
    { "fft_transform", checkasm_check_fft_transform },

    /* unit micro-benchmarks of the MDF echo canceller spectral kernels */
    { "mdf_kernels", checkasm_check_mdf_kernels },

    /* unit micro-benchmarks of the preprocessor (denoiser) per-bin kernels */
    { "preproc_kernels", checkasm_check_preproc_kernels },

    /* per-stage micro-benchmarks of the smallft radix stages */
    { "smallft_stage", checkasm_check_smallft_stage },

    /* whole-transform integration benchmark */
    { "smallft_transform", checkasm_check_smallft_transform },

    /* fftwrap's forward-FFT 1/n pre-scale loop */
    { "fftwrap_scale", checkasm_check_fftwrap_scale },

    /* filterbank's psd-from-bands gather + weighted sum */
    { "fbank_psd16", checkasm_check_fbank_psd16 },

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
#if defined(USE_RVV)
    { "RVV",  "rvv",  SPEEXDSP_CPU_FLAG_RVV },
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
