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

/* RVV build of the resampler kernels for checkasm. Base-ISA C; the V
 * instructions live only in resample_rvv_asm.S (linked in alongside). Keeps the
 * native USE_RVV so resample.c pulls in resample_rvv.h.
 *
 * RESAMPLE_RVV_FORCE_ON pins SPX_RESAMPLE_ARCH to SPX_ARCH_RVV so the asm is
 * tested without reading the state's arch. checkasm still runs this variant
 * only when spx_rvv_detect() passes. */
#define CKA_PREFIX  ckarvv_
#define CKA_VARIANT rvv
#include "wrap_resample_rename.h"
#include "wrap.h"

#ifdef USE_RVV
#  define RESAMPLE_RVV_FORCE_ON
#endif

#include "../wrap_impl.h"

#ifdef USE_RVV
#include "wrap_resample_shims.h"

CKA_RESAMPLER_BASIC_SHIM(direct_single)
CKA_RESAMPLER_BASIC_SHIM(interpolate_single)
#ifndef FIXED_POINT
CKA_RESAMPLER_BASIC_SHIM(direct_double)
CKA_RESAMPLER_BASIC_SHIM(interpolate_double)
#endif
#endif /* USE_RVV */
