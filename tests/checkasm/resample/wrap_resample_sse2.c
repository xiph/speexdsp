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

/* SSE2 build of the double-precision resampler functions. Compiled only when
 * has_sse2 && !fixed-point. The native USE_SSE2 makes resample_sse.h define
 * OVERRIDE_INNER_PRODUCT_DOUBLE / OVERRIDE_INTERPOLATE_PRODUCT_DOUBLE, so the
 * double kernels use the SSE2 intrinsics. No full-pipeline wrappers (the SSE
 * TU's state already runs the SSE2 double kernels), hence
 * CKA_KERNEL_SHIMS_ONLY. */
#define CKA_PREFIX  ckasse2_
#define CKA_VARIANT sse2
#include "wrap_resample_rename.h"
#include "wrap.h"
#include "../wrap_impl.h"

#define CKA_KERNEL_SHIMS_ONLY
#include "wrap_resample_shims.h"

CKA_RESAMPLER_BASIC_SHIM(direct_double)
CKA_RESAMPLER_BASIC_SHIM(interpolate_double)
