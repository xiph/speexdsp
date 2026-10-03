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

/* NEON build of the resampler functions. Compiled only when has_neon. Keeps the
 * native USE_NEON so resample.c pulls in resample_neon.h. Only the kernels that
 * header overrides (the HAVE_NEON_* gates in wrap.h) differ from the C
 * reference and are exposed; the double-precision ones would be byte-identical
 * to C. */
#define CKA_PREFIX  ckaneon_
#define CKA_VARIANT neon
#include "wrap_resample_rename.h"
#include "wrap.h"

#ifdef USE_NEON
#  if !defined(__aarch64__) && !defined(__ARM_NEON)
#    error "wrap_resample_neon.c built without aarch64 / ARM NEON ISA; check toolchain flags"
#  endif
#endif

#include "../wrap_impl.h"
#include "wrap_resample_shims.h"

#ifdef HAVE_NEON_DIRECT_SINGLE
CKA_RESAMPLER_BASIC_SHIM(direct_single)
#endif
#ifdef HAVE_NEON_INTERPOLATE_SINGLE
CKA_RESAMPLER_BASIC_SHIM(interpolate_single)
#endif
