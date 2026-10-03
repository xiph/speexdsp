/* Copyright (C) 2026 Tristan Matthews */
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

#ifndef SPEEXDSP_TESTS_CHECKASM_WRAP_COMMON_H
#define SPEEXDSP_TESTS_CHECKASM_WRAP_COMMON_H

/* Shared plumbing for the wrap_<mod>_<variant>.c translation units. Each
 * one #includes a whole libspeexdsp module (see wrap_impl.h) to get a C or
 * SIMD build of the module's static kernels, so the module's public API is
 * renamed to a per-TU prefix (the module's wrap_<mod>_rename.h) and the
 * shims exposing the kernels carry a per-TU suffix (wrap_<mod>_shims.h).
 *
 *   #define CKA_PREFIX  ckamc_        unique token, before the rename header
 *   #define CKA_VARIANT c             c / sse / neon / rvv, before the shims
 */

#ifndef CKA_PREFIX
#  error "define CKA_PREFIX (a unique token) before including the wrap headers"
#endif

#define CKA_CAT2(a, b) a ## b
#define CKA_CAT(a, b)  CKA_CAT2(a, b)

/* CKA_SHIM(foo) -> foo_c, foo_rvv, ... per CKA_VARIANT. The name is pasted
 * in the outer macro so a kernel name that the module's RVV header has
 * #define'd to its dispatch wrapper is not expanded first. */
#define CKA_SHIM(name)    CKA_SHIM_1(name ## _, CKA_VARIANT)
#define CKA_SHIM_1(a, b)  CKA_CAT2(a, b)

#endif /* SPEEXDSP_TESTS_CHECKASM_WRAP_COMMON_H */
