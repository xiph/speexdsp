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

/* Rename preprocess.c's public API (and the fftwrap/filterbank/mdf
 * symbols it references, which each TU stubs out) to a per-TU-unique
 * prefix so each wrap_preproc_*.c can #include the whole file without
 * link-time collisions -- same scheme as ../mdf/wrap_mdf_rename.h.
 * #define CKA_PREFIX and include this BEFORE wrap.h; it also names the
 * module source for ../wrap_impl.h. */
#include "../wrap_common.h"
#define CKA_MODULE_SRC "preprocess.c"

#define speex_preprocess_state_init      CKA_CAT(CKA_PREFIX, speex_preprocess_state_init)
#define speex_preprocess_state_destroy   CKA_CAT(CKA_PREFIX, speex_preprocess_state_destroy)
#define speex_preprocess                 CKA_CAT(CKA_PREFIX, speex_preprocess)
#define speex_preprocess_run             CKA_CAT(CKA_PREFIX, speex_preprocess_run)
#define speex_preprocess_estimate_update CKA_CAT(CKA_PREFIX, speex_preprocess_estimate_update)
#define speex_preprocess_ctl             CKA_CAT(CKA_PREFIX, speex_preprocess_ctl)

/* fftwrap/filterbank/mdf symbols preprocess.c links against; stubbed in
 * wrap_preproc_stubs.h. */
#define spx_fft_init                     CKA_CAT(CKA_PREFIX, spx_fft_init)
#define spx_fft_destroy                  CKA_CAT(CKA_PREFIX, spx_fft_destroy)
#define spx_fft                          CKA_CAT(CKA_PREFIX, spx_fft)
#define spx_ifft                         CKA_CAT(CKA_PREFIX, spx_ifft)
#define filterbank_new                   CKA_CAT(CKA_PREFIX, filterbank_new)
#define filterbank_destroy               CKA_CAT(CKA_PREFIX, filterbank_destroy)
#define filterbank_compute_bank32        CKA_CAT(CKA_PREFIX, filterbank_compute_bank32)
#define filterbank_compute_psd16         CKA_CAT(CKA_PREFIX, filterbank_compute_psd16)
#define filterbank_compute_bank          CKA_CAT(CKA_PREFIX, filterbank_compute_bank)
#define filterbank_compute_psd           CKA_CAT(CKA_PREFIX, filterbank_compute_psd)
#define speex_echo_get_residual          CKA_CAT(CKA_PREFIX, speex_echo_get_residual)
