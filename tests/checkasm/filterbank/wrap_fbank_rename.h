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

/* Rename filterbank.c's public API to a per-TU-unique prefix so each
 * wrap_fbank_*.c can #include the whole file without link-time
 * collisions -- same scheme as ../preprocess/wrap_preproc_rename.h.
 * #define CKA_PREFIX and include this BEFORE wrap.h; it also names the
 * module source for ../wrap_impl.h. */
#include "../wrap_common.h"
#define CKA_MODULE_SRC "filterbank.c"

#define filterbank_new            CKA_CAT(CKA_PREFIX, filterbank_new)
#define filterbank_destroy        CKA_CAT(CKA_PREFIX, filterbank_destroy)
#define filterbank_compute_bank32 CKA_CAT(CKA_PREFIX, filterbank_compute_bank32)
#define filterbank_compute_psd16  CKA_CAT(CKA_PREFIX, filterbank_compute_psd16)
#define filterbank_compute_bank   CKA_CAT(CKA_PREFIX, filterbank_compute_bank)
#define filterbank_compute_psd    CKA_CAT(CKA_PREFIX, filterbank_compute_psd)
#define filterbank_psy_smooth     CKA_CAT(CKA_PREFIX, filterbank_psy_smooth)
