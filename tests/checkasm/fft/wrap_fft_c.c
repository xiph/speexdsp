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

/* Pure-C reference build of the kf_bfly* butterflies plus the shared cfg and
 * stage helpers. #undef USE_RVV keeps the butterflies scalar; this TU is the
 * benchmark baseline, so it is built with the no-autovec flags
 * (checkasm_c_ref_args in tests/meson.build). */
#define CKA_PREFIX  ckafc_
#define CKA_VARIANT c
#include "wrap_fft_rename.h"
#include "wrap.h"

#undef USE_RVV

#include "../wrap_impl.h"
#include "wrap_fft_shims.h"

/* ------------- cfg + stage helpers ------------- */

kiss_fft_cfg fft_make_cfg(int nfft, int inverse)
{
    return kiss_fft_alloc(nfft, inverse, NULL, NULL);
}

void fft_free_cfg(kiss_fft_cfg cfg)
{
    speex_free(cfg);
}

/* Mirrors kf_work's traversal: level i has fstride = N = product of the
 * earlier radices and mm = the parent stage length (1 for the outermost
 * level, where N==1 leaves it unused). */
int fft_stages(int nfft, struct fft_stage *out, int max)
{
    int facbuf[2 * MAXFACTORS];
    int n = 0, fstride = 1, prev_m = 1, i;

    kf_factor(nfft, facbuf);
    for (i = 0; n < max; i++) {
        int p = facbuf[2 * i], m = facbuf[2 * i + 1];
        out[n].p       = p;
        out[n].m       = m;
        out[n].fstride = fstride;
        out[n].N       = fstride;
        out[n].mm      = (i == 0) ? 1 : prev_m;
        n++;
        if (m == 1)
            break;
        prev_m = m;
        fstride *= p;
    }
    return n;
}
