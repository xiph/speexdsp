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

/* Pure-C reference build of the smallft radix stages plus the shared
 * lookup and stage helpers. #undef USE_RVV keeps the stages scalar; this
 * TU is the benchmark baseline, so it is built with the no-autovec flags
 * (checkasm_c_ref_args in tests/meson.build). */
#define CKA_PREFIX  ckasc_
#define CKA_VARIANT c
#include "wrap_smallft_rename.h"
#include "wrap.h"

#undef USE_RVV

#include "../wrap_impl.h"
#include "wrap_smallft_shims.h"

/* ------------- lookup + stage helpers ------------- */

struct drft_lookup *smallft_make_lookup(int n)
{
    struct drft_lookup *l = (struct drft_lookup *) speex_alloc(sizeof(*l));
    if (l)
        spx_drft_init(l, n);
    return l;
}

void smallft_free_lookup(struct drft_lookup *l)
{
    spx_drft_clear(l);
    speex_free(l);
}

/* Mirrors drftf1's (backward==0) / drftb1's (backward==1) stage walk,
 * recording the radix-2/4 stages with the arguments the driver passes.
 * dradfg/dradbg stages are walked (their iw advance matters) but not
 * recorded. */
int smallft_stages(struct drft_lookup *l, int backward,
                   struct drft_stage *out, int max)
{
    const int *ifac = l->splitcache;
    const int n = l->n, nf = ifac[1];
    int cnt = 0, k1;

    if (!backward) {
        int l2 = n, iw = n;
        for (k1 = 0; k1 < nf; k1++) {
            const int ip = ifac[nf - k1 + 1];
            const int l1 = l2 / ip;
            const int ido = n / l2;
            iw -= (ip - 1) * ido;
            if ((ip == 2 || ip == 4) && cnt < max) {
                out[cnt].ip = ip; out[cnt].ido = ido;
                out[cnt].l1 = l1; out[cnt].iw = iw;
                cnt++;
            }
            l2 = l1;
        }
    } else {
        int l1 = 1, iw = 1;
        for (k1 = 0; k1 < nf; k1++) {
            const int ip = ifac[k1 + 2];
            const int l2 = ip * l1;
            const int ido = n / l2;
            if ((ip == 2 || ip == 4) && cnt < max) {
                out[cnt].ip = ip; out[cnt].ido = ido;
                out[cnt].l1 = l1; out[cnt].iw = iw;
                cnt++;
            }
            l1 = l2;
            iw += (ip - 1) * ido;
        }
    }
    return cnt;
}

/* ------------- fftwrap.c's forward-FFT pre-scale, C baseline -------------
 * Kept in this no-autovec TU so the fftwrap_scale benchmark compares the
 * RVV kernel against genuinely scalar code. */
void fftwrap_scale_c(float *out, const float *in, const float *scale, int n)
{
    int i;
    for (i = 0; i < n; i++)
        out[i] = *scale * in[i];
}
