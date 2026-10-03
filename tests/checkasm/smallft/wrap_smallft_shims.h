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

/* External-linkage shims around smallft.c's static radix stages and the
 * whole transform. Included by wrap_smallft_c.c and wrap_smallft_rvv.c
 * after the module, with CKA_VARIANT naming the TU; wrap.h carries the
 * matching prototypes. The stage shim derives the wa pointers from the
 * lookup exactly as drftf1/drftb1 do and calls the (SPX_DRAD*-dispatched)
 * stage. */

void CKA_SHIM(smallft_stage)(struct drft_lookup *l, const struct drft_stage *st,
                             int backward, float *cc, float *ch)
{
    float *wa = l->trigcache + l->n;
    float *wa1 = wa + st->iw - 1;
    float *wa2 = wa + st->iw + st->ido - 1;
    float *wa3 = wa + st->iw + 2 * st->ido - 1;

    if (st->ip == 4) {
        if (backward)
            SPX_DRADB4(l->arch, st->ido, st->l1, cc, ch, wa1, wa2, wa3);
        else
            SPX_DRADF4(l->arch, st->ido, st->l1, cc, ch, wa1, wa2, wa3);
    } else {
        if (backward)
            SPX_DRADB2(l->arch, st->ido, st->l1, cc, ch, wa1);
        else
            SPX_DRADF2(l->arch, st->ido, st->l1, cc, ch, wa1);
    }
}

/* Whole transform, in place; st/unused only pad the signature to the
 * stage shim's. */
void CKA_SHIM(smallft_forward)(struct drft_lookup *l, const struct drft_stage *st,
                               int backward, float *data, float *unused)
{
    (void) st;
    (void) unused;
    if (backward)
        spx_drft_backward(l, data);
    else
        spx_drft_forward(l, data);
}
