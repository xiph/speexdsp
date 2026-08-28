/* Copyright (C) 2006 Jean-Marc Valin
 * Copyright (C) 2026 Tristan Matthews
 */
/*
   filterbank_rvv_asm.S -- hand-written RVV 1.0 kernel for the
   psd-from-filterbank interpolation.

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

/* RVV 1.0 kernel (VLEN >= 128) for filterbank_compute_psd16, float
 * builds only. The per-bin loop is a two-index gather + weighted sum:
 *
 *   ps[i] = mel[bank_left[i]]*filter_left[i]
 *         + mel[bank_right[i]]*filter_right[i]
 *
 * The int32 index streams are unit-stride loads scaled to byte offsets
 * for two vluxei32 gathers into the (tiny, cache-resident) mel table;
 * the weights are plain vfmul/vfmacc. LMUL=4 keeps five live register
 * groups while amortizing the strip-mine overhead.
 *
 * `.option arch, +v` keeps the file self-contained (no -march needed).
 * Leaf functions, caller-saved registers only. The wrapper guarantees
 * len >= 8 and in-bounds indices (filterbank_new builds them). */

    .option arch, +v
    .text

/* ------------------------------------------------------------------ */
/* unsigned int spx_fbank_rvv_vlenb(void) -- VLEN in bytes. */
    .balign 4
    .globl  spx_fbank_rvv_vlenb
    .hidden spx_fbank_rvv_vlenb
    .type   spx_fbank_rvv_vlenb, @function
spx_fbank_rvv_vlenb:
    vsetvli a0, zero, e8, m1, ta, ma
    ret
    .size spx_fbank_rvv_vlenb, . - spx_fbank_rvv_vlenb

/* ------------------------------------------------------------------ */
/* int spx_fbank_rvv_compliant(void)
 * Reject draft RVV 0.7.1 hardware (sets VILL on a ta,ma vtype config).
 * Mirrors spx_preproc_rvv_compliant. */
    .balign 4
    .globl  spx_fbank_rvv_compliant
    .hidden spx_fbank_rvv_compliant
    .type   spx_fbank_rvv_compliant, @function
spx_fbank_rvv_compliant:
    vsetvli t0, zero, e8, m1, ta, ma
    csrr    a0, vtype
    sgtz    a0, a0
    ret
    .size spx_fbank_rvv_compliant, . - spx_fbank_rvv_compliant

/* ------------------------------------------------------------------ */
/* void spx_fbank_rvv_psd16_f32(const int *bank_left, const int *bank_right,
 *                              const float *filter_left,
 *                              const float *filter_right,
 *                              const float *mel, float *ps, int len)
 * a0 = bank_left, a1 = bank_right, a2 = filter_left, a3 = filter_right,
 * a4 = mel, a5 = ps, a6 = len */
    .balign 4
    .globl  spx_fbank_rvv_psd16_f32
    .hidden spx_fbank_rvv_psd16_f32
    .type   spx_fbank_rvv_psd16_f32, @function
spx_fbank_rvv_psd16_f32:
1:  vsetvli t0, a6, e32, m4, ta, ma
    vle32.v    v8,  (a0)           /* bank_left[i] */
    vle32.v    v4,  (a1)           /* bank_right[i] */
    vsll.vi    v8,  v8, 2          /* int index -> byte offset */
    vluxei32.v v16, (a4), v8       /* mel[id1] */
    vle32.v    v12, (a2)           /* filter_left[i] */
    vsll.vi    v4,  v4, 2
    vfmul.vv   v16, v16, v12
    vluxei32.v v24, (a4), v4       /* mel[id2] */
    vle32.v    v28, (a3)           /* filter_right[i] */
    slli    t1, t0, 2
    vfmacc.vv  v16, v24, v28
    add     a0, a0, t1
    add     a1, a1, t1
    vse32.v    v16, (a5)
    add     a2, a2, t1
    add     a3, a3, t1
    add     a5, a5, t1
    sub     a6, a6, t0
    bnez    a6, 1b
    ret
    .size spx_fbank_rvv_psd16_f32, . - spx_fbank_rvv_psd16_f32
