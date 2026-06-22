/* Copyright (C) 2007-2008 Jean-Marc Valin
 * Copyright (C) 2008 Thorvald Natvig
 * Copyright (C) 2026 Tyler Kim
 */
/*
   resample_rvv_asm.S -- hand-written RVV 1.0 kernels (vector loops only).

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

/* Vector-length-agnostic RVV 1.0 kernels, correct for any VLEN >= 128. Each
 * strip-mines with vsetvl: accumulators are pre-zeroed over a full vlmax group,
 * partial tail strips use the tail-undisturbed policy, and one reduction over
 * vlmax ends the kernel. Tail-agnostic accumulation would be wrong -- the spec
 * lets tail lanes hold arbitrary values. The interpolate kernels gather their
 * four sinc taps with one segment load (vlsseg4).
 *
 * `.option arch, +v` makes the file self-contained -- no -march needed: the
 * scalar mul/remu come from the base -march's M extension and flw/fcvt.d.s from
 * the d/f that +v implies. The meson probe just checks that the assembler
 * accepts `.option arch, +v`; if not, RVV is left disabled rather than failing
 * the build.
 *
 * All kernels are leaf functions using only caller-saved registers (t*, ft*,
 * a*, v*), so there is no stack frame. len (a2) is >= 8 per resample.c's
 * filt_len rounding; len == 0 also returns zero.
 *
 * The float kernels return in fa0 (the lp64d/ilp32d hard-double convention), so
 * the caller binds them only under __riscv_float_abi_double; on other float ABIs
 * the caller would read a0, so it stays scalar there. The int16 kernels return
 * integers and are float-ABI-agnostic.
 */

    .option arch, +v
    .text

/* ------------------------------------------------------------------ */
/* unsigned int spx_resample_rvv_vlenb(void)
 * VLEN in bytes (vlmax for e8m1). Call only once the 'V' HWCAP bit confirms V. */
    .balign 4
    .globl  spx_resample_rvv_vlenb
    .hidden spx_resample_rvv_vlenb
    .type   spx_resample_rvv_vlenb, @function
spx_resample_rvv_vlenb:
    vsetvli a0, zero, e8, m1, ta, ma
    ret
    .size spx_resample_rvv_vlenb, . - spx_resample_rvv_vlenb

/* ------------------------------------------------------------------ */
/* int spx_resample_rvv_compliant(void)
 * Reject draft RVV 0.7.1 hardware that also reports 'V' through HWCAP. The
 * tail/mask-agnostic (ta, ma) vtype config was only introduced in RVV 0.9, so
 * 0.7.1 hardware flags it illegal by setting the VILL bit -- vtype's sign bit
 * (XLEN-1) -- while compliant 1.0 hardware leaves it clear, giving a positive
 * vtype. A signed > 0 test thus returns 1 on 1.0 and 0 on 0.7.1. Mirrors
 * dav1d's has_compliant_rvv. Call only once the 'V' HWCAP bit confirms V. */
    .balign 4
    .globl  spx_resample_rvv_compliant
    .hidden spx_resample_rvv_compliant
    .type   spx_resample_rvv_compliant, @function
spx_resample_rvv_compliant:
    vsetvli t0, zero, e8, m1, ta, ma
    csrr    a0, vtype
    sgtz    a0, a0
    ret
    .size spx_resample_rvv_compliant, . - spx_resample_rvv_compliant

/* ------------------------------------------------------------------ */
/* float spx_resample_rvv_ip_f32(const float *a, const float *b,
 *                               unsigned int len)
 * Plain f32 dot product. Single-strip small-LMUL paths for short filters,
 * else one m4 accumulator.                                             */
    .balign 4
    .globl  spx_resample_rvv_ip_f32
    .hidden spx_resample_rvv_ip_f32
    .type   spx_resample_rvv_ip_f32, @function
spx_resample_rvv_ip_f32:
    vsetvli t0, zero, e32, m1, ta, ma
    bgeu    t0, a2, .Lip_f32_one_m1     /* len <= vlmax(e32m1): one m1 strip */
    slli    t1, t0, 1
    bgeu    t1, a2, .Lip_f32_one_m2     /* len <= vlmax(e32m2): one m2 strip */

    vsetvli t2, zero, e32, m4, ta, ma   /* vl = vlmax, t2 = vlmax(e32m4) */
    vmv.v.i v8, 0
    vmv.s.x v0, zero
    remu    t5, a2, t2
    sub     a2, a2, t5
    slli    t4, t2, 2
    beqz    a2, .Lip_f32_tail
.Lip_f32_loop:
    vle32.v v16, (a0)
    vle32.v v20, (a1)
    add     a0, a0, t4
    add     a1, a1, t4
    sub     a2, a2, t2
    vfmacc.vv v8, v16, v20
    bnez    a2, .Lip_f32_loop
.Lip_f32_tail:
    beqz    t5, .Lip_f32_red
    vsetvli zero, t5, e32, m4, tu, ma
    vle32.v v16, (a0)
    vle32.v v20, (a1)
    vfmacc.vv v8, v16, v20
.Lip_f32_red:
    vsetvli zero, t2, e32, m4, ta, ma
    vfredusum.vs v0, v8, v0
    vfmv.f.s fa0, v0
    ret

.Lip_f32_one_m1:
    vmv.s.x v0, zero                    /* zero v0 at vl>0 so len==0 reduces to 0 */
    vsetvli zero, a2, e32, m1, ta, ma
    vle32.v v8, (a0)
    vle32.v v9, (a1)
    vfmul.vv v8, v8, v9
    vfredusum.vs v0, v8, v0
    vfmv.f.s fa0, v0
    ret

.Lip_f32_one_m2:
    vsetvli zero, a2, e32, m2, ta, ma
    vle32.v v8, (a0)
    vle32.v v10, (a1)
    vfmul.vv v8, v8, v10
    vmv.s.x v0, zero
    vfredusum.vs v0, v8, v0
    vfmv.f.s fa0, v0
    ret
    .size spx_resample_rvv_ip_f32, . - spx_resample_rvv_ip_f32

/* ------------------------------------------------------------------ */
/* double spx_resample_rvv_ipd_f32(const float *a, const float *b,
 *                                 unsigned int len)
 * f32 multiply (product rounded to f32, like the C reference), f64
 * widening add into one m4 accumulator group.                          */
    .balign 4
    .globl  spx_resample_rvv_ipd_f32
    .hidden spx_resample_rvv_ipd_f32
    .type   spx_resample_rvv_ipd_f32, @function
spx_resample_rvv_ipd_f32:
    vsetvli t2, zero, e32, m2, ta, ma   /* t2 = vlmax(e32m2) == vlmax(e64m4) */
    vsetvli t1, zero, e64, m4, ta, ma
    vmv.v.i v8, 0                       /* f64m4 accumulator = +0.0 */
    vmv.s.x v0, zero
    remu    t5, a2, t2
    sub     a2, a2, t5
    slli    t4, t2, 2
    vsetvli zero, t2, e32, m2, ta, ma
    beqz    a2, .Lipd_tail
.Lipd_loop:
    vle32.v v16, (a0)
    vle32.v v18, (a1)
    add     a0, a0, t4
    add     a1, a1, t4
    sub     a2, a2, t2
    vfmul.vv v16, v16, v18
    vfwadd.wv v8, v8, v16
    bnez    a2, .Lipd_loop
.Lipd_tail:
    beqz    t5, .Lipd_reduce
    vsetvli zero, t5, e32, m2, tu, ma
    vle32.v v16, (a0)
    vle32.v v18, (a1)
    vfmul.vv v16, v16, v18
    vfwadd.wv v8, v8, v16
.Lipd_reduce:
    vsetvli zero, t2, e64, m4, ta, ma
    vfredusum.vs v0, v8, v0
    vfmv.f.s fa0, v0
    ret
    .size spx_resample_rvv_ipd_f32, . - spx_resample_rvv_ipd_f32

/* ------------------------------------------------------------------ */
/* float spx_resample_rvv_interp_f32(const float *a, const float *b,
 *          unsigned int len, unsigned int oversample, const float *frac)
 * ret = sum_k frac[k] * sum_j a[j] * b[j*oversample + k]. One vlsseg4
 * gathers all four taps; four m2 accumulators. Full strips run at a
 * constant vl == vlmax (strip-size address increments hoisted), one
 * tail-undisturbed partial strip handles len % vlmax. The frac combine is
 * fused into the epilogue (vfmacc.vf), so there is a single reduction
 * instead of four -- numerically a per-lane FMA combine instead of
 * scale-after-reduce, well inside the float tolerances.               */
    .balign 4
    .globl  spx_resample_rvv_interp_f32
    .hidden spx_resample_rvv_interp_f32
    .type   spx_resample_rvv_interp_f32, @function
spx_resample_rvv_interp_f32:
    slli    t6, a3, 2                   /* byte stride = oversample * 4 */
    flw     ft0, 0(a4)
    flw     ft1, 4(a4)
    flw     ft2, 8(a4)
    flw     ft3, 12(a4)
    vsetvli t2, zero, e32, m2, ta, ma   /* t2 = vlmax(e32m2) */
    vsetvli t1, zero, e32, m8, ta, ma
    vmv.v.i v8, 0                       /* zero v8..v15: all four m2 accs */
    vmv.s.x v0, zero
    remu    t5, a2, t2
    sub     a2, a2, t5
    slli    t4, t2, 2                   /* a-bytes per full strip */
    mul     t3, t2, t6                  /* b-bytes per full strip */
    vsetvli zero, t2, e32, m2, ta, ma   /* vl = vlmax for the full strips */
    beqz    a2, .Lif_tail
.Lif_loop:
    vle32.v v4, (a0)
    vlsseg4e32.v v16, (a1), t6          /* taps: v16v17 v18v19 v20v21 v22v23 */
    add     a0, a0, t4
    add     a1, a1, t3
    sub     a2, a2, t2
    vfmacc.vv v8,  v4, v16
    vfmacc.vv v10, v4, v18
    vfmacc.vv v12, v4, v20
    vfmacc.vv v14, v4, v22
    bnez    a2, .Lif_loop
.Lif_tail:
    beqz    t5, .Lif_reduce
    vsetvli zero, t5, e32, m2, tu, ma
    vle32.v v4, (a0)
    vlsseg4e32.v v16, (a1), t6
    vfmacc.vv v8,  v4, v16
    vfmacc.vv v10, v4, v18
    vfmacc.vv v12, v4, v20
    vfmacc.vv v14, v4, v22
.Lif_reduce:
    vsetvli zero, t2, e32, m2, ta, ma
    vfmul.vf  v8, v8, ft0
    vfmacc.vf v8, ft1, v10
    vfmacc.vf v8, ft2, v12
    vfmacc.vf v8, ft3, v14
    vfredusum.vs v0, v8, v0
    vfmv.f.s fa0, v0
    ret
    .size spx_resample_rvv_interp_f32, . - spx_resample_rvv_interp_f32

/* ------------------------------------------------------------------ */
/* double spx_resample_rvv_interpd_f32(const float *a, const float *b,
 *          unsigned int len, unsigned int oversample, const float *frac)
 * Same gather; vfwmacc (exact f64 products) into four f64m2 accumulators,
 * frac combine fused into the f64 epilogue.                            */
    .balign 4
    .globl  spx_resample_rvv_interpd_f32
    .hidden spx_resample_rvv_interpd_f32
    .type   spx_resample_rvv_interpd_f32, @function
spx_resample_rvv_interpd_f32:
    slli    t6, a3, 2
    flw     ft0, 0(a4)
    flw     ft1, 4(a4)
    flw     ft2, 8(a4)
    flw     ft3, 12(a4)
    fcvt.d.s ft0, ft0
    fcvt.d.s ft1, ft1
    fcvt.d.s ft2, ft2
    fcvt.d.s ft3, ft3
    vsetvli t2, zero, e32, m1, ta, ma   /* t2 = vlmax(e32m1) == vlmax(e64m2) */
    vsetvli t1, zero, e64, m8, ta, ma
    vmv.v.i v8, 0                       /* zero v8..v15: four f64m2 accs */
    vmv.s.x v0, zero
    remu    t5, a2, t2
    sub     a2, a2, t5
    slli    t4, t2, 2
    mul     t3, t2, t6
    vsetvli zero, t2, e32, m1, ta, ma
    beqz    a2, .Lid_tail
.Lid_loop:
    vle32.v v4, (a0)
    vlsseg4e32.v v16, (a1), t6          /* taps: v16 v17 v18 v19 */
    add     a0, a0, t4
    add     a1, a1, t3
    sub     a2, a2, t2
    vfwmacc.vv v8,  v4, v16
    vfwmacc.vv v10, v4, v17
    vfwmacc.vv v12, v4, v18
    vfwmacc.vv v14, v4, v19
    bnez    a2, .Lid_loop
.Lid_tail:
    beqz    t5, .Lid_reduce
    vsetvli zero, t5, e32, m1, tu, ma
    vle32.v v4, (a0)
    vlsseg4e32.v v16, (a1), t6
    vfwmacc.vv v8,  v4, v16
    vfwmacc.vv v10, v4, v17
    vfwmacc.vv v12, v4, v18
    vfwmacc.vv v14, v4, v19
.Lid_reduce:
    vsetvli zero, t2, e64, m2, ta, ma
    vfmul.vf  v8, v8, ft0
    vfmacc.vf v8, ft1, v10
    vfmacc.vf v8, ft2, v12
    vfmacc.vf v8, ft3, v14
    vfredusum.vs v0, v8, v0
    vfmv.f.s fa0, v0
    ret
    .size spx_resample_rvv_interpd_f32, . - spx_resample_rvv_interpd_f32

/* ------------------------------------------------------------------ */
/* spx_int32_t spx_resample_rvv_ip_i16(const spx_int16_t *a,
 *          const spx_int16_t *b, unsigned int len)
 * Raw wrapping-int32 widening-MAC dot product (exact products, wrap-around
 * accumulation: associative under any lane split, so bit-exact vs the C
 * loop). The caller applies SATURATE32PSHR.                            */
    .balign 4
    .globl  spx_resample_rvv_ip_i16
    .hidden spx_resample_rvv_ip_i16
    .type   spx_resample_rvv_ip_i16, @function
spx_resample_rvv_ip_i16:
    vsetvli t0, zero, e16, m1, ta, ma
    bgeu    t0, a2, .Lip_i16_one        /* len <= vlmax(e16m1): one strip */

    vsetvli t2, zero, e32, m4, ta, ma   /* t2 = vlmax(e32m4) == vlmax(e16m2) */
    vmv.v.i v8, 0
    vmv.s.x v0, zero
    remu    t5, a2, t2
    sub     a2, a2, t5
    slli    t4, t2, 1
    vsetvli zero, t2, e16, m2, ta, ma
    beqz    a2, .Lip_i16_tail
.Lip_i16_loop:
    vle16.v v16, (a0)
    vle16.v v18, (a1)
    add     a0, a0, t4
    add     a1, a1, t4
    sub     a2, a2, t2
    vwmacc.vv v8, v16, v18
    bnez    a2, .Lip_i16_loop
.Lip_i16_tail:
    beqz    t5, .Lip_i16_red
    vsetvli zero, t5, e16, m2, tu, ma
    vle16.v v16, (a0)
    vle16.v v18, (a1)
    vwmacc.vv v8, v16, v18
.Lip_i16_red:
    vsetvli zero, t2, e32, m4, ta, ma
    vredsum.vs v0, v8, v0
    vmv.x.s a0, v0
    ret

.Lip_i16_one:
    vsetivli zero, 1, e32, m1, ta, ma   /* zero v0 at vl>0 so len==0 reduces to 0 */
    vmv.s.x v0, zero
    vsetvli t0, a2, e16, m1, ta, ma
    vle16.v v8, (a0)
    vle16.v v9, (a1)
    vwmul.vv v10, v8, v9                /* exact i32 products, m2 */
    vsetvli zero, t0, e32, m2, ta, ma
    vredsum.vs v0, v10, v0
    vmv.x.s a0, v0
    ret
    .size spx_resample_rvv_ip_i16, . - spx_resample_rvv_ip_i16

/* ------------------------------------------------------------------ */
/* void spx_resample_rvv_interp4_i16(const spx_int16_t *a,
 *          const spx_int16_t *b, unsigned int len,
 *          unsigned int oversample, spx_int32_t accum[4])
 * Bit-exact four-way widening-MAC gather; the caller does the scalar
 * MULT16_32_Q15 / SATURATE32PSHR combine. Full strips at constant
 * vl == vlmax (any strip partition is bit-exact: wrapping int32
 * accumulation is associative), one _tu tail strip.                    */
    .balign 4
    .globl  spx_resample_rvv_interp4_i16
    .hidden spx_resample_rvv_interp4_i16
    .type   spx_resample_rvv_interp4_i16, @function
spx_resample_rvv_interp4_i16:
    slli    t6, a3, 1                   /* byte stride = oversample * 2 */
    vsetvli t2, zero, e16, m1, ta, ma   /* t2 = vlmax(e16m1) == vlmax(e32m2) */
    vsetvli t1, zero, e32, m8, ta, ma
    vmv.v.i v8, 0                       /* zero v8..v15: four i32m2 accs */
    vmv.s.x v0, zero
    remu    t5, a2, t2
    sub     a2, a2, t5
    slli    t4, t2, 1
    mul     t3, t2, t6
    vsetvli zero, t2, e16, m1, ta, ma
    beqz    a2, .Lii_tail
.Lii_loop:
    vle16.v v4, (a0)
    vlsseg4e16.v v16, (a1), t6          /* taps: v16 v17 v18 v19 */
    add     a0, a0, t4
    add     a1, a1, t3
    sub     a2, a2, t2
    vwmacc.vv v8,  v4, v16
    vwmacc.vv v10, v4, v17
    vwmacc.vv v12, v4, v18
    vwmacc.vv v14, v4, v19
    bnez    a2, .Lii_loop
.Lii_tail:
    beqz    t5, .Lii_reduce
    vsetvli zero, t5, e16, m1, tu, ma
    vle16.v v4, (a0)
    vlsseg4e16.v v16, (a1), t6
    vwmacc.vv v8,  v4, v16
    vwmacc.vv v10, v4, v17
    vwmacc.vv v12, v4, v18
    vwmacc.vv v14, v4, v19
.Lii_reduce:
    vsetvli zero, t2, e32, m2, ta, ma
    vredsum.vs v1, v8,  v0
    vredsum.vs v2, v10, v0
    vredsum.vs v3, v12, v0
    vredsum.vs v4, v14, v0
    vmv.x.s t0, v1
    sw      t0, 0(a4)
    vmv.x.s t0, v2
    sw      t0, 4(a4)
    vmv.x.s t0, v3
    sw      t0, 8(a4)
    vmv.x.s t0, v4
    sw      t0, 12(a4)
    ret
    .size spx_resample_rvv_interp4_i16, . - spx_resample_rvv_interp4_i16

#if defined(__ELF__)
/* Non-executable stack note */
    .section .note.GNU-stack, "", @progbits
#endif
