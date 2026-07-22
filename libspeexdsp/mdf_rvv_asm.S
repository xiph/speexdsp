/* Copyright (C) 2003-2008 Jean-Marc Valin
 * Copyright (C) 2026 Tristan Matthews
 */
/*
   mdf_rvv_asm.S -- hand-written RVV 1.0 kernels for the MDF echo canceller.

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

/* RVV 1.0 kernels (VLEN >= 128) for mdf.c's per-frame O(M*N) spectral
 * loops. All of them walk the half-complex (packed) spectrum layout
 * [DC, r1,i1, ..., rK,iK, Nyquist]: the interior complex bins live at odd
 * offsets, so every kernel starts one element in, loads (re,im) pairs with
 * vlseg2, and leaves the two real-only edge bins (index 0 and N-1) to the
 * C wrapper in mdf_rvv.h. Bins strip-mine with LMUL=4; the M filter blocks
 * (N elements apart) accumulate into vector registers inside the strip.
 *
 * Fixed point matches C bit for bit: vwmul/vwmacc exact e16xe16->e32
 * products with mod-2^32 accumulation (== the int32 wrap the C MAC16_16
 * chains compile to), vadd+vsra for PSHR32's add-then-shift rounding (a
 * plain vssra would differ from C when the rounding add wraps), truncating
 * vnsrl/vnsra narrows (== the wrapping int16 casts / EXTRACT16(SHR32)).
 * Wrapping adds commute, so splitting tmp1 into separate +/- accumulators
 * preserves bit-exactness. Float kernels use FMAs and reordered sums, so
 * they pair with a small checkasm tolerance.
 *
 * Reduction kernels accumulate across strips in a register group, so their
 * strip-body vsetvlis use tu (a shorter tail strip must not clobber the
 * partial sums sitting in the tail lanes); everything else is ta.
 *
 * `.option arch, +v` keeps the file self-contained (no -march needed).
 * Leaf functions, caller-saved registers only. The wrappers guarantee
 * N >= 4 even, M >= 1, len >= 0 and in-bounds loads. */

    .option arch, +v
    .text

/* ------------------------------------------------------------------ */
/* unsigned int spx_mdf_rvv_vlenb(void) -- VLEN in bytes. */
    .balign 4
    .globl  spx_mdf_rvv_vlenb
    .hidden spx_mdf_rvv_vlenb
    .type   spx_mdf_rvv_vlenb, @function
spx_mdf_rvv_vlenb:
    vsetvli a0, zero, e8, m1, ta, ma
    ret
    .size spx_mdf_rvv_vlenb, . - spx_mdf_rvv_vlenb

/* ------------------------------------------------------------------ */
/* int spx_mdf_rvv_compliant(void)
 * Reject draft RVV 0.7.1 hardware (sets VILL on a ta,ma vtype config).
 * Mirrors spx_resample_rvv_compliant. */
    .balign 4
    .globl  spx_mdf_rvv_compliant
    .hidden spx_mdf_rvv_compliant
    .type   spx_mdf_rvv_compliant, @function
spx_mdf_rvv_compliant:
    vsetvli t0, zero, e8, m1, ta, ma
    csrr    a0, vtype
    sgtz    a0, a0
    ret
    .size spx_mdf_rvv_compliant, . - spx_mdf_rvv_compliant

/* ================================================================== */
/*                        float (f32) kernels                         */
/* ================================================================== */

/* ------------------------------------------------------------------ */
/* void spx_mdf_rvv_smul_accum_f32(const float *X, const float *Y,
 *                                 float *acc, int N, int M)
 * acc[1..N-2] = sum over M blocks of the complex product X*Y (interior
 * bins only). Registers: X (v8,v12) Y (v16,v20) accr v24 acci v28. */
    .balign 4
    .globl  spx_mdf_rvv_smul_accum_f32
    .hidden spx_mdf_rvv_smul_accum_f32
    .type   spx_mdf_rvv_smul_accum_f32, @function
spx_mdf_rvv_smul_accum_f32:
    slli    t4, a3, 2              /* block stride = N floats */
    addi    a0, a0, 4              /* skip the DC bin */
    addi    a1, a1, 4
    addi    a2, a2, 4
    addi    t5, a3, -2
    srai    t5, t5, 1              /* interior bin count */
1:  /* strip loop */
    vsetvli t0, t5, e32, m4, ta, ma
    vmv.v.i v24, 0                 /* accr */
    vmv.v.i v28, 0                 /* acci */
    mv      t1, a0
    mv      t2, a1
    mv      t3, a4
2:  /* block loop */
    vlseg2e32.v v8, (t1)           /* xr v8, xi v12 */
    vlseg2e32.v v16, (t2)          /* yr v16, yi v20 */
    add     t1, t1, t4
    add     t2, t2, t4
    vfmacc.vv  v24, v8, v16        /* accr += xr*yr */
    vfnmsac.vv v24, v12, v20       /* accr -= xi*yi */
    vfmacc.vv  v28, v12, v16       /* acci += xi*yr */
    vfmacc.vv  v28, v8, v20        /* acci += xr*yi */
    addi    t3, t3, -1
    bnez    t3, 2b
    vsseg2e32.v v24, (a2)
    slli    t6, t0, 3              /* vl bins = vl*8 bytes */
    add     a0, a0, t6
    add     a1, a1, t6
    add     a2, a2, t6
    sub     t5, t5, t0
    bnez    t5, 1b
    ret
    .size spx_mdf_rvv_smul_accum_f32, . - spx_mdf_rvv_smul_accum_f32

/* ------------------------------------------------------------------ */
/* void spx_mdf_rvv_wsmul_conj_f32(const float *w, const float *X,
 *                                 const float *Y, float *prod, int N,
 *                                 float p)
 * prod[1..N-2]: per interior bin j, W = p*w[j];
 *   prod_r = W*(xr*yr + xi*yi), prod_i = W*(xr*yi - xi*yr). */
    .balign 4
    .globl  spx_mdf_rvv_wsmul_conj_f32
    .hidden spx_mdf_rvv_wsmul_conj_f32
    .type   spx_mdf_rvv_wsmul_conj_f32, @function
spx_mdf_rvv_wsmul_conj_f32:
    addi    a0, a0, 4              /* w[1] (unit stride) */
    addi    a1, a1, 4
    addi    a2, a2, 4
    addi    a3, a3, 4
    addi    t5, a4, -2
    srai    t5, t5, 1
1:
    vsetvli t0, t5, e32, m4, ta, ma
    vle32.v v4, (a0)
    vfmul.vf v4, v4, fa0           /* W = p*w[j] */
    vlseg2e32.v v8, (a1)           /* xr v8, xi v12 */
    vlseg2e32.v v16, (a2)          /* yr v16, yi v20 */
    vfmul.vv   v24, v8, v16        /* t1 = xr*yr */
    vfmacc.vv  v24, v12, v20       /* t1 += xi*yi */
    vfmul.vv   v28, v8, v20        /* t2 = xr*yi */
    vfnmsac.vv v28, v12, v16       /* t2 -= xi*yr */
    vfmul.vv   v24, v24, v4
    vfmul.vv   v28, v28, v4
    vsseg2e32.v v24, (a3)
    slli    t6, t0, 3
    slli    t1, t0, 2
    add     a0, a0, t1
    add     a1, a1, t6
    add     a2, a2, t6
    add     a3, a3, t6
    sub     t5, t5, t0
    bnez    t5, 1b
    ret
    .size spx_mdf_rvv_wsmul_conj_f32, . - spx_mdf_rvv_wsmul_conj_f32

/* ------------------------------------------------------------------ */
/* power_spectrum(_accum): ps[1..(N-2)/2] (+)= xr^2 + xi^2, interior bins. */
.macro mdf_f32_power_spectrum name, accum
    .balign 4
    .globl  \name
    .hidden \name
    .type   \name, @function
\name:
    addi    a0, a0, 4
    addi    a1, a1, 4
    addi    t5, a2, -2
    srai    t5, t5, 1
1:
    vsetvli t0, t5, e32, m4, ta, ma
    vlseg2e32.v v8, (a0)           /* xr v8, xi v12 */
.if \accum
    vle32.v v16, (a1)
    vfmacc.vv v16, v8, v8
    vfmacc.vv v16, v12, v12
.else
    vfmul.vv  v16, v8, v8
    vfmacc.vv v16, v12, v12
.endif
    vse32.v v16, (a1)
    slli    t6, t0, 3
    add     a0, a0, t6
    slli    t6, t0, 2
    add     a1, a1, t6
    sub     t5, t5, t0
    bnez    t5, 1b
    ret
    .size \name, . - \name
.endm

/* void spx_mdf_rvv_power_spectrum_f32(const float *X, float *ps, int N) */
mdf_f32_power_spectrum spx_mdf_rvv_power_spectrum_f32, 0
/* void spx_mdf_rvv_power_spectrum_accum_f32(const float *X, float *ps, int N) */
mdf_f32_power_spectrum spx_mdf_rvv_power_spectrum_accum_f32, 1

/* ------------------------------------------------------------------ */
/* float spx_mdf_rvv_inner_prod_f32(const float *x, const float *y, int len)
 * Plain dot product (mdf_inner_prod's shifts are no-ops in float; also
 * reused with x==y for mdf_adjust_prop's sum of squares). */
    .balign 4
    .globl  spx_mdf_rvv_inner_prod_f32
    .hidden spx_mdf_rvv_inner_prod_f32
    .type   spx_mdf_rvv_inner_prod_f32, @function
spx_mdf_rvv_inner_prod_f32:
    vsetvli t1, zero, e32, m4, ta, ma
    vmv.v.i v24, 0
    beqz    a2, 2f
1:  /* tu: tail lanes hold live partial sums from fuller strips */
    vsetvli t0, a2, e32, m4, tu, ma
    vle32.v v8, (a0)
    vle32.v v16, (a1)
    vfmacc.vv v24, v8, v16
    slli    t2, t0, 2
    add     a0, a0, t2
    add     a1, a1, t2
    sub     a2, a2, t0
    bnez    a2, 1b
2:
    vsetvli t0, zero, e32, m4, ta, ma
    vmv.s.x v0, zero
    vfredusum.vs v0, v24, v0
    vfmv.f.s fa0, v0
    ret
    .size spx_mdf_rvv_inner_prod_f32, . - spx_mdf_rvv_inner_prod_f32

/* ================================================================== */
/*                     fixed-point (i16) kernels                      */
/* ================================================================== */

/* ------------------------------------------------------------------ */
/* spectral_mul_accum(16): acc[1..N-2] = (int16)PSHR32(sum, shift) of the
 * complex products over M blocks. X is int16; Y is int32 (top halves
 * used, is32=1) or int16 (is32=0). tmp1's interleaved add/sub chain is
 * split into separate + and - accumulators (wrapping adds commute):
 *   P v20 = sum xr*tyr, Q v24 = sum xi*tyi, R v28 = sum xi*tyr + xr*tyi
 * Loads: xr v4, xi v6 (e16m2); Y32 v8/v12 (e32m4); tyr v16, tyi v18. */
.macro mdf_i16_smul_accum name, is32
    .balign 4
    .globl  \name
    .hidden \name
    .type   \name, @function
\name:
    slli    t4, a3, 1              /* X block stride = N*2 bytes */
.if \is32
    slli    t5, a3, 2              /* Y block stride = N*4 bytes */
    addi    a1, a1, 4
.else
    mv      t5, t4
    addi    a1, a1, 2
.endif
    addi    a0, a0, 2
    addi    a2, a2, 2
    addi    t6, a3, -2
    srai    t6, t6, 1              /* interior bin count */
    addi    t0, a5, -1
    li      t1, 1
    sll     a6, t1, t0             /* PSHR32 rounding bias (shift >= 1) */
1:  /* strip loop */
    vsetvli t0, t6, e32, m4, ta, ma
    vmv.v.i v20, 0                 /* P */
    vmv.v.i v24, 0                 /* Q */
    vmv.v.i v28, 0                 /* R */
    vsetvli zero, zero, e16, m2, ta, ma
    mv      t1, a0
    mv      t2, a1
    mv      t3, a4
2:  /* block loop */
    vlseg2e16.v v4, (t1)           /* xr v4, xi v6 */
.if \is32
    vlseg2e32.v v8, (t2)           /* yr v8m4, yi v12m4 */
    vnsrl.wi v16, v8, 16           /* tyr = TOP16(yr) */
    vnsrl.wi v18, v12, 16          /* tyi = TOP16(yi) */
.else
    vlseg2e16.v v16, (t2)          /* tyr v16, tyi v18 */
.endif
    add     t1, t1, t4
    add     t2, t2, t5
    vwmacc.vv v20, v4, v16         /* P += xr*tyr */
    vwmacc.vv v24, v6, v18         /* Q += xi*tyi */
    vwmacc.vv v28, v6, v16         /* R += xi*tyr */
    vwmacc.vv v28, v4, v18         /* R += xr*tyi */
    addi    t3, t3, -1
    bnez    t3, 2b
    vsetvli zero, zero, e32, m4, ta, ma
    vsub.vv v20, v20, v24          /* tmp1 = P - Q */
    vadd.vx v20, v20, a6           /* PSHR32: wrapping add-then-shift, */
    vsra.vx v20, v20, a5           /* exactly as the C compiles */
    vadd.vx v28, v28, a6
    vsra.vx v28, v28, a5
    vsetvli zero, zero, e16, m2, ta, ma
    vnsrl.wi v4, v20, 0            /* truncating int16 casts */
    vnsrl.wi v6, v28, 0
    vsseg2e16.v v4, (a2)
    slli    t1, t0, 2              /* vl bins = vl*4 bytes of X/acc */
    add     a0, a0, t1
    add     a2, a2, t1
.if \is32
    slli    t1, t0, 3
.endif
    add     a1, a1, t1
    sub     t6, t6, t0
    bnez    t6, 1b
    ret
    .size \name, . - \name
.endm

/* void spx_mdf_rvv_smul_accum_i16(const spx_int16_t *X, const spx_int32_t *Y,
 *                                 spx_int16_t *acc, int N, int M, int shift) */
mdf_i16_smul_accum spx_mdf_rvv_smul_accum_i16, 1
/* void spx_mdf_rvv_smul_accum16_i16(const spx_int16_t *X, const spx_int16_t *Y,
 *                                   spx_int16_t *acc, int N, int M, int shift) */
mdf_i16_smul_accum spx_mdf_rvv_smul_accum16_i16, 0

/* ------------------------------------------------------------------ */
/* power_spectrum(_accum): ps[1..(N-2)/2] (+)= xr^2 + xi^2 in int32,
 * interior bins, bit-exact (exact products, wrapping adds). */
.macro mdf_i16_power_spectrum name, accum
    .balign 4
    .globl  \name
    .hidden \name
    .type   \name, @function
\name:
    addi    a0, a0, 2
    addi    a1, a1, 4
    addi    t5, a2, -2
    srai    t5, t5, 1
1:
    vsetvli t0, t5, e16, m2, ta, ma
    vlseg2e16.v v4, (a0)           /* xr v4, xi v6 */
    vwmul.vv  v8, v4, v4           /* e32m4 */
    vwmacc.vv v8, v6, v6
.if \accum
    vle32.v v12, (a1)
    vsetvli zero, zero, e32, m4, ta, ma
    vadd.vv v8, v8, v12
.endif
    vse32.v v8, (a1)
    slli    t6, t0, 2
    add     a0, a0, t6
    add     a1, a1, t6
    sub     t5, t5, t0
    bnez    t5, 1b
    ret
    .size \name, . - \name
.endm

/* void spx_mdf_rvv_power_spectrum_i16(const spx_int16_t *X, spx_int32_t *ps, int N) */
mdf_i16_power_spectrum spx_mdf_rvv_power_spectrum_i16, 0
/* void spx_mdf_rvv_power_spectrum_accum_i16(const spx_int16_t *X, spx_int32_t *ps, int N) */
mdf_i16_power_spectrum spx_mdf_rvv_power_spectrum_accum_i16, 1

/* ------------------------------------------------------------------ */
/* spx_int32_t spx_mdf_rvv_inner_prod_i16(const spx_int16_t *x,
 *                                        const spx_int16_t *y,
 *                                        int len, int shift)
 * mdf_inner_prod: per pair, part = x0*y0 + x1*y1 (exact e32), then
 * sum += part >> shift, all wrapping -> order-free, bit-exact. An odd
 * trailing element is dropped, as in C. */
    .balign 4
    .globl  spx_mdf_rvv_inner_prod_i16
    .hidden spx_mdf_rvv_inner_prod_i16
    .type   spx_mdf_rvv_inner_prod_i16, @function
spx_mdf_rvv_inner_prod_i16:
    srai    a2, a2, 1              /* pair count */
    vsetvli t1, zero, e32, m4, ta, ma
    vmv.v.i v24, 0
    beqz    a2, 2f
1:
    vsetvli t0, a2, e16, m2, tu, ma
    vlseg2e16.v v4, (a0)           /* xe v4, xo v6 */
    vlseg2e16.v v8, (a1)           /* ye v8, yo v10 */
    vwmul.vv  v12, v4, v8          /* part e32m4 */
    vwmacc.vv v12, v6, v10
    vsetvli zero, zero, e32, m4, tu, ma
    vsra.vx v12, v12, a3
    vadd.vv v24, v24, v12
    slli    t2, t0, 2              /* vl pairs = vl*4 bytes */
    add     a0, a0, t2
    add     a1, a1, t2
    sub     a2, a2, t0
    bnez    a2, 1b
2:
    vsetvli t0, zero, e32, m4, ta, ma
    vmv.s.x v0, zero
    vredsum.vs v0, v24, v0
    vmv.x.s a0, v0
    ret
    .size spx_mdf_rvv_inner_prod_i16, . - spx_mdf_rvv_inner_prod_i16

/* ------------------------------------------------------------------ */
/* spx_int32_t spx_mdf_rvv_prop_sumsq_i16(const spx_int32_t *W, int len,
 *                                        int shift)
 * mdf_adjust_prop's inner sum: sq(EXTRACT16(SHR32(W, shift))) summed
 * wrapping. vnsra = arithmetic shift + truncating narrow in one op. */
    .balign 4
    .globl  spx_mdf_rvv_prop_sumsq_i16
    .hidden spx_mdf_rvv_prop_sumsq_i16
    .type   spx_mdf_rvv_prop_sumsq_i16, @function
spx_mdf_rvv_prop_sumsq_i16:
    vsetvli t1, zero, e32, m4, ta, ma
    vmv.v.i v24, 0
    beqz    a1, 2f
1:
    vsetvli t0, a1, e16, m2, tu, ma
    vle32.v v8, (a0)               /* e32m4 */
    vnsra.wx v4, v8, a2            /* EXTRACT16(SHR32(W, shift)) */
    vwmacc.vv v24, v4, v4
    slli    t2, t0, 2
    add     a0, a0, t2
    sub     a1, a1, t0
    bnez    a1, 1b
2:
    vsetvli t0, zero, e32, m4, ta, ma
    vmv.s.x v0, zero
    vredsum.vs v0, v24, v0
    vmv.x.s a0, v0
    ret
    .size spx_mdf_rvv_prop_sumsq_i16, . - spx_mdf_rvv_prop_sumsq_i16

#if defined(__ELF__)
/* Non-executable stack note */
    .section .note.GNU-stack, "", @progbits
#endif
