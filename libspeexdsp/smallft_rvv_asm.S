/* Copyright (C) 1994-1996 Paul N. Swarztrauber (FFTPACK)
 * Copyright (C) The OggSquish source code, Xiph.Org
 * Copyright (C) 2026 Tristan Matthews
 */
/*
   smallft_rvv_asm.S -- hand-written RVV 1.0 kernels for the smallft
   (FFTPACK real FFT) radix-4 and radix-2 stages.

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

/* RVV 1.0 kernels (VLEN >= 128) for dradf4/dradb4 (full replacements) and
 * dradf2/dradb2 (the l1==1 shape they always take here, since drfti1 moves
 * the single factor 2 to the front of ifac). Each keeps FFTPACK's
 * three-part structure:
 *
 *  - part 1 (the first element of each of the l1 rows) and part 3 (the
 *    even-ido Nyquist column) vectorize across the rows with strided
 *    accesses; the ido==1 shape of part 1 (a quarter of a 4^k transform's
 *    work) gets a unit-stride vlseg4/vle fast path.
 *  - part 2 (the twiddled interior columns, (ido-2)/2 complex pairs per
 *    row) picks its axis per shape like the kiss bfly kernels: across the
 *    pair index when npairs >= l1 (vector twiddles via vlseg2, the
 *    mirrored column walked with negative-stride segment accesses), else
 *    across the rows (scalar twiddles broadcast per column, everything a
 *    positive-stride segment access).
 *
 * cc and ch never alias (drftf1/drftb1 ping-pong two buffers), so parts
 * read only cc and write only ch. Float only -- smallft is only built for
 * FLOATING_POINT -- but the kernels take pointer/integer args exclusively
 * and build float constants with vmv.v.x bit patterns, so they carry no
 * float-ABI dependency (mirroring the kiss kernels). Segment ops keep
 * NFIELDS*EMUL <= 8: seg2 runs at LMUL <= 2 (or m1 across rows), seg4 at
 * LMUL=2.
 *
 * `.option arch, +v` keeps the file self-contained (no -march needed).
 * Callers guarantee ido >= 1, l1 >= 1 and in-bounds twiddle loads. */

    .option arch, +v
    .text

#define HSQT2_BITS  0x3F3504F3   /* 0.70710678118654752f */
#define NHSQT2_BITS 0xBF3504F3
#define SQRT2_BITS  0x3FB504F3   /* 1.414213562373095f */
#define NSQRT2_BITS 0xBFB504F3

/* ------------------------------------------------------------------ */
/* unsigned int spx_drft_rvv_vlenb(void) -- VLEN in bytes. */
    .balign 4
    .globl  spx_drft_rvv_vlenb
    .hidden spx_drft_rvv_vlenb
    .type   spx_drft_rvv_vlenb, @function
spx_drft_rvv_vlenb:
    vsetvli a0, zero, e8, m1, ta, ma
    ret
    .size spx_drft_rvv_vlenb, . - spx_drft_rvv_vlenb

/* ------------------------------------------------------------------ */
/* int spx_drft_rvv_compliant(void)
 * Reject draft RVV 0.7.1 hardware (sets VILL on a ta,ma vtype config). */
    .balign 4
    .globl  spx_drft_rvv_compliant
    .hidden spx_drft_rvv_compliant
    .type   spx_drft_rvv_compliant, @function
spx_drft_rvv_compliant:
    vsetvli t0, zero, e8, m1, ta, ma
    csrr    a0, vtype
    sgtz    a0, a0
    ret
    .size spx_drft_rvv_compliant, . - spx_drft_rvv_compliant

/* Forward complex twiddle (conjugate style used by dradf*):
 *   cr = wr*yr + wi*yi ; ci = wr*yi - wi*yr */
.macro drft_cmulf cr, ci, wr, wi, yr, yi
    vfmul.vv   \cr, \wr, \yr
    vfmacc.vv  \cr, \wi, \yi
    vfmul.vv   \ci, \wr, \yi
    vfnmsac.vv \ci, \wi, \yr
.endm

/* Backward complex twiddle (used by dradb*):
 *   or = wr*cr - wi*ci ; oi = wr*ci + wi*cr */
.macro drft_cmulb or, oi, wr, wi, cr, ci
    vfmul.vv   \or, \wr, \cr
    vfnmsac.vv \or, \wi, \ci
    vfmul.vv   \oi, \wr, \ci
    vfmacc.vv  \oi, \wi, \cr
.endm

/* ================================================================== */
/* void spx_drft_rvv_radf4_f32(int ido, int l1, const float *cc,
 *                             float *ch, const float *wa1,
 *                             const float *wa2, const float *wa3)      */
/* ================================================================== */
    .balign 4
    .globl  spx_drft_rvv_radf4_f32
    .hidden spx_drft_rvv_radf4_f32
    .type   spx_drft_rvv_radf4_f32, @function
spx_drft_rvv_radf4_f32:
    slli    t0, a0, 2              /* idoB = ido*4 */
    mul     a7, t0, a1             /* t0B = l1*ido*4 (plane stride) */

/* ---- part 1: rows' first elements, lanes = rows ----
 *  x0..x3 = cc[k*ido + m*t0]; ta = x1+x3, tb = x0+x2;
 *  ch[4k*ido + {0, 2ido-1, 2ido, 4ido-1}] = {ta+tb, x0-x2, x3-x1, tb-ta} */
    li      t1, 1
    bne     a0, t1, .Lf4_p1_strided

    /* ido==1: unit-stride loads, one vsseg4 (outputs are 4k..4k+3) */
    mv      t1, a1
    mv      t2, a2
    mv      t3, a3
1:
    vsetvli t4, t1, e32, m2, ta, ma
    add     t5, t2, a7
    vle32.v v0, (t2)               /* x0 */
    vle32.v v4, (t5)               /* x1 */
    add     t5, t5, a7
    vle32.v v8, (t5)               /* x2 */
    add     t5, t5, a7
    vle32.v v12, (t5)              /* x3 */
    vfadd.vv v24, v4, v12          /* ta */
    vfadd.vv v26, v0, v8           /* tb */
    vfadd.vv v16, v24, v26         /* f0 = ta+tb */
    vfsub.vv v18, v0, v8           /* f1 = x0-x2 */
    vfsub.vv v20, v12, v4          /* f2 = x3-x1 */
    vfsub.vv v22, v26, v24         /* f3 = tb-ta */
    vsseg4e32.v v16, (t3)
    slli    t5, t4, 2
    add     t2, t2, t5
    slli    t5, t4, 4
    add     t3, t3, t5
    sub     t1, t1, t4
    bnez    t1, 1b
    ret                            /* ido < 2: done */

.Lf4_p1_strided:
    mv      t1, a1
    mv      t2, a2
    mv      t3, a3
    slli    t6, t0, 2              /* 4*idoB = ch lane stride */
1:
    vsetvli t4, t1, e32, m4, ta, ma
    add     t5, t2, a7
    vlse32.v v0, (t2), t0          /* x0 */
    vlse32.v v4, (t5), t0          /* x1 */
    add     t5, t5, a7
    vlse32.v v8, (t5), t0          /* x2 */
    add     t5, t5, a7
    vlse32.v v12, (t5), t0         /* x3 */
    vfadd.vv v24, v4, v12          /* ta */
    vfadd.vv v28, v0, v8           /* tb */
    vfadd.vv v16, v24, v28         /* f0 */
    vfsub.vv v20, v0, v8           /* f1 */
    vfsub.vv v28, v28, v24         /* f3 (before ta dies) */
    vfsub.vv v24, v12, v4          /* f2 */
    vsse32.v v16, (t3), t6         /* @ +0 */
    slli    t5, t0, 1
    add     t5, t5, t3
    addi    t5, t5, -4
    vsse32.v v20, (t5), t6         /* @ +2ido-1 */
    addi    t5, t5, 4
    vsse32.v v24, (t5), t6         /* @ +2ido */
    slli    t5, t0, 2
    add     t5, t5, t3
    addi    t5, t5, -4
    vsse32.v v28, (t5), t6         /* @ +4ido-1 */
    mul     t5, t4, t0
    add     t2, t2, t5
    slli    t5, t5, 2
    add     t3, t3, t5
    sub     t1, t1, t4
    bnez    t1, 1b

    li      t1, 2
    beq     a0, t1, .Lf4_p3
    addi    t1, a0, -1
    srai    t1, t1, 1              /* npairs = (ido-1)>>1, both parities */

/* ---- part 2: twiddled interior columns ----
 * per (row k, pair i): x0 = cc[k*ido + i-1..i], ym = cc[.. + m*t0]
 * twiddled by wam -> crm/cim; radix-4 butterfly; 4 output pair streams,
 * two of them walking down from the row's mirror point 2ido-i. */
    addi    sp, sp, -64
    sd      s0, 0(sp)
    sd      s1, 8(sp)
    sd      s2, 16(sp)
    sd      s3, 24(sp)
    sd      s4, 32(sp)
    sd      s5, 40(sp)
    sd      a2, 48(sp)
    sd      a3, 56(sp)

    blt     t1, a1, .Lf4_p2_k      /* npairs < l1 -> vectorize across rows */

/* -- across pairs: outer row loop, strip-mined pairs (m2) -- */
    addi    s0, a2, 4              /* x0 base   = cc + 1 elt */
    addi    s1, a3, 4              /* outA base = ch + 1 elt */
    slli    t2, t0, 1
    add     s2, a3, t2
    addi    s2, s2, -12            /* outB base = ch + 2ido-3 elts */
    mv      s3, a1                 /* row counter */
    slli    s4, t0, 1              /* 2idoB (C/D offset from A/B) */
.Lf4_i_row:
    mv      t2, s0                 /* x0 walk */
    mv      t3, s1                 /* outA walk */
    mv      t4, s2                 /* outB walk (desc) */
    li      t5, 0                  /* twiddle byte offset */
    mv      t6, t1                 /* pairs remaining */
.Lf4_i_strip:
    vsetvli s5, t6, e32, m2, ta, ma
    vlseg2e32.v v0, (t2)           /* x0r v0, x0i v2 */
    add     a2, t2, a7
    vlseg2e32.v v4, (a2)           /* y1 */
    add     a2, a2, a7
    vlseg2e32.v v8, (a2)           /* y2 */
    add     a2, a2, a7
    vlseg2e32.v v12, (a2)          /* y3 */
    add     a2, a4, t5
    vlseg2e32.v v16, (a2)          /* w1 */
    drft_cmulf v20, v22, v16, v18, v4, v6    /* cr2/ci2 */
    add     a2, a5, t5
    vlseg2e32.v v16, (a2)          /* w2 */
    drft_cmulf v24, v26, v16, v18, v8, v10   /* cr3/ci3 */
    add     a2, a6, t5
    vlseg2e32.v v16, (a2)          /* w3 */
    drft_cmulf v28, v30, v16, v18, v12, v14  /* cr4/ci4 */
    vfadd.vv v4, v20, v28          /* tr1 = cr2+cr4 */
    vfsub.vv v6, v28, v20          /* tr4 = cr4-cr2 */
    vfadd.vv v8, v22, v30          /* ti1 = ci2+ci4 */
    vfsub.vv v10, v22, v30         /* ti4 = ci2-ci4 */
    vfadd.vv v12, v2, v26          /* ti2 = x0i+ci3 */
    vfsub.vv v14, v2, v26          /* ti3 = x0i-ci3 */
    vfadd.vv v16, v0, v24          /* tr2 = x0r+cr3 */
    vfsub.vv v18, v0, v24          /* tr3 = x0r-cr3 */
    vfadd.vv v20, v4, v16          /* outA = (tr1+tr2, ti1+ti2) */
    vfadd.vv v22, v8, v12
    vsseg2e32.v v20, (t3)
    vfadd.vv v20, v10, v18         /* outC = (ti4+tr3, tr4+ti3) */
    vfadd.vv v22, v6, v14
    add     a2, t3, s4
    vsseg2e32.v v20, (a2)
    /* Descending streams: reverse the lanes and store unit-stride from
     * the strip's lowest pair address (strided stores run element-at-a-
     * time on the K1; vrgather + vsseg is much cheaper). */
    vid.v   v28
    addi    a3, s5, -1
    vrsub.vx v28, v28, a3          /* lane reversal indices */
    vfsub.vv v20, v18, v10         /* outB = (tr3-ti4, tr4-ti3), desc */
    vfsub.vv v22, v6, v14
    vrgather.vv v24, v20, v28
    vrgather.vv v26, v22, v28
    slli    a3, s5, 3
    addi    a3, a3, -8
    sub     a2, t4, a3             /* lowest pair of the strip */
    vsseg2e32.v v24, (a2)
    vfsub.vv v20, v16, v4          /* outD = (tr2-tr1, ti1-ti2), desc */
    vfsub.vv v22, v8, v12
    vrgather.vv v24, v20, v28
    vrgather.vv v26, v22, v28
    add     a2, a2, s4
    vsseg2e32.v v24, (a2)
    slli    a2, s5, 3
    add     t2, t2, a2
    add     t3, t3, a2
    add     t5, t5, a2
    sub     t4, t4, a2
    sub     t6, t6, s5
    bnez    t6, .Lf4_i_strip
    add     s0, s0, t0             /* cc row stride ido */
    slli    a2, t0, 2
    add     s1, s1, a2             /* ch row stride 4ido */
    add     s2, s2, a2
    addi    s3, s3, -1
    bnez    s3, .Lf4_i_row
    j       .Lf4_p2_done

/* -- across rows: outer pair loop (m1, broadcast twiddles) --
 * cc lane stride idoB, ch lane stride 4*idoB; the mirrored streams sit
 * at fixed per-pair offsets, so all accesses stride forward. */
.Lf4_p2_k:
    mv      s0, t1                 /* pair countdown */
    addi    s1, a2, 4              /* cc pair base (+8 per pair) */
    addi    s2, a3, 4              /* outA base (+8 per pair) */
    slli    t2, t0, 1
    add     s3, a3, t2
    addi    s3, s3, -12            /* outB base (-8 per pair) */
    mv      s4, a4                 /* w1 walk */
    mv      s5, a5                 /* w2 walk (w3 walks in a6) */
    slli    t6, t0, 2              /* 4*idoB = ch lane stride */
.Lf4_k_pair:
    vsetvli t2, a1, e32, m1, ta, ma
    lw      a2, 0(s4)
    vmv.v.x v24, a2                /* w1r */
    lw      a2, 4(s4)
    vmv.v.x v25, a2                /* w1i */
    lw      a2, 0(s5)
    vmv.v.x v26, a2                /* w2r */
    lw      a2, 4(s5)
    vmv.v.x v27, a2                /* w2i */
    lw      a2, 0(a6)
    vmv.v.x v28, a2                /* w3r */
    lw      a2, 4(a6)
    vmv.v.x v29, a2                /* w3i */
    mv      t1, a1                 /* rows remaining */
    mv      t3, s1                 /* cc walk */
    mv      t4, s2                 /* outA walk */
    mv      t5, s3                 /* outB walk */
.Lf4_k_strip:
    vsetvli t2, t1, e32, m1, ta, ma
    vlsseg2e32.v v0, (t3), t0      /* x0r v0, x0i v1 */
    add     a2, t3, a7
    vlsseg2e32.v v2, (a2), t0      /* y1 */
    add     a2, a2, a7
    vlsseg2e32.v v4, (a2), t0      /* y2 */
    add     a2, a2, a7
    vlsseg2e32.v v6, (a2), t0      /* y3 */
    drft_cmulf v8, v9, v24, v25, v2, v3      /* cr2/ci2 */
    drft_cmulf v10, v11, v26, v27, v4, v5    /* cr3/ci3 */
    drft_cmulf v12, v13, v28, v29, v6, v7    /* cr4/ci4 */
    vfadd.vv v2, v8, v12           /* tr1 */
    vfsub.vv v3, v12, v8           /* tr4 */
    vfadd.vv v4, v9, v13           /* ti1 */
    vfsub.vv v5, v9, v13           /* ti4 */
    vfadd.vv v6, v1, v11           /* ti2 */
    vfsub.vv v7, v1, v11           /* ti3 */
    vfadd.vv v8, v0, v10           /* tr2 */
    vfsub.vv v9, v0, v10           /* tr3 */
    vfadd.vv v14, v2, v8           /* outA */
    vfadd.vv v15, v4, v6
    vssseg2e32.v v14, (t4), t6
    vfadd.vv v14, v5, v9           /* outC */
    vfadd.vv v15, v3, v7
    slli    a2, t0, 1
    add     a2, a2, t4
    vssseg2e32.v v14, (a2), t6
    vfsub.vv v14, v9, v5           /* outB */
    vfsub.vv v15, v3, v7
    vssseg2e32.v v14, (t5), t6
    vfsub.vv v14, v8, v2           /* outD */
    vfsub.vv v15, v4, v6
    slli    a2, t0, 1
    add     a2, a2, t5
    vssseg2e32.v v14, (a2), t6
    mul     a2, t2, t0
    add     t3, t3, a2
    slli    a2, a2, 2
    add     t4, t4, a2
    add     t5, t5, a2
    sub     t1, t1, t2
    bnez    t1, .Lf4_k_strip
    addi    s1, s1, 8
    addi    s2, s2, 8
    addi    s3, s3, -8
    addi    s4, s4, 8
    addi    s5, s5, 8
    addi    a6, a6, 8
    addi    s0, s0, -1
    bnez    s0, .Lf4_k_pair

.Lf4_p2_done:
    ld      s0, 0(sp)
    ld      s1, 8(sp)
    ld      s2, 16(sp)
    ld      s3, 24(sp)
    ld      s4, 32(sp)
    ld      s5, 40(sp)
    ld      a2, 48(sp)
    ld      a3, 56(sp)
    addi    sp, sp, 64
    andi    t1, a0, 1
    bnez    t1, .Lf4_ret           /* odd ido: no part 3 */

/* ---- part 3: Nyquist column (even ido), lanes = rows ----
 * a = cc[t0+ido-1 + k*ido], b = a+2t0, c = a+t0, d = cc[ido-1 + k*ido]
 * ti1 = -hsqt2*(a+b); tr1 = hsqt2*(a-b)
 * ch[4k*ido + {ido-1, ido, 3ido-1, 3ido}] = {tr1+d, ti1-c, d-tr1, ti1+c} */
.Lf4_p3:
    li      t2, HSQT2_BITS
    li      t3, NHSQT2_BITS
    vsetvli t1, zero, e32, m4, ta, ma
    vmv.v.x v24, t2                /* hsqt2 */
    vmv.v.x v28, t3                /* -hsqt2 */
    mv      t1, a1
    add     t2, a2, a7
    add     t2, t2, t0
    addi    t2, t2, -4             /* &a */
    add     t3, a3, t0
    addi    t3, t3, -4             /* &f0 */
    slli    t6, t0, 2              /* ch lane stride 4*idoB */
1:
    vsetvli t4, t1, e32, m4, ta, ma
    slli    t5, a7, 1
    add     t5, t5, t2
    vlse32.v v0, (t2), t0          /* a */
    vlse32.v v4, (t5), t0          /* b */
    add     t5, t2, a7
    vlse32.v v8, (t5), t0          /* c */
    sub     t5, t2, a7
    vlse32.v v12, (t5), t0         /* d */
    vfadd.vv v16, v0, v4
    vfmul.vv v16, v16, v28         /* ti1 */
    vfsub.vv v20, v0, v4
    vfmul.vv v20, v20, v24         /* tr1 */
    vfadd.vv v0, v20, v12          /* f0 = tr1+d */
    vsse32.v v0, (t3), t6
    vfsub.vv v0, v16, v8           /* f1 = ti1-c */
    addi    t5, t3, 4
    vsse32.v v0, (t5), t6
    vfsub.vv v0, v12, v20          /* f2 = d-tr1 */
    slli    t5, t0, 1
    add     t5, t5, t3
    vsse32.v v0, (t5), t6
    vfadd.vv v0, v16, v8           /* f3 = ti1+c */
    addi    t5, t5, 4
    vsse32.v v0, (t5), t6
    mul     t5, t4, t0
    add     t2, t2, t5
    slli    t5, t5, 2
    add     t3, t3, t5
    sub     t1, t1, t4
    bnez    t1, 1b
.Lf4_ret:
    ret
    .size spx_drft_rvv_radf4_f32, . - spx_drft_rvv_radf4_f32

/* ================================================================== */
/* void spx_drft_rvv_radb4_f32(int ido, int l1, const float *cc,
 *                             float *ch, const float *wa1,
 *                             const float *wa2, const float *wa3)      */
/* ================================================================== */
    .balign 4
    .globl  spx_drft_rvv_radb4_f32
    .hidden spx_drft_rvv_radb4_f32
    .type   spx_drft_rvv_radb4_f32, @function
spx_drft_rvv_radb4_f32:
    slli    t0, a0, 2              /* idoB */
    mul     a7, t0, a1             /* t0B */

/* ---- part 1: lanes = rows ----
 * r0 = cc[4k*ido], r1 = cc[..+2ido-1], r2 = cc[..+2ido], r3 = cc[..+4ido-1]
 * tr3 = 2r1; tr4 = 2r2; tr1 = r0-r3; tr2 = r0+r3
 * ch[k*ido + m*t0] = {tr2+tr3, tr1-tr4, tr2-tr3, tr1+tr4} */
    li      t1, 1
    bne     a0, t1, .Lb4_p1_strided

    /* ido==1: r0..r3 consecutive -> vlseg4; unit-stride stores */
    mv      t1, a1
    mv      t2, a2
    mv      t3, a3
1:
    vsetvli t4, t1, e32, m2, ta, ma
    vlseg4e32.v v0, (t2)           /* r0 v0, r1 v2, r2 v4, r3 v6 */
    vfadd.vv v8, v2, v2            /* tr3 */
    vfadd.vv v10, v4, v4           /* tr4 */
    vfsub.vv v12, v0, v6           /* tr1 */
    vfadd.vv v14, v0, v6           /* tr2 */
    vfadd.vv v16, v14, v8          /* f0 */
    vfsub.vv v18, v12, v10         /* f1 */
    vfsub.vv v20, v14, v8          /* f2 */
    vfadd.vv v22, v12, v10         /* f3 */
    vse32.v v16, (t3)
    add     t5, t3, a7
    vse32.v v18, (t5)
    add     t5, t5, a7
    vse32.v v20, (t5)
    add     t5, t5, a7
    vse32.v v22, (t5)
    slli    t5, t4, 4
    add     t2, t2, t5
    slli    t5, t4, 2
    add     t3, t3, t5
    sub     t1, t1, t4
    bnez    t1, 1b
    ret

.Lb4_p1_strided:
    mv      t1, a1
    mv      t2, a2
    mv      t3, a3
    slli    t6, t0, 2              /* cc lane stride 4*idoB */
1:
    vsetvli t4, t1, e32, m4, ta, ma
    vlse32.v v0, (t2), t6          /* r0 */
    slli    t5, t0, 1
    add     t5, t5, t2
    addi    t5, t5, -4
    vlse32.v v4, (t5), t6          /* r1 @ +2ido-1 */
    addi    t5, t5, 4
    vlse32.v v8, (t5), t6          /* r2 @ +2ido */
    slli    t5, t0, 2
    add     t5, t5, t2
    addi    t5, t5, -4
    vlse32.v v12, (t5), t6         /* r3 @ +4ido-1 */
    vfadd.vv v16, v4, v4           /* tr3 */
    vfadd.vv v20, v8, v8           /* tr4 */
    vfsub.vv v24, v0, v12          /* tr1 */
    vfadd.vv v28, v0, v12          /* tr2 */
    vfadd.vv v0, v28, v16          /* f0 */
    vfsub.vv v4, v24, v20          /* f1 */
    vfsub.vv v8, v28, v16          /* f2 */
    vfadd.vv v12, v24, v20         /* f3 */
    vsse32.v v0, (t3), t0
    add     t5, t3, a7
    vsse32.v v4, (t5), t0
    add     t5, t5, a7
    vsse32.v v8, (t5), t0
    add     t5, t5, a7
    vsse32.v v12, (t5), t0
    mul     t5, t4, t0
    add     t3, t3, t5             /* ch += vl*idoB */
    slli    t5, t5, 2
    add     t2, t2, t5             /* cc += vl*4*idoB */
    sub     t1, t1, t4
    bnez    t1, 1b

    li      t1, 2
    beq     a0, t1, .Lb4_p3
    addi    t1, a0, -1
    srai    t1, t1, 1              /* npairs = (ido-1)>>1, both parities */

/* ---- part 2 ---- */
    addi    sp, sp, -64
    sd      s0, 0(sp)
    sd      s1, 8(sp)
    sd      s2, 16(sp)
    sd      s3, 24(sp)
    sd      s4, 32(sp)
    sd      s5, 40(sp)
    sd      a2, 48(sp)
    sd      a3, 56(sp)

    blt     t1, a1, .Lb4_p2_k

/* -- across pairs: outer row loop (m2) --
 * A asc = cc[4k*ido+1..], B asc = A+2ido, C desc from 4k*ido+2ido-3,
 * D desc from 4k*ido+4ido-3; outputs asc at ch[k*ido+1 + m*t0]. */
    addi    s0, a2, 4              /* A base */
    addi    s1, a3, 4              /* out0 base */
    slli    t2, t0, 1
    add     s2, a2, t2
    addi    s2, s2, -12            /* C base */
    mv      s3, a1                 /* row counter */
    slli    s4, t0, 1              /* 2idoB */
.Lb4_i_row:
    mv      t2, s0                 /* A walk */
    mv      t3, s1                 /* out0 walk */
    mv      t4, s2                 /* C walk (desc) */
    li      t5, 0                  /* twiddle offset */
    mv      t6, t1
.Lb4_i_strip:
    vsetvli s5, t6, e32, m2, ta, ma
    /* Descending streams: load unit-stride from the strip's lowest pair
     * address and reverse the lanes (strided loads run element-at-a-time
     * on the K1). v28 = reversal indices, a3 = (vl-1)*8. */
    vid.v   v28
    addi    a3, s5, -1
    vrsub.vx v28, v28, a3
    slli    a3, s5, 3
    addi    a3, a3, -8
    vlseg2e32.v v0, (t2)           /* Ar v0, Ai v2 */
    add     a2, t4, s4
    sub     a2, a2, a3
    vlseg2e32.v v16, (a2)          /* D, ascending raw */
    vrgather.vv v4, v16, v28       /* Dr */
    vrgather.vv v6, v18, v28       /* Di */
    vfadd.vv v8, v2, v6            /* ti1 = Ai+Di */
    vfsub.vv v10, v2, v6           /* ti2 = Ai-Di */
    vfsub.vv v12, v0, v4           /* tr1 = Ar-Dr */
    vfadd.vv v14, v0, v4           /* tr2 = Ar+Dr */
    add     a2, t2, s4
    vlseg2e32.v v0, (a2)           /* Br v0, Bi v2 */
    sub     a2, t4, a3
    vlseg2e32.v v16, (a2)          /* C, ascending raw */
    vrgather.vv v4, v16, v28       /* Cr */
    vrgather.vv v6, v18, v28       /* Ci */
    vfsub.vv v16, v2, v6           /* ti3 = Bi-Ci */
    vfadd.vv v18, v2, v6           /* tr4 = Bi+Ci */
    vfsub.vv v20, v0, v4           /* ti4 = Br-Cr */
    vfadd.vv v22, v0, v4           /* tr3 = Br+Cr */
    vfadd.vv v24, v14, v22         /* out0r = tr2+tr3 */
    vfsub.vv v14, v14, v22         /* cr3 */
    vfadd.vv v26, v10, v16         /* out0i = ti2+ti3 */
    vfsub.vv v10, v10, v16         /* ci3 */
    vsseg2e32.v v24, (t3)
    vfsub.vv v28, v12, v18         /* cr2 = tr1-tr4 */
    vfadd.vv v12, v12, v18         /* cr4 = tr1+tr4 */
    vfadd.vv v30, v8, v20          /* ci2 = ti1+ti4 */
    vfsub.vv v8, v8, v20           /* ci4 = ti1-ti4 */
    add     a2, a4, t5
    vlseg2e32.v v0, (a2)           /* w1 */
    drft_cmulb v4, v6, v0, v2, v28, v30
    add     a2, t3, a7
    vsseg2e32.v v4, (a2)
    add     a3, a5, t5
    vlseg2e32.v v0, (a3)           /* w2 */
    drft_cmulb v4, v6, v0, v2, v14, v10
    add     a2, a2, a7
    vsseg2e32.v v4, (a2)
    add     a3, a6, t5
    vlseg2e32.v v0, (a3)           /* w3 */
    drft_cmulb v4, v6, v0, v2, v12, v8
    add     a2, a2, a7
    vsseg2e32.v v4, (a2)
    slli    a2, s5, 3
    add     t2, t2, a2
    add     t3, t3, a2
    add     t5, t5, a2
    sub     t4, t4, a2
    sub     t6, t6, s5
    bnez    t6, .Lb4_i_strip
    slli    a2, t0, 2
    add     s0, s0, a2             /* cc row stride 4ido */
    add     s2, s2, a2
    add     s1, s1, t0             /* ch row stride ido */
    addi    s3, s3, -1
    bnez    s3, .Lb4_i_row
    j       .Lb4_p2_done

/* -- across rows: outer pair loop (m1, broadcast twiddles) --
 * cc lane stride 4*idoB, ch lane stride idoB. */
.Lb4_p2_k:
    mv      s0, t1
    addi    s1, a2, 4              /* A base (+8 per pair) */
    slli    t2, t0, 1
    add     s2, a2, t2
    addi    s2, s2, -12            /* C base (-8 per pair) */
    addi    s3, a3, 4              /* out0 base (+8 per pair) */
    mv      s4, a4
    mv      s5, a5                 /* w3 walks in a6 */
    slli    t6, t0, 2              /* cc lane stride */
.Lb4_k_pair:
    vsetvli t2, a1, e32, m1, ta, ma
    lw      a2, 0(s4)
    vmv.v.x v24, a2                /* w1r */
    lw      a2, 4(s4)
    vmv.v.x v25, a2                /* w1i */
    lw      a2, 0(s5)
    vmv.v.x v26, a2                /* w2r */
    lw      a2, 4(s5)
    vmv.v.x v27, a2                /* w2i */
    lw      a2, 0(a6)
    vmv.v.x v28, a2                /* w3r */
    lw      a2, 4(a6)
    vmv.v.x v29, a2                /* w3i */
    mv      t1, a1
    mv      t3, s1                 /* A walk */
    mv      t4, s2                 /* C walk */
    mv      t5, s3                 /* out0 walk */
.Lb4_k_strip:
    vsetvli t2, t1, e32, m1, ta, ma
    vlsseg2e32.v v0, (t3), t6      /* A */
    slli    a2, t0, 1
    add     a2, a2, t4
    vlsseg2e32.v v2, (a2), t6      /* D @ C+2ido */
    vfadd.vv v4, v1, v3            /* ti1 */
    vfsub.vv v5, v1, v3            /* ti2 */
    vfsub.vv v6, v0, v2            /* tr1 */
    vfadd.vv v7, v0, v2            /* tr2 */
    slli    a2, t0, 1
    add     a2, a2, t3
    vlsseg2e32.v v0, (a2), t6      /* B @ A+2ido */
    vlsseg2e32.v v2, (t4), t6      /* C */
    vfsub.vv v8, v1, v3            /* ti3 */
    vfadd.vv v9, v1, v3            /* tr4 */
    vfsub.vv v10, v0, v2           /* ti4 */
    vfadd.vv v11, v0, v2           /* tr3 */
    vfadd.vv v14, v7, v11          /* out0r */
    vfsub.vv v7, v7, v11           /* cr3 */
    vfadd.vv v15, v5, v8           /* out0i */
    vfsub.vv v5, v5, v8            /* ci3 */
    vssseg2e32.v v14, (t5), t0
    vfsub.vv v12, v6, v9           /* cr2 */
    vfadd.vv v13, v6, v9           /* cr4 */
    vfadd.vv v8, v4, v10           /* ci2 */
    vfsub.vv v9, v4, v10           /* ci4 */
    drft_cmulb v14, v15, v24, v25, v12, v8
    add     a2, t5, a7
    vssseg2e32.v v14, (a2), t0
    drft_cmulb v14, v15, v26, v27, v7, v5
    add     a2, a2, a7
    vssseg2e32.v v14, (a2), t0
    drft_cmulb v14, v15, v28, v29, v13, v9
    add     a2, a2, a7
    vssseg2e32.v v14, (a2), t0
    mul     a2, t2, t0
    add     t5, t5, a2             /* ch += vl*idoB */
    slli    a2, a2, 2
    add     t3, t3, a2             /* cc += vl*4*idoB */
    add     t4, t4, a2
    sub     t1, t1, t2
    bnez    t1, .Lb4_k_strip
    addi    s1, s1, 8
    addi    s2, s2, -8
    addi    s3, s3, 8
    addi    s4, s4, 8
    addi    s5, s5, 8
    addi    a6, a6, 8
    addi    s0, s0, -1
    bnez    s0, .Lb4_k_pair

.Lb4_p2_done:
    ld      s0, 0(sp)
    ld      s1, 8(sp)
    ld      s2, 16(sp)
    ld      s3, 24(sp)
    ld      s4, 32(sp)
    ld      s5, 40(sp)
    ld      a2, 48(sp)
    ld      a3, 56(sp)
    addi    sp, sp, 64
    andi    t1, a0, 1
    bnez    t1, .Lb4_ret

/* ---- part 3: lanes = rows ----
 * a = cc[ido + 4k*ido], b = cc[3ido + ..], c = cc[ido-1 + ..],
 * d = cc[3ido-1 + ..]
 * ti1 = a+b; ti2 = b-a; tr1 = c-d; tr2 = c+d
 * ch[ido-1 + k*ido + m*t0] = {2tr2, sqrt2*(tr1-ti1), 2ti2,
 *                             -sqrt2*(tr1+ti1)} */
.Lb4_p3:
    li      t2, SQRT2_BITS
    li      t3, NSQRT2_BITS
    vsetvli t1, zero, e32, m4, ta, ma
    vmv.v.x v24, t2                /* sqrt2 */
    vmv.v.x v28, t3                /* -sqrt2 */
    mv      t1, a1
    add     t2, a2, t0             /* &a */
    add     t3, a3, t0
    addi    t3, t3, -4             /* &f0 */
    slli    t6, t0, 2              /* cc lane stride 4*idoB */
1:
    vsetvli t4, t1, e32, m4, ta, ma
    vlse32.v v0, (t2), t6          /* a */
    slli    t5, t0, 1
    add     t5, t5, t2
    vlse32.v v4, (t5), t6          /* b @ a+2ido */
    addi    t5, t2, -4
    vlse32.v v8, (t5), t6          /* c */
    slli    t5, t0, 1
    add     t5, t5, t2
    addi    t5, t5, -4
    vlse32.v v12, (t5), t6         /* d */
    vfadd.vv v16, v0, v4           /* ti1 */
    vfsub.vv v20, v4, v0           /* ti2 */
    vfsub.vv v0, v8, v12           /* tr1 */
    vfadd.vv v8, v8, v12           /* tr2 */
    vfadd.vv v8, v8, v8            /* f0 = 2tr2 */
    vsse32.v v8, (t3), t0
    vfsub.vv v8, v0, v16
    vfmul.vv v8, v8, v24           /* f1 = sqrt2*(tr1-ti1) */
    add     t5, t3, a7
    vsse32.v v8, (t5), t0
    vfadd.vv v20, v20, v20         /* f2 = 2ti2 */
    add     t5, t5, a7
    vsse32.v v20, (t5), t0
    vfadd.vv v8, v0, v16
    vfmul.vv v8, v8, v28           /* f3 = -sqrt2*(tr1+ti1) */
    add     t5, t5, a7
    vsse32.v v8, (t5), t0
    mul     t5, t4, t0
    add     t3, t3, t5             /* ch += vl*idoB */
    slli    t5, t5, 2
    add     t2, t2, t5             /* cc += vl*4*idoB */
    sub     t1, t1, t4
    bnez    t1, 1b
.Lb4_ret:
    ret
    .size spx_drft_rvv_radb4_f32, . - spx_drft_rvv_radb4_f32

/* ================================================================== */
/* void spx_drft_rvv_radf2_f32(int ido, const float *cc, float *ch,
 *                             const float *wa1)
 * dradf2 restricted to l1==1 (the only shape drftf1 produces here).    */
/* ================================================================== */
    .balign 4
    .globl  spx_drft_rvv_radf2_f32
    .hidden spx_drft_rvv_radf2_f32
    .type   spx_drft_rvv_radf2_f32, @function
spx_drft_rvv_radf2_f32:
    slli    t0, a0, 2              /* idoB (== t0B since l1==1) */

    /* part 1 (vl=1): ch[0] = cc[0]+cc[ido]; ch[2ido-1] = cc[0]-cc[ido] */
    vsetivli zero, 1, e32, m1, ta, ma
    add     t1, a1, t0
    vle32.v v0, (a1)
    vle32.v v1, (t1)
    vfadd.vv v2, v0, v1
    vse32.v v2, (a2)
    vfsub.vv v2, v0, v1
    slli    t1, t0, 1
    add     t1, t1, a2
    addi    t1, t1, -4
    vse32.v v2, (t1)

    li      t1, 2
    blt     a0, t1, .Lf2_ret
    beq     a0, t1, .Lf2_p3

    /* part 2: pairs, strip-mined (m2)
     * A = cc[i-1..i], B = cc[ido+i-1..i] twiddled -> tr2/ti2
     * asc  ch[i-1..i]        = (Ar+tr2, Ai+ti2)
     * desc ch[2ido-i-1..-i]  = (Ar-tr2, ti2-Ai) */
    addi    t1, a0, -1
    srai    t1, t1, 1              /* npairs = (ido-1)>>1, both parities */
    addi    t2, a1, 4              /* A walk */
    add     t3, a1, t0
    addi    t3, t3, 4              /* B walk */
    mv      t4, a3                 /* wa1 walk */
    addi    t5, a2, 4              /* asc walk */
    slli    t6, t0, 1
    add     t6, t6, a2
    addi    t6, t6, -12            /* desc walk */
1:
    vsetvli a4, t1, e32, m2, ta, ma
    vlseg2e32.v v0, (t2)           /* Ar v0, Ai v2 */
    vlseg2e32.v v4, (t3)           /* Br v4, Bi v6 */
    vlseg2e32.v v8, (t4)           /* w1 */
    drft_cmulf v12, v14, v8, v10, v4, v6     /* tr2/ti2 */
    vfadd.vv v16, v0, v12          /* asc = (Ar+tr2, Ai+ti2) */
    vfadd.vv v18, v2, v14
    vsseg2e32.v v16, (t5)
    vid.v   v24
    addi    a5, a4, -1
    vrsub.vx v24, v24, a5          /* lane reversal indices */
    vfsub.vv v16, v0, v12          /* desc = (Ar-tr2, ti2-Ai) */
    vfsub.vv v18, v14, v2
    vrgather.vv v20, v16, v24
    vrgather.vv v22, v18, v24
    slli    a5, a4, 3
    addi    a5, a5, -8
    sub     a5, t6, a5             /* lowest pair of the strip */
    vsseg2e32.v v20, (a5)
    slli    a5, a4, 3
    add     t2, t2, a5
    add     t3, t3, a5
    add     t4, t4, a5
    add     t5, t5, a5
    sub     t6, t6, a5
    sub     t1, t1, a4
    bnez    t1, 1b
    andi    t1, a0, 1
    bnez    t1, .Lf2_ret

    /* part 3 (vl=1): ch[ido-1] = cc[ido-1]; ch[ido] = -cc[2ido-1] */
.Lf2_p3:
    vsetivli zero, 1, e32, m1, ta, ma
    add     t1, a1, t0
    addi    t1, t1, -4
    vle32.v v0, (t1)               /* cc[ido-1] */
    slli    t2, t0, 1
    add     t2, t2, a1
    addi    t2, t2, -4
    vle32.v v1, (t2)               /* cc[2ido-1] */
    add     t1, a2, t0
    addi    t1, t1, -4
    vse32.v v0, (t1)
    vfsgnjn.vv v1, v1, v1
    addi    t1, t1, 4
    vse32.v v1, (t1)
.Lf2_ret:
    ret
    .size spx_drft_rvv_radf2_f32, . - spx_drft_rvv_radf2_f32

/* ================================================================== */
/* void spx_drft_rvv_radb2_f32(int ido, const float *cc, float *ch,
 *                             const float *wa1)
 * dradb2 restricted to l1==1.                                          */
/* ================================================================== */
    .balign 4
    .globl  spx_drft_rvv_radb2_f32
    .hidden spx_drft_rvv_radb2_f32
    .type   spx_drft_rvv_radb2_f32, @function
spx_drft_rvv_radb2_f32:
    slli    t0, a0, 2              /* idoB (== t0B since l1==1) */

    /* part 1 (vl=1): ch[0] = cc[0]+cc[2ido-1]; ch[ido] = cc[0]-cc[2ido-1] */
    vsetivli zero, 1, e32, m1, ta, ma
    slli    t1, t0, 1
    add     t1, t1, a1
    addi    t1, t1, -4
    vle32.v v0, (a1)
    vle32.v v1, (t1)
    vfadd.vv v2, v0, v1
    vse32.v v2, (a2)
    vfsub.vv v2, v0, v1
    add     t1, a2, t0
    vse32.v v2, (t1)

    li      t1, 2
    blt     a0, t1, .Lb2_ret
    beq     a0, t1, .Lb2_p3

    /* part 2: pairs, strip-mined (m2)
     * A asc = cc[i-1..i]; B desc = cc[2ido-i-1..-i]
     * ch[i-1..i]     = (Ar+Br, Ai-Bi)
     * ch[ido+i-1..i] = w1 * (Ar-Br, Ai+Bi)  [backward twiddle] */
    addi    t1, a0, -1
    srai    t1, t1, 1              /* npairs = (ido-1)>>1, both parities */
    addi    t2, a1, 4              /* A walk */
    slli    t3, t0, 1
    add     t3, t3, a1
    addi    t3, t3, -12            /* B walk (desc) */
    mv      t4, a3                 /* wa1 walk */
    addi    t5, a2, 4              /* out0 walk */
    add     t6, a2, t0
    addi    t6, t6, 4              /* out1 walk */
1:
    vsetvli a4, t1, e32, m2, ta, ma
    vlseg2e32.v v0, (t2)           /* Ar v0, Ai v2 */
    vid.v   v20
    addi    a5, a4, -1
    vrsub.vx v20, v20, a5          /* lane reversal indices */
    slli    a5, a4, 3
    addi    a5, a5, -8
    sub     a5, t3, a5
    vlseg2e32.v v16, (a5)          /* B, ascending raw */
    vrgather.vv v4, v16, v20       /* Br */
    vrgather.vv v6, v18, v20       /* Bi */
    vfadd.vv v8, v0, v4            /* out0 = (Ar+Br, Ai-Bi) */
    vfsub.vv v10, v2, v6
    vsseg2e32.v v8, (t5)
    vfsub.vv v12, v0, v4           /* tr2 */
    vfadd.vv v14, v2, v6           /* ti2 */
    vlseg2e32.v v8, (t4)           /* w1 */
    drft_cmulb v16, v18, v8, v10, v12, v14
    vsseg2e32.v v16, (t6)
    slli    a5, a4, 3
    add     t2, t2, a5
    add     t4, t4, a5
    add     t5, t5, a5
    add     t6, t6, a5
    sub     t3, t3, a5
    sub     t1, t1, a4
    bnez    t1, 1b
    andi    t1, a0, 1
    bnez    t1, .Lb2_ret

    /* part 3 (vl=1): ch[ido-1] = 2*cc[ido-1]; ch[2ido-1] = -2*cc[ido] */
.Lb2_p3:
    vsetivli zero, 1, e32, m1, ta, ma
    add     t1, a1, t0
    addi    t1, t1, -4
    vle32.v v0, (t1)               /* cc[ido-1] */
    addi    t1, t1, 4
    vle32.v v1, (t1)               /* cc[ido] */
    vfadd.vv v0, v0, v0
    vfadd.vv v1, v1, v1
    vfsgnjn.vv v1, v1, v1
    add     t1, a2, t0
    addi    t1, t1, -4
    vse32.v v0, (t1)
    slli    t1, t0, 1
    add     t1, t1, a2
    addi    t1, t1, -4
    vse32.v v1, (t1)
.Lb2_ret:
    ret
    .size spx_drft_rvv_radb2_f32, . - spx_drft_rvv_radb2_f32

#if defined(__ELF__)
/* Non-executable stack note */
    .section .note.GNU-stack, "", @progbits
#endif
