/* Copyright (C) 2003-2004 Mark Borgerding
 * Copyright (C) 2005-2007 Jean-Marc Valin
 * Copyright (C) 2026 Tristan Matthews
 */
/*
   kiss_fft_rvv_asm.S -- hand-written RVV 1.0 KISS FFT butterflies.

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

/* RVV 1.0 kernels (VLEN >= 128) for kf_work's radix-2/3/4/5 butterfly
 * stages: N sub-FFTs (mm apart), each with m butterflies whose legs sit m
 * apart and whose leg-k twiddles stride the table by k*fstride. All work is
 * elementwise, so each kernel strip-mines one axis, tail-agnostic: the m
 * butterflies when m >= N, or -- radix 2/4 only, which receive N/mm -- the
 * N sub-FFTs when N > m, keeping the small-m tail stages of power-of-2
 * FFTs at vl=N. Leg k's twiddle pointer is derived per strip from leg 1's
 * as base + k*(tw1 - base).
 *
 * Fixed point matches C bit for bit: vwmul exact products, vssra with
 * vxrm=rnu (== sround/PSHR add-then-shift rounding), truncating vnsrl
 * narrows (== the wrapping int16 casts; vnclip would saturate), wrapping
 * e16 adds. SEW toggles between e16m1 and e32m2 use vl-preserving
 * `vsetvli zero, zero`. Float kernels use FMAs and reordered sums, so they
 * pair with a small checkasm tolerance; C_FIXDIV is a no-op in float, so
 * only radix-4 depends on the direction there (a store swap).
 *
 * `.option arch, +v` keeps the file self-contained (no -march needed).
 * Leaf functions, caller-saved registers only, integer args only (no
 * float-ABI dependency). kf_work guarantees m, N >= 1 and in-bounds
 * twiddle loads. */

    .option arch, +v
    .text

/* ------------------------------------------------------------------ */
/* unsigned int spx_kf_rvv_vlenb(void) -- VLEN in bytes. */
    .balign 4
    .globl  spx_kf_rvv_vlenb
    .hidden spx_kf_rvv_vlenb
    .type   spx_kf_rvv_vlenb, @function
spx_kf_rvv_vlenb:
    vsetvli a0, zero, e8, m1, ta, ma
    ret
    .size spx_kf_rvv_vlenb, . - spx_kf_rvv_vlenb

/* ------------------------------------------------------------------ */
/* int spx_kf_rvv_compliant(void)
 * Reject draft RVV 0.7.1 hardware (sets VILL on a ta,ma vtype config).
 * Mirrors spx_resample_rvv_compliant. */
    .balign 4
    .globl  spx_kf_rvv_compliant
    .hidden spx_kf_rvv_compliant
    .type   spx_kf_rvv_compliant, @function
spx_kf_rvv_compliant:
    vsetvli t0, zero, e8, m1, ta, ma
    csrr    a0, vtype
    sgtz    a0, a0
    ret
    .size spx_kf_rvv_compliant, . - spx_kf_rvv_compliant

/* ================================================================== */
/*                        float (f32) kernels                         */
/* ================================================================== */

/* ------------------------------------------------------------------ */
/* Float complex multiply, vector twiddle: (\dr,\di) = (\fr,\fi)*(\tr,\ti). */
.macro kf_f32_cmul dr, di, fr, fi, tr, ti
    vfmul.vv   \dr, \fr, \tr
    vfnmsac.vv \dr, \fi, \ti
    vfmul.vv   \di, \fr, \ti
    vfmacc.vv  \di, \fi, \tr
.endm

/* Float complex multiply, scalar twiddle in FP regs (\tr,\ti). */
.macro kf_f32_cmul_vf dr, di, fr, fi, tr, ti
    vfmul.vf   \dr, \fr, \tr
    vfnmsac.vf \dr, \ti, \fi
    vfmul.vf   \di, \fr, \ti
    vfmacc.vf  \di, \tr, \fi
.endm

/* void spx_kf_rvv_bfly2_f32(kiss_fft_cpx *Fout, const kiss_fft_cpx *tw,
 *                           size_t fstride, int m, int N, int mm, int inverse)
 * t = Fout[m]*tw; Fout[m] = Fout - t; Fout += t. Direction-independent in
 * float, so a6 is ignored. */
    .balign 4
    .globl  spx_kf_rvv_bfly2_f32
    .hidden spx_kf_rvv_bfly2_f32
    .type   spx_kf_rvv_bfly2_f32, @function
spx_kf_rvv_bfly2_f32:
    slli    t0, a3, 3               /* m*8 = byte offset of Fout[m] */
    slli    t1, a2, 3               /* fstride*8 = tw byte stride */
    bltu    a3, a4, .Lb2f_i         /* N > m: vectorize across sub-FFTs */

    slli    a2, a5, 3               /* mm*8: sub-FFT advance */
.Lb2f_j_outer:                      /* for each of the N sub-FFTs */
    mv      t2, a0                  /* Fout run */
    add     t3, a0, t0              /* Fout2 run */
    mv      t4, a1                  /* tw run */
    mv      t5, a3                  /* rem = m */
.Lb2f_j_strip:
    vsetvli t6, t5, e32, m1, ta, ma
    vlseg2e32.v  v0, (t2)           /* v0=F0.r v1=F0.i */
    vlseg2e32.v  v2, (t3)           /* v2=F2.r v3=F2.i */
    vlsseg2e32.v v4, (t4), t1       /* v4=tw.r v5=tw.i */
    kf_f32_cmul v6,v7, v2,v3, v4,v5    /* t = F2*tw */
    vfsub.vv   v2, v0, v6           /* Fout2 = Fout - t */
    vfsub.vv   v3, v1, v7
    vfadd.vv   v0, v0, v6           /* Fout  += t */
    vfadd.vv   v1, v1, v7
    vsseg2e32.v  v0, (t2)
    vsseg2e32.v  v2, (t3)
    slli    a7, t6, 3
    add     t2, t2, a7
    add     t3, t3, a7
    mul     a7, t6, t1
    add     t4, t4, a7
    sub     t5, t5, t6
    bnez    t5, .Lb2f_j_strip
    add     a0, a0, a2
    addi    a4, a4, -1
    bnez    a4, .Lb2f_j_outer
    ret

.Lb2f_i:                            /* per-j scalar twiddle, vector across i */
    slli    t2, a5, 3               /* mm*8: element stride */
    mv      t3, a0                  /* Fout + j */
    mv      t5, a1                  /* tw ptr */
.Lb2f_i_outer:                      /* for each of the m butterflies */
    flw     fa0, 0(t5)              /* tw.r */
    flw     fa1, 4(t5)              /* tw.i */
    mv      a6, t3                  /* Fout strip run */
    add     a5, t3, t0              /* Fout2 strip run */
    mv      t6, a4                  /* rem = N */
.Lb2f_i_strip:
    vsetvli a2, t6, e32, m1, ta, ma
    vlsseg2e32.v v0, (a6), t2
    vlsseg2e32.v v2, (a5), t2
    kf_f32_cmul_vf v6,v7, v2,v3, fa0,fa1
    vfsub.vv   v2, v0, v6
    vfsub.vv   v3, v1, v7
    vfadd.vv   v0, v0, v6
    vfadd.vv   v1, v1, v7
    vssseg2e32.v v0, (a6), t2
    vssseg2e32.v v2, (a5), t2
    mul     a7, a2, t2
    add     a6, a6, a7
    add     a5, a5, a7
    sub     t6, t6, a2
    bnez    t6, .Lb2f_i_strip
    addi    t3, t3, 8
    add     t5, t5, t1
    addi    a3, a3, -1
    bnez    a3, .Lb2f_i_outer
    ret
    .size spx_kf_rvv_bfly2_f32, . - spx_kf_rvv_bfly2_f32

/* ------------------------------------------------------------------ */
/* Radix-4 float body. In: F0=v0v1 F1=v2v3 F2=v4v5 F3=v6v7, T1..T3=v8..v13.
 * Out: Fout=v0v1, Fout[2m]=v4v5, rotated pair P=v2v3 / Q=v6v7 (forward
 * stores P to Fout[m], inverse swaps). */
.macro kf_b4f_body
    kf_f32_cmul v14,v15, v2,v3, v8,v9    /* s0 = F1*T1 */
    kf_f32_cmul v16,v17, v4,v5, v10,v11    /* s1 */
    kf_f32_cmul v18,v19, v6,v7, v12,v13    /* s2 */
    vfsub.vv   v20, v0, v16         /* s5 = F0 - s1 */
    vfsub.vv   v21, v1, v17
    vfadd.vv   v0, v0, v16          /* F0 += s1 */
    vfadd.vv   v1, v1, v17
    vfadd.vv   v22, v14, v18        /* s3 = s0 + s2 */
    vfadd.vv   v23, v15, v19
    vfsub.vv   v24, v14, v18        /* s4 = s0 - s2 */
    vfsub.vv   v25, v15, v19
    vfsub.vv   v4, v0, v22          /* Fout[2m] = F0 - s3 */
    vfsub.vv   v5, v1, v23
    vfadd.vv   v0, v0, v22          /* F0 += s3 */
    vfadd.vv   v1, v1, v23
    vfadd.vv   v2, v20, v25         /* P.r = s5.r + s4.i */
    vfsub.vv   v3, v21, v24         /* P.i = s5.i - s4.r */
    vfsub.vv   v6, v20, v25         /* Q.r = s5.r - s4.i */
    vfadd.vv   v7, v21, v24         /* Q.i = s5.i + s4.r */
.endm

/* void spx_kf_rvv_bfly4_f32(kiss_fft_cpx *Fout, const kiss_fft_cpx *tw,
 *                           size_t fstride, int m, int N, int mm, int inverse) */
    .balign 4
    .globl  spx_kf_rvv_bfly4_f32
    .hidden spx_kf_rvv_bfly4_f32
    .type   spx_kf_rvv_bfly4_f32, @function
spx_kf_rvv_bfly4_f32:
    slli    t0, a3, 3               /* m*8 */
    slli    t1, a2, 3               /* fstride*8 */
    bltu    a3, a4, .Lb4f_i

    slli    a2, a5, 3               /* mm*8 */
.Lb4f_j_outer:
    mv      t2, a0                  /* Fout run */
    mv      t3, a1                  /* tw1 run */
    mv      t4, a3                  /* rem = m */
.Lb4f_j_strip:
    vsetvli a5, t4, e32, m1, ta, ma
    vlseg2e32.v  v0, (t2)
    add     a7, t2, t0
    vlseg2e32.v  v2, (a7)
    add     a7, a7, t0
    vlseg2e32.v  v4, (a7)
    add     a7, a7, t0
    vlseg2e32.v  v6, (a7)
    vlsseg2e32.v v8, (t3), t1       /* T1, stride fstride */
    sub     t5, t3, a1              /* d = tw1 - base */
    add     t6, t3, t5              /* tw2 = base + 2d */
    slli    a7, t1, 1
    vlsseg2e32.v v10, (t6), a7      /* T2, stride 2*fstride */
    add     t6, t6, t5              /* tw3 = base + 3d */
    add     a7, a7, t1
    vlsseg2e32.v v12, (t6), a7      /* T3, stride 3*fstride */
    kf_b4f_body
    vsseg2e32.v  v0, (t2)
    add     a7, t2, t0
    beqz    a6, 1f
    vsseg2e32.v  v6, (a7)           /* inverse: Fout[m] = Q */
    add     a7, a7, t0
    vsseg2e32.v  v4, (a7)
    add     a7, a7, t0
    vsseg2e32.v  v2, (a7)           /* Fout[3m] = P */
    j       2f
1:
    vsseg2e32.v  v2, (a7)           /* forward: Fout[m] = P */
    add     a7, a7, t0
    vsseg2e32.v  v4, (a7)
    add     a7, a7, t0
    vsseg2e32.v  v6, (a7)           /* Fout[3m] = Q */
2:
    slli    a7, a5, 3
    add     t2, t2, a7
    mul     a7, a5, t1
    add     t3, t3, a7
    sub     t4, t4, a5
    bnez    t4, .Lb4f_j_strip
    add     a0, a0, a2
    addi    a4, a4, -1
    bnez    a4, .Lb4f_j_outer
    ret

.Lb4f_i:                            /* per-j scalar twiddles, vector across i */
    slli    a5, a5, 3               /* mm*8: element stride */
    mv      t3, a1                  /* tw1 ptr */
.Lb4f_i_outer:
    vsetvli a2, a4, e32, m1, ta, ma /* broadcast at the widest vl of this j */
    flw     fa0, 0(t3)
    flw     fa1, 4(t3)
    vfmv.v.f v8, fa0
    vfmv.v.f v9, fa1
    sub     t5, t3, a1              /* d */
    add     t6, t3, t5              /* tw2 */
    flw     fa0, 0(t6)
    flw     fa1, 4(t6)
    vfmv.v.f v10, fa0
    vfmv.v.f v11, fa1
    add     t6, t6, t5              /* tw3 */
    flw     fa0, 0(t6)
    flw     fa1, 4(t6)
    vfmv.v.f v12, fa0
    vfmv.v.f v13, fa1
    mv      t2, a0                  /* Fout + j strip run */
    mv      t4, a4                  /* rem = N */
.Lb4f_i_strip:
    vsetvli a2, t4, e32, m1, ta, ma
    vlsseg2e32.v v0, (t2), a5
    add     a7, t2, t0
    vlsseg2e32.v v2, (a7), a5
    add     a7, a7, t0
    vlsseg2e32.v v4, (a7), a5
    add     a7, a7, t0
    vlsseg2e32.v v6, (a7), a5
    kf_b4f_body
    vssseg2e32.v v0, (t2), a5
    add     a7, t2, t0
    beqz    a6, 1f
    vssseg2e32.v v6, (a7), a5
    add     a7, a7, t0
    vssseg2e32.v v4, (a7), a5
    add     a7, a7, t0
    vssseg2e32.v v2, (a7), a5
    j       2f
1:
    vssseg2e32.v v2, (a7), a5
    add     a7, a7, t0
    vssseg2e32.v v4, (a7), a5
    add     a7, a7, t0
    vssseg2e32.v v6, (a7), a5
2:
    mul     a7, a2, a5
    add     t2, t2, a7
    sub     t4, t4, a2
    bnez    t4, .Lb4f_i_strip
    addi    a0, a0, 8
    add     t3, t3, t1
    addi    a3, a3, -1
    bnez    a3, .Lb4f_i_outer
    ret
    .size spx_kf_rvv_bfly4_f32, . - spx_kf_rvv_bfly4_f32

/* ------------------------------------------------------------------ */
/* One radix-3 strip: legs at (t2 + k*t0), twiddles strided from t3, stores
   back in place, then advances t2/t3 by vl.  Register names are parameters
   so the same body serves the LMUL=1 and LMUL=2 instantiations. */
.macro kf_b3f_strip lmul, f0r,f0i, f1r,f1i, f2r,f2i, t1r,t1i, t2r,t2i, s1r,s1i, s2r,s2i, s3r,s3i
    vsetvli a2, t4, e32, \lmul, ta, ma
    vlseg2e32.v  \f0r, (t2)         /* F0 */
    add     a7, t2, t0
    vlseg2e32.v  \f1r, (a7)         /* F1 = Fout[m] */
    add     a7, a7, t0
    vlseg2e32.v  \f2r, (a7)         /* F2 = Fout[2m] */
    vlsseg2e32.v \t1r, (t3), t1     /* T1 */
    sub     t5, t3, a1
    add     t6, t3, t5              /* tw2 = base + 2d */
    slli    a7, t1, 1
    vlsseg2e32.v \t2r, (t6), a7     /* T2, stride 2*fstride */
    kf_f32_cmul \s1r,\s1i, \f1r,\f1i, \t1r,\t1i    /* s1 */
    kf_f32_cmul \s2r,\s2i, \f2r,\f2i, \t2r,\t2i    /* s2 */
    vfadd.vv   \s3r, \s1r, \s2r     /* s3 = s1 + s2 */
    vfadd.vv   \s3i, \s1i, \s2i
    vfsub.vv   \s1r, \s1r, \s2r     /* s0 = s1 - s2 */
    vfsub.vv   \s1i, \s1i, \s2i
    vfmul.vf   \f1r, \s3r, fa1
    vfsub.vv   \f1r, \f0r, \f1r     /* Fm.r = F0.r - 0.5*s3.r */
    vfmul.vf   \f1i, \s3i, fa1
    vfsub.vv   \f1i, \f0i, \f1i     /* Fm.i */
    vfmul.vf   \s1r, \s1r, fa0      /* s0 *= epi3.i */
    vfmul.vf   \s1i, \s1i, fa0
    vfadd.vv   \f0r, \f0r, \s3r     /* F0 += s3 */
    vfadd.vv   \f0i, \f0i, \s3i
    vfadd.vv   \f2r, \f1r, \s1i     /* Fout[2m] = (Fm.r + s0.i, Fm.i - s0.r) */
    vfsub.vv   \f2i, \f1i, \s1r
    vfsub.vv   \f1r, \f1r, \s1i     /* Fout[m]  = (Fm.r - s0.i, Fm.i + s0.r) */
    vfadd.vv   \f1i, \f1i, \s1r
    vsseg2e32.v  \f0r, (t2)
    add     a7, t2, t0
    vsseg2e32.v  \f1r, (a7)
    add     a7, a7, t0
    vsseg2e32.v  \f2r, (a7)
    slli    a7, a2, 3
    add     t2, t2, a7
    mul     a7, a2, t1
    add     t3, t3, a7
    sub     t4, t4, a2
.endm

/* void spx_kf_rvv_bfly3_f32(kiss_fft_cpx *Fout, const kiss_fft_cpx *tw,
 *                           size_t fstride, int m, int inverse)
 * epi3.i = tw[fstride*m].i; direction-independent in float (a4 ignored). */
    .balign 4
    .globl  spx_kf_rvv_bfly3_f32
    .hidden spx_kf_rvv_bfly3_f32
    .type   spx_kf_rvv_bfly3_f32, @function
spx_kf_rvv_bfly3_f32:
    mul     t2, a2, a3
    slli    t2, t2, 3
    add     t2, a1, t2
    flw     fa0, 4(t2)              /* epi3.i */
    li      a7, 0x3f000000
    fmv.w.x fa1, a7                 /* 0.5f (HALF_OF) */
    slli    t0, a3, 3               /* m*8 */
    slli    t1, a2, 3               /* fstride*8 */
    mv      t2, a0                  /* Fout run */
    mv      t3, a1                  /* tw1 run */
    mv      t4, a3                  /* rem = m */
    vsetvli a6, zero, e32, m1, ta, ma  /* a6 = VLMAX at m1 */
    bltu    a6, a3, .Lb3f_m2        /* m needs multiple m1 strips: LMUL=2 */
.Lb3f_strip:
    kf_b3f_strip m1, v0,v1, v2,v3, v4,v5, v6,v7, v8,v9, v10,v11, v12,v13, v14,v15
    bnez    t4, .Lb3f_strip
    ret
/* LMUL=2 body, taken when m*4 > vlenb so the across-u strip does not split
   into a full strip plus a short tail.  Same body, every vN renamed v2N. */
.Lb3f_m2:
    kf_b3f_strip m2, v0,v2, v4,v6, v8,v10, v12,v14, v16,v18, v20,v22, v24,v26, v28,v30
    bnez    t4, .Lb3f_m2
    ret
    .size spx_kf_rvv_bfly3_f32, . - spx_kf_rvv_bfly3_f32

/* ------------------------------------------------------------------ */
/* Segment store at \addr: unit-stride unless \stride is given. */
.macro kf_b5f_st reg, addr, stride
.ifb \stride
    vsseg2e32.v  \reg, (\addr)
.else
    vssseg2e32.v \reg, (\addr), \stride
.endif
.endm

/* Radix-5 combine: (F0, s1..s4) -> the five outputs at \base + k*t0.
   s9 lands in the s2 pair and s10 in the s1 pair; s7/s8 and the tA/tB/tC
   scratch pairs are free registers of the caller's choosing (tC also
   stages the stores).  ya/yb stay in fa0-fa3; clobbers a7.  Passing
   \stride makes the stores element-strided (the across-i path). */
.macro kf_b5f_combine f0r,f0i, s1r,s1i, s2r,s2i, s3r,s3i, s4r,s4i, s7r,s7i, s8r,s8i, tAr,tAi, tBr,tBi, tCr,tCi, base, stride
    vfadd.vv   \s8r, \s2r, \s3r     /* s8 = s2 + s3 */
    vfadd.vv   \s8i, \s2i, \s3i
    vfsub.vv   \s2r, \s2r, \s3r     /* s9 = s2 - s3 */
    vfsub.vv   \s2i, \s2i, \s3i
    vfadd.vv   \s7r, \s1r, \s4r     /* s7 = s1 + s4 */
    vfadd.vv   \s7i, \s1i, \s4i
    vfsub.vv   \s1r, \s1r, \s4r     /* s10 = s1 - s4 */
    vfsub.vv   \s1i, \s1i, \s4i
    vfadd.vv   \tAr, \f0r, \s7r     /* Fout = F0 + s7 + s8 */
    vfadd.vv   \tAr, \tAr, \s8r
    vfadd.vv   \tAi, \f0i, \s7i
    vfadd.vv   \tAi, \tAi, \s8i
    kf_b5f_st  \tAr, \base, \stride
    vmv.v.v    \tAr, \f0r           /* s5.r = F0.r + ya.r*s7.r + yb.r*s8.r */
    vfmacc.vf  \tAr, fa0, \s7r
    vfmacc.vf  \tAr, fa2, \s8r
    vmv.v.v    \tAi, \f0i
    vfmacc.vf  \tAi, fa0, \s7i
    vfmacc.vf  \tAi, fa2, \s8i      /* s5.i */
    vfmul.vf   \tBr, \s1i, fa1      /* s6.r = ya.i*s10.i + yb.i*s9.i */
    vfmacc.vf  \tBr, fa3, \s2i
    vfmul.vf   \tBi, \s1r, fa1      /* s6.i = -(ya.i*s10.r + yb.i*s9.r) */
    vfmacc.vf  \tBi, fa3, \s2r
    vfneg.v    \tBi, \tBi
    vfsub.vv   \tCr, \tAr, \tBr     /* Fout[m] = s5 - s6 */
    vfsub.vv   \tCi, \tAi, \tBi
    add     a7, \base, t0
    kf_b5f_st  \tCr, a7, \stride
    vfadd.vv   \tCr, \tAr, \tBr     /* Fout[4m] = s5 + s6 */
    vfadd.vv   \tCi, \tAi, \tBi
    slli    a7, t0, 2
    add     a7, \base, a7
    kf_b5f_st  \tCr, a7, \stride
    vmv.v.v    \tAr, \f0r           /* s11.r = F0.r + yb.r*s7.r + ya.r*s8.r */
    vfmacc.vf  \tAr, fa2, \s7r
    vfmacc.vf  \tAr, fa0, \s8r
    vmv.v.v    \tAi, \f0i
    vfmacc.vf  \tAi, fa2, \s7i
    vfmacc.vf  \tAi, fa0, \s8i      /* s11.i */
    vfmul.vf   \tBr, \s2i, fa1      /* s12.r = ya.i*s9.i - yb.i*s10.i */
    vfnmsac.vf \tBr, fa3, \s1i
    vfmul.vf   \tBi, \s1r, fa3      /* s12.i = yb.i*s10.r - ya.i*s9.r */
    vfnmsac.vf \tBi, fa1, \s2r
    vfadd.vv   \tCr, \tAr, \tBr     /* Fout[2m] = s11 + s12 */
    vfadd.vv   \tCi, \tAi, \tBi
    slli    a7, t0, 1
    add     a7, \base, a7
    kf_b5f_st  \tCr, a7, \stride
    vfsub.vv   \tCr, \tAr, \tBr     /* Fout[3m] = s11 - s12 */
    vfsub.vv   \tCi, \tAi, \tBi
    add     a7, a7, t0
    kf_b5f_st  \tCr, a7, \stride
.endm

/* void spx_kf_rvv_bfly5_f32(kiss_fft_cpx *Fout, const kiss_fft_cpx *tw,
 *                           size_t fstride, int m, int N, int mm, int inverse)
 * ya = tw[fstride*m], yb = tw[2*fstride*m]; direction-independent in float
 * (a6 ignored, reused as scratch).  N > m vectorizes across the N sub-FFTs
 * with scalar twiddles (this is what makes m==1..4 stages worthwhile);
 * otherwise across u within each sub-FFT, at LMUL=2 when one m1 register
 * cannot hold m lanes. */
    .balign 4
    .globl  spx_kf_rvv_bfly5_f32
    .hidden spx_kf_rvv_bfly5_f32
    .type   spx_kf_rvv_bfly5_f32, @function
spx_kf_rvv_bfly5_f32:
    mul     t2, a2, a3
    slli    t2, t2, 3
    add     t3, a1, t2
    flw     fa0, 0(t3)              /* ya.r */
    flw     fa1, 4(t3)              /* ya.i */
    add     t3, t3, t2
    flw     fa2, 0(t3)              /* yb.r */
    flw     fa3, 4(t3)              /* yb.i */
    slli    t0, a3, 3               /* m*8 */
    slli    t1, a2, 3               /* fstride*8 */
    bltu    a3, a4, .Lb5f_i         /* N > m: vectorize across sub-FFTs */
    slli    a6, a5, 3               /* mm*8: sub-FFT advance */
    vsetvli a5, zero, e32, m1, ta, ma  /* a5 = VLMAX at m1 */
    bltu    a5, a3, .Lb5f_j2_outer  /* m needs multiple m1 strips: LMUL=2 */
.Lb5f_j_outer:                      /* for each of the N sub-FFTs */
    mv      t2, a0                  /* Fout run */
    mv      t3, a1                  /* tw1 run */
    mv      t4, a3                  /* rem = m */
.Lb5f_strip:
    vsetvli a2, t4, e32, m1, ta, ma
    vlseg2e32.v  v0, (t2)           /* F0 (= scratch[0], kept in v0v1) */
    add     a7, t2, t0
    vlseg2e32.v  v2, (a7)
    add     a7, a7, t0
    vlseg2e32.v  v4, (a7)
    add     a7, a7, t0
    vlseg2e32.v  v6, (a7)
    add     a7, a7, t0
    vlseg2e32.v  v8, (a7)
    vlsseg2e32.v v10, (t3), t1      /* T1 */
    sub     t5, t3, a1              /* d */
    add     t6, t3, t5              /* tw2 */
    slli    a7, t1, 1
    vlsseg2e32.v v12, (t6), a7      /* T2, stride 2*fstride */
    add     t6, t6, t5              /* tw3 */
    add     a7, a7, t1
    vlsseg2e32.v v14, (t6), a7      /* T3, stride 3*fstride */
    add     t6, t6, t5              /* tw4 */
    add     a7, a7, t1
    vlsseg2e32.v v16, (t6), a7      /* T4, stride 4*fstride */
    kf_f32_cmul v18,v19, v2,v3, v10,v11    /* s1 */
    kf_f32_cmul v20,v21, v4,v5, v12,v13    /* s2 */
    kf_f32_cmul v22,v23, v6,v7, v14,v15    /* s3 */
    kf_f32_cmul v24,v25, v8,v9, v16,v17    /* s4 */
    kf_b5f_combine v0,v1, v18,v19, v20,v21, v22,v23, v24,v25, v2,v3, v6,v7, v26,v27, v30,v31, v28,v29, t2
    slli    a7, a2, 3
    add     t2, t2, a7
    mul     a7, a2, t1
    add     t3, t3, a7
    sub     t4, t4, a2
    bnez    t4, .Lb5f_strip
    add     a0, a0, a6
    addi    a4, a4, -1
    bnez    a4, .Lb5f_j_outer
    ret
/* LMUL=2 body, taken when m*4 > vlenb so the m1 strip would split into a
   full strip plus a short tail (e.g. m==5 at VLEN=128).  Only 16 register
   groups exist at m2, so each leg's twiddle product is folded right after
   its loads to retire Fk/Tk instead of keeping all loads live at once. */
.Lb5f_j2_outer:                     /* for each of the N sub-FFTs */
    mv      t2, a0                  /* Fout run */
    mv      t3, a1                  /* tw1 run */
    mv      t4, a3                  /* rem = m */
.Lb5f_m2:
    vsetvli a2, t4, e32, m2, ta, ma
    vlseg2e32.v  v0, (t2)           /* F0.r v0v1, F0.i v2v3 (kept) */
    add     a7, t2, t0
    vlseg2e32.v  v20, (a7)          /* F1 */
    vlsseg2e32.v v24, (t3), t1      /* T1 */
    kf_f32_cmul v4,v6, v20,v22, v24,v26    /* s1 */
    add     a7, a7, t0
    vlseg2e32.v  v20, (a7)          /* F2 */
    sub     t5, t3, a1              /* d */
    add     t6, t3, t5              /* tw2 */
    slli    a5, t1, 1
    vlsseg2e32.v v24, (t6), a5      /* T2, stride 2*fstride */
    kf_f32_cmul v8,v10, v20,v22, v24,v26    /* s2 */
    add     a7, a7, t0
    vlseg2e32.v  v20, (a7)          /* F3 */
    add     t6, t6, t5              /* tw3 */
    add     a5, a5, t1
    vlsseg2e32.v v24, (t6), a5      /* T3, stride 3*fstride */
    kf_f32_cmul v12,v14, v20,v22, v24,v26    /* s3 */
    add     a7, a7, t0
    vlseg2e32.v  v20, (a7)          /* F4 */
    add     t6, t6, t5              /* tw4 */
    add     a5, a5, t1
    vlsseg2e32.v v24, (t6), a5      /* T4, stride 4*fstride */
    kf_f32_cmul v16,v18, v20,v22, v24,v26    /* s4 */
    kf_b5f_combine v0,v2, v4,v6, v8,v10, v12,v14, v16,v18, v20,v22, v24,v26, v12,v14, v16,v18, v28,v30, t2
    slli    a7, a2, 3
    add     t2, t2, a7
    mul     a7, a2, t1
    add     t3, t3, a7
    sub     t4, t4, a2
    bnez    t4, .Lb5f_m2
    add     a0, a0, a6
    addi    a4, a4, -1
    bnez    a4, .Lb5f_j2_outer
    ret
/* Across-i path: for each of the m butterfly indices u, process the N
   sub-FFTs (mm apart) as vectors; the four leg twiddles depend only on u,
   so they stay in FP scalars for the whole strip loop.  This is the only
   path that pays off at m==1..4 (the across-u paths degenerate there). */
.Lb5f_i:
    slli    t2, a5, 3               /* mm*8: element stride */
    mv      t3, a1                  /* tw1 ptr */
.Lb5f_i_outer:                      /* for each of the m butterflies */
    flw     ft0, 0(t3)              /* T1.r */
    flw     ft1, 4(t3)              /* T1.i */
    sub     t5, t3, a1              /* d */
    add     t6, t3, t5              /* tw2 */
    flw     ft2, 0(t6)
    flw     ft3, 4(t6)
    add     t6, t6, t5              /* tw3 */
    flw     ft4, 0(t6)
    flw     ft5, 4(t6)
    add     t6, t6, t5              /* tw4 */
    flw     ft6, 0(t6)
    flw     ft7, 4(t6)
    mv      a6, a0                  /* F0 strip run */
    mv      t4, a4                  /* rem = N */
.Lb5f_i_strip:
    vsetvli a2, t4, e32, m1, ta, ma
    vlsseg2e32.v v0, (a6), t2       /* F0 */
    add     a7, a6, t0
    vlsseg2e32.v v2, (a7), t2       /* F1 */
    add     a7, a7, t0
    vlsseg2e32.v v4, (a7), t2       /* F2 */
    add     a7, a7, t0
    vlsseg2e32.v v6, (a7), t2       /* F3 */
    add     a7, a7, t0
    vlsseg2e32.v v8, (a7), t2       /* F4 */
    kf_f32_cmul_vf v10,v11, v2,v3, ft0,ft1    /* s1 */
    kf_f32_cmul_vf v12,v13, v4,v5, ft2,ft3    /* s2 */
    kf_f32_cmul_vf v14,v15, v6,v7, ft4,ft5    /* s3 */
    kf_f32_cmul_vf v16,v17, v8,v9, ft6,ft7    /* s4 */
    kf_b5f_combine v0,v1, v10,v11, v12,v13, v14,v15, v16,v17, v14,v15, v18,v19, v20,v21, v22,v23, v2,v3, a6, t2
    mul     a7, a2, t2
    add     a6, a6, a7
    sub     t4, t4, a2
    bnez    t4, .Lb5f_i_strip
    addi    a0, a0, 8
    add     t3, t3, t1
    addi    a3, a3, -1
    bnez    a3, .Lb5f_i_outer
    ret
    .size spx_kf_rvv_bfly5_f32, . - spx_kf_rvv_bfly5_f32

/* ================================================================== */
/*                      fixed-point (i16) kernels                     */
/* ================================================================== */

/* Radix-2 inverse body (C_MUL path). In: F0=v0v1 F2=v2v3 T=v4v5.
 * Out: Fout=v0v1, Fout[m]=v2v3. */
.macro kf_b2x_inv_body
    vwmul.vv v8, v2, v4             /* F2.r*tw.r */
    vwmul.vv v10, v3, v5            /* F2.i*tw.i */
    vwmul.vv v12, v2, v5            /* F2.r*tw.i */
    vwmul.vv v14, v3, v4            /* F2.i*tw.r */
    vsetvli zero, zero, e32, m2, ta, ma
    vsub.vv  v8, v8, v10
    vadd.vv  v12, v12, v14
    vssra.vi v8, v8, 15             /* sround */
    vssra.vi v12, v12, 15
    vsetvli zero, zero, e16, m1, ta, ma
    vnsrl.wi v6, v8, 0              /* t.r (truncating cast) */
    vnsrl.wi v7, v12, 0             /* t.i */
    vsub.vv  v2, v0, v6             /* Fout2 = Fout - t */
    vsub.vv  v3, v1, v7
    vadd.vv  v0, v0, v6             /* Fout += t */
    vadd.vv  v1, v1, v7
.endm

/* Radix-2 forward body (half-input scaling; a6 = 16384): t >>= 1, then
 * Fout2 = PSHR32((Fout<<14) - t, 15), Fout = PSHR32((Fout<<14) + t, 15). */
.macro kf_b2x_fwd_body
    vwmul.vv v8, v2, v4
    vwmul.vv v10, v3, v5
    vwmul.vv v12, v2, v5
    vwmul.vv v14, v3, v4
    vwmul.vx v16, v0, a6            /* Fout.r << 14 (exact: * 16384) */
    vwmul.vx v18, v1, a6            /* Fout.i << 14 */
    vsetvli zero, zero, e32, m2, ta, ma
    vsub.vv  v8, v8, v10
    vadd.vv  v12, v12, v14
    vsra.vi  v8, v8, 1              /* SHR32(., 1): no rounding */
    vsra.vi  v12, v12, 1
    vsub.vv  v10, v16, v8
    vadd.vv  v16, v16, v8
    vsub.vv  v14, v18, v12
    vadd.vv  v18, v18, v12
    vssra.vi v10, v10, 15           /* PSHR32(., 15) */
    vssra.vi v16, v16, 15
    vssra.vi v14, v14, 15
    vssra.vi v18, v18, 15
    vsetvli zero, zero, e16, m1, ta, ma
    vnsrl.wi v2, v10, 0             /* Fout2.r */
    vnsrl.wi v3, v14, 0
    vnsrl.wi v0, v16, 0             /* Fout.r */
    vnsrl.wi v1, v18, 0
.endm

/* Radix-2 across-j loop (m >= N): unit segment loads, strided twiddles. */
.macro kf_b2x_pathj body
    slli    a2, a5, 2               /* mm*4 */
1:
    mv      t2, a0
    add     t3, a0, t0
    mv      t4, a1
    mv      t5, a3
2:
    vsetvli t6, t5, e16, m1, ta, ma
    vlseg2e16.v  v0, (t2)
    vlseg2e16.v  v2, (t3)
    vlsseg2e16.v v4, (t4), t1
    \body
    vsseg2e16.v  v0, (t2)
    vsseg2e16.v  v2, (t3)
    slli    a7, t6, 2
    add     t2, t2, a7
    add     t3, t3, a7
    mul     a7, t6, t1
    add     t4, t4, a7
    sub     t5, t5, t6
    bnez    t5, 2b
    add     a0, a0, a2
    addi    a4, a4, -1
    bnez    a4, 1b
    ret
.endm

/* Radix-2 across-i loop (N > m): strided segment loads, broadcast twiddle. */
.macro kf_b2x_pathi body
    slli    t2, a5, 2               /* mm*4: element stride */
    mv      t3, a0                  /* Fout + j */
    mv      t5, a1                  /* tw ptr */
1:
    lh      a7, 0(t5)               /* tw.r */
    lh      a2, 2(t5)               /* tw.i */
    vsetvli a5, a4, e16, m1, ta, ma
    vmv.v.x v4, a7
    vmv.v.x v5, a2
    mv      t4, t3                  /* Fout strip run */
    add     a5, t3, t0              /* Fout2 strip run */
    mv      t6, a4                  /* rem = N */
2:
    vsetvli a2, t6, e16, m1, ta, ma
    vlsseg2e16.v v0, (t4), t2
    vlsseg2e16.v v2, (a5), t2
    \body
    vssseg2e16.v v0, (t4), t2
    vssseg2e16.v v2, (a5), t2
    mul     a7, a2, t2
    add     t4, t4, a7
    add     a5, a5, a7
    sub     t6, t6, a2
    bnez    t6, 2b
    addi    t3, t3, 4
    add     t5, t5, t1
    addi    a3, a3, -1
    bnez    a3, 1b
    ret
.endm

/* void spx_kf_rvv_bfly2_i16(kiss_fft_cpx *Fout, const kiss_fft_cpx *tw,
 *                           size_t fstride, int m, int N, int mm, int inverse) */
    .balign 4
    .globl  spx_kf_rvv_bfly2_i16
    .hidden spx_kf_rvv_bfly2_i16
    .type   spx_kf_rvv_bfly2_i16, @function
spx_kf_rvv_bfly2_i16:
    csrwi   vxrm, 0                 /* rnu: vssra == sround/PSHR rounding */
    slli    t0, a3, 2               /* m*4 */
    slli    t1, a2, 2               /* fstride*4 */
    bnez    a6, .Lb2x_inv
    li      a6, 16384               /* forward bodies scale Fout by <<14 */
    bltu    a3, a4, .Lb2x_fwd_i
    kf_b2x_pathj kf_b2x_fwd_body
.Lb2x_fwd_i:
    kf_b2x_pathi kf_b2x_fwd_body
.Lb2x_inv:
    bltu    a3, a4, .Lb2x_inv_i
    kf_b2x_pathj kf_b2x_inv_body
.Lb2x_inv_i:
    kf_b2x_pathi kf_b2x_inv_body
    .size spx_kf_rvv_bfly2_i16, . - spx_kf_rvv_bfly2_i16

/* ------------------------------------------------------------------ */
/* Radix-4 fixed-point C_MUL: dest = sround_a6(F*T); a6 = 17 fwd (C_MUL4),
 * 15 inv (C_MUL). */
.macro kf_b4x_cmul fr, fi, tr, ti, dr, di
    vwmul.vv v16, \fr, \tr
    vwmul.vv v18, \fi, \ti
    vwmul.vv v20, \fr, \ti
    vwmul.vv v22, \fi, \tr
    vsetvli zero, zero, e32, m2, ta, ma
    vsub.vv  v16, v16, v18
    vadd.vv  v20, v20, v22
    vssra.vx v16, v16, a6
    vssra.vx v20, v20, a6
    vsetvli zero, zero, e16, m1, ta, ma
    vnsrl.wi \dr, v16, 0
    vnsrl.wi \di, v20, 0
.endm

/* Radix-4 fixed-point body, both directions; forward also prescales Fout
 * by PSHR16(.,2). Register roles as in kf_b4f_body. (The C forward path's
 * PSHR16 of Fout[2m] is dead -- overwritten by C_SUB -- and skipped.) */
.macro kf_b4x_body
    kf_b4x_cmul v2, v3, v8, v9, v24, v25      /* s0 */
    kf_b4x_cmul v4, v5, v10, v11, v26, v27    /* s1 */
    kf_b4x_cmul v6, v7, v12, v13, v28, v29    /* s2 */
    xori    a7, a6, 17
    bnez    a7, 3f                  /* inverse: no input scaling */
    vssra.vi v0, v0, 2              /* Fout = PSHR16(Fout, 2) */
    vssra.vi v1, v1, 2
3:
    vsub.vv  v14, v0, v26           /* s5 = F0 - s1 */
    vsub.vv  v15, v1, v27
    vadd.vv  v0, v0, v26            /* F0 += s1 */
    vadd.vv  v1, v1, v27
    vadd.vv  v16, v24, v28          /* s3 = s0 + s2 */
    vadd.vv  v17, v25, v29
    vsub.vv  v18, v24, v28          /* s4 = s0 - s2 */
    vsub.vv  v19, v25, v29
    vsub.vv  v4, v0, v16            /* Fout[2m] = F0 - s3 */
    vsub.vv  v5, v1, v17
    vadd.vv  v0, v0, v16            /* F0 += s3 */
    vadd.vv  v1, v1, v17
    vadd.vv  v2, v14, v19           /* P.r = s5.r + s4.i */
    vsub.vv  v3, v15, v18           /* P.i = s5.i - s4.r */
    vsub.vv  v6, v14, v19           /* Q.r = s5.r - s4.i */
    vadd.vv  v7, v15, v18           /* Q.i = s5.i + s4.r */
.endm

/* void spx_kf_rvv_bfly4_i16(kiss_fft_cpx *Fout, const kiss_fft_cpx *tw,
 *                           size_t fstride, int m, int N, int mm, int inverse) */
    .balign 4
    .globl  spx_kf_rvv_bfly4_i16
    .hidden spx_kf_rvv_bfly4_i16
    .type   spx_kf_rvv_bfly4_i16, @function
spx_kf_rvv_bfly4_i16:
    csrwi   vxrm, 0
    slli    t0, a3, 2               /* m*4 */
    slli    t1, a2, 2               /* fstride*4 */
    beqz    a6, 1f
    li      a6, 15                  /* inverse: C_MUL shift */
    j       2f
1:
    li      a6, 17                  /* forward: C_MUL4 shift */
2:
    bltu    a3, a4, .Lb4x_i

    slli    a2, a5, 2               /* mm*4 */
.Lb4x_j_outer:
    mv      t2, a0
    mv      t3, a1
    mv      t4, a3
.Lb4x_j_strip:
    vsetvli a5, t4, e16, m1, ta, ma
    vlseg2e16.v  v0, (t2)
    add     a7, t2, t0
    vlseg2e16.v  v2, (a7)
    add     a7, a7, t0
    vlseg2e16.v  v4, (a7)
    add     a7, a7, t0
    vlseg2e16.v  v6, (a7)
    vlsseg2e16.v v8, (t3), t1
    sub     t5, t3, a1              /* d */
    add     t6, t3, t5              /* tw2 */
    slli    a7, t1, 1
    vlsseg2e16.v v10, (t6), a7
    add     t6, t6, t5              /* tw3 */
    add     a7, a7, t1
    vlsseg2e16.v v12, (t6), a7
    kf_b4x_body
    vsseg2e16.v  v0, (t2)
    add     a7, t2, t0
    xori    t6, a6, 17
    beqz    t6, 1f
    vsseg2e16.v  v6, (a7)           /* inverse: Fout[m] = Q */
    add     a7, a7, t0
    vsseg2e16.v  v4, (a7)
    add     a7, a7, t0
    vsseg2e16.v  v2, (a7)
    j       2f
1:
    vsseg2e16.v  v2, (a7)           /* forward: Fout[m] = P */
    add     a7, a7, t0
    vsseg2e16.v  v4, (a7)
    add     a7, a7, t0
    vsseg2e16.v  v6, (a7)
2:
    slli    a7, a5, 2
    add     t2, t2, a7
    mul     a7, a5, t1
    add     t3, t3, a7
    sub     t4, t4, a5
    bnez    t4, .Lb4x_j_strip
    add     a0, a0, a2
    addi    a4, a4, -1
    bnez    a4, .Lb4x_j_outer
    ret

.Lb4x_i:                            /* per-j scalar twiddles, vector across i */
    slli    a2, a5, 2               /* mm*4: element stride */
    mv      t3, a1                  /* tw1 ptr */
.Lb4x_i_outer:
    vsetvli a5, a4, e16, m1, ta, ma /* broadcast at the widest vl of this j */
    lh      a7, 0(t3)
    vmv.v.x v8, a7
    lh      a7, 2(t3)
    vmv.v.x v9, a7
    sub     t5, t3, a1              /* d */
    add     t6, t3, t5              /* tw2 */
    lh      a7, 0(t6)
    vmv.v.x v10, a7
    lh      a7, 2(t6)
    vmv.v.x v11, a7
    add     t6, t6, t5              /* tw3 */
    lh      a7, 0(t6)
    vmv.v.x v12, a7
    lh      a7, 2(t6)
    vmv.v.x v13, a7
    mv      t2, a0                  /* Fout + j strip run */
    mv      t4, a4                  /* rem = N */
.Lb4x_i_strip:
    vsetvli a5, t4, e16, m1, ta, ma
    vlsseg2e16.v v0, (t2), a2
    add     a7, t2, t0
    vlsseg2e16.v v2, (a7), a2
    add     a7, a7, t0
    vlsseg2e16.v v4, (a7), a2
    add     a7, a7, t0
    vlsseg2e16.v v6, (a7), a2
    kf_b4x_body
    vssseg2e16.v v0, (t2), a2
    add     a7, t2, t0
    xori    t6, a6, 17
    beqz    t6, 1f
    vssseg2e16.v v6, (a7), a2
    add     a7, a7, t0
    vssseg2e16.v v4, (a7), a2
    add     a7, a7, t0
    vssseg2e16.v v2, (a7), a2
    j       2f
1:
    vssseg2e16.v v2, (a7), a2
    add     a7, a7, t0
    vssseg2e16.v v4, (a7), a2
    add     a7, a7, t0
    vssseg2e16.v v6, (a7), a2
2:
    mul     a7, a5, a2
    add     t2, t2, a7
    sub     t4, t4, a5
    bnez    t4, .Lb4x_i_strip
    addi    a0, a0, 4
    add     t3, t3, t1
    addi    a3, a3, -1
    bnez    a3, .Lb4x_i_outer
    ret
    .size spx_kf_rvv_bfly4_i16, . - spx_kf_rvv_bfly4_i16

/* ------------------------------------------------------------------ */
/* LMUL=2 fixed-point C_MUL half: \d = sround(\ar*\br <op> \ai*\bi), run
   through the two e32,m4 wide groups v24-27 / v28-31 (\op is vsub for the
   real part, vadd for the imaginary part). */
.macro kf_i16_cmulh_m2 op, d, ar, br, ai, bi
    vwmul.vv v24, \ar, \br
    vwmul.vv v28, \ai, \bi
    vsetvli zero, zero, e32, m4, ta, ma
    \op\().vv v24, v24, v28
    vssra.vi v24, v24, 15
    vsetvli zero, zero, e16, m2, ta, ma
    vnsrl.wi \d, v24, 0
.endm

/* void spx_kf_rvv_bfly3_i16(kiss_fft_cpx *Fout, const kiss_fft_cpx *tw,
 *                           size_t fstride, int m, int inverse)
 * epi3.i = tw[fstride*m].i (a6); forward applies C_FIXDIV(.,3). */
    .balign 4
    .globl  spx_kf_rvv_bfly3_i16
    .hidden spx_kf_rvv_bfly3_i16
    .type   spx_kf_rvv_bfly3_i16, @function
spx_kf_rvv_bfly3_i16:
    csrwi   vxrm, 0
    mul     t2, a2, a3
    slli    t2, t2, 2
    add     t2, a1, t2
    lh      a6, 2(t2)               /* epi3.i */
    li      a5, 10922               /* SAMP_MAX/3 */
    slli    t0, a3, 2               /* m*4 */
    slli    t1, a2, 2               /* fstride*4 */
    mv      t2, a0                  /* Fout run */
    mv      t3, a1                  /* tw1 run */
    mv      t4, a3                  /* rem = m */
    vsetvli a7, zero, e16, m1, ta, ma  /* a7 = VLMAX at m1 */
    bgeu    a7, a3, .Lb3x_strip     /* m fits one m1 register */
    slli    t5, a7, 2
    bgeu    t5, a3, .Lb3x_m2        /* <= 2 m2 strips: LMUL=2 wins; beyond
                                       that the batched m1 pipeline is faster
                                       (measured on SpacemiT X60 at m=75) */
.Lb3x_strip:
    vsetvli a2, t4, e16, m1, ta, ma
    vlseg2e16.v  v0, (t2)           /* F0 */
    add     a7, t2, t0
    vlseg2e16.v  v2, (a7)           /* F1 */
    add     a7, a7, t0
    vlseg2e16.v  v4, (a7)           /* F2 */
    vlsseg2e16.v v6, (t3), t1       /* T1 */
    sub     t5, t3, a1
    add     t6, t3, t5              /* tw2 */
    slli    a7, t1, 1
    vlsseg2e16.v v8, (t6), a7       /* T2, stride 2*fstride */
    bnez    a4, 3f                  /* inverse: no C_FIXDIV */
    vwmul.vx v16, v0, a5
    vwmul.vx v18, v1, a5
    vwmul.vx v20, v2, a5
    vwmul.vx v22, v3, a5
    vwmul.vx v24, v4, a5
    vwmul.vx v26, v5, a5
    vsetvli zero, zero, e32, m2, ta, ma
    vssra.vi v16, v16, 15
    vssra.vi v18, v18, 15
    vssra.vi v20, v20, 15
    vssra.vi v22, v22, 15
    vssra.vi v24, v24, 15
    vssra.vi v26, v26, 15
    vsetvli zero, zero, e16, m1, ta, ma
    vnsrl.wi v0, v16, 0
    vnsrl.wi v1, v18, 0
    vnsrl.wi v2, v20, 0
    vnsrl.wi v3, v22, 0
    vnsrl.wi v4, v24, 0
    vnsrl.wi v5, v26, 0
3:
    vwmul.vv v16, v2, v6            /* s1 = C_MUL(F1, T1) */
    vwmul.vv v18, v3, v7
    vwmul.vv v20, v2, v7
    vwmul.vv v22, v3, v6
    vwmul.vv v24, v4, v8            /* s2 = C_MUL(F2, T2) */
    vwmul.vv v26, v5, v9
    vwmul.vv v28, v4, v9
    vwmul.vv v30, v5, v8
    vsetvli zero, zero, e32, m2, ta, ma
    vsub.vv  v16, v16, v18
    vadd.vv  v20, v20, v22
    vsub.vv  v24, v24, v26
    vadd.vv  v28, v28, v30
    vssra.vi v16, v16, 15
    vssra.vi v20, v20, 15
    vssra.vi v24, v24, 15
    vssra.vi v28, v28, 15
    vsetvli zero, zero, e16, m1, ta, ma
    vnsrl.wi v10, v16, 0            /* s1.r */
    vnsrl.wi v11, v20, 0            /* s1.i */
    vnsrl.wi v12, v24, 0            /* s2.r */
    vnsrl.wi v13, v28, 0            /* s2.i */
    vadd.vv  v14, v10, v12          /* s3 = s1 + s2 */
    vadd.vv  v15, v11, v13
    vsub.vv  v10, v10, v12          /* s0 = s1 - s2 */
    vsub.vv  v11, v11, v13
    vsra.vi  v2, v14, 1             /* Fm = F0 - HALF_OF(s3) */
    vsub.vv  v2, v0, v2
    vsra.vi  v3, v15, 1
    vsub.vv  v3, v1, v3
    vwmul.vx v16, v10, a6           /* s0 = C_MULBYSCALAR(s0, epi3.i) */
    vwmul.vx v18, v11, a6
    vsetvli zero, zero, e32, m2, ta, ma
    vssra.vi v16, v16, 15
    vssra.vi v18, v18, 15
    vsetvli zero, zero, e16, m1, ta, ma
    vnsrl.wi v10, v16, 0
    vnsrl.wi v11, v18, 0
    vadd.vv  v0, v0, v14            /* F0 += s3 */
    vadd.vv  v1, v1, v15
    vadd.vv  v4, v2, v11            /* Fout[2m] = (Fm.r + s0.i, Fm.i - s0.r) */
    vsub.vv  v5, v3, v10
    vsub.vv  v2, v2, v11            /* Fout[m]  = (Fm.r - s0.i, Fm.i + s0.r) */
    vadd.vv  v3, v3, v10
    vsseg2e16.v  v0, (t2)
    add     a7, t2, t0
    vsseg2e16.v  v2, (a7)
    add     a7, a7, t0
    vsseg2e16.v  v4, (a7)
    slli    a7, a2, 2
    add     t2, t2, a7
    mul     a7, a2, t1
    add     t3, t3, a7
    sub     t4, t4, a2
    bnez    t4, .Lb3x_strip
    ret
/* LMUL=2 body, taken when m*2 > vlenb.  Wide temporaries are e32,m4 here,
   so the widen/round/narrow batches run through two m4 groups (v24-27,
   v28-31; plus v20-23 while the s1 slots are still free) instead of the m1
   body's v16-v31 sweep.  Element arithmetic is identical (bit-exact). */
.Lb3x_m2:
    vsetvli a2, t4, e16, m2, ta, ma
    vlseg2e16.v  v0, (t2)           /* F0 */
    add     a7, t2, t0
    vlseg2e16.v  v4, (a7)           /* F1 */
    add     a7, a7, t0
    vlseg2e16.v  v8, (a7)           /* F2 */
    vlsseg2e16.v v12, (t3), t1      /* T1 */
    sub     t5, t3, a1
    add     t6, t3, t5              /* tw2 */
    slli    a7, t1, 1
    vlsseg2e16.v v16, (t6), a7      /* T2, stride 2*fstride */
    bnez    a4, 31f                 /* inverse: no C_FIXDIV */
    vwmul.vx v20, v0, a5
    vwmul.vx v24, v2, a5
    vwmul.vx v28, v4, a5
    vsetvli zero, zero, e32, m4, ta, ma
    vssra.vi v20, v20, 15
    vssra.vi v24, v24, 15
    vssra.vi v28, v28, 15
    vsetvli zero, zero, e16, m2, ta, ma
    vnsrl.wi v0, v20, 0
    vnsrl.wi v2, v24, 0
    vnsrl.wi v4, v28, 0
    vwmul.vx v20, v6, a5
    vwmul.vx v24, v8, a5
    vwmul.vx v28, v10, a5
    vsetvli zero, zero, e32, m4, ta, ma
    vssra.vi v20, v20, 15
    vssra.vi v24, v24, 15
    vssra.vi v28, v28, 15
    vsetvli zero, zero, e16, m2, ta, ma
    vnsrl.wi v6, v20, 0
    vnsrl.wi v8, v24, 0
    vnsrl.wi v10, v28, 0
31:
    kf_i16_cmulh_m2 vsub, v20, v4, v12, v6, v14   /* s1.r = C_MUL(F1,T1).r */
    kf_i16_cmulh_m2 vadd, v22, v4, v14, v6, v12   /* s1.i */
    kf_i16_cmulh_m2 vsub, v4, v8, v16, v10, v18   /* s2.r = C_MUL(F2,T2).r */
    kf_i16_cmulh_m2 vadd, v6, v8, v18, v10, v16   /* s2.i */
    vadd.vv  v8, v20, v4            /* s3 = s1 + s2 */
    vadd.vv  v10, v22, v6
    vsub.vv  v20, v20, v4           /* s0 = s1 - s2 */
    vsub.vv  v22, v22, v6
    vsra.vi  v4, v8, 1              /* Fm = F0 - HALF_OF(s3) */
    vsub.vv  v4, v0, v4
    vsra.vi  v6, v10, 1
    vsub.vv  v6, v2, v6
    vwmul.vx v24, v20, a6           /* s0 = C_MULBYSCALAR(s0, epi3.i) */
    vwmul.vx v28, v22, a6
    vsetvli zero, zero, e32, m4, ta, ma
    vssra.vi v24, v24, 15
    vssra.vi v28, v28, 15
    vsetvli zero, zero, e16, m2, ta, ma
    vnsrl.wi v20, v24, 0
    vnsrl.wi v22, v28, 0
    vadd.vv  v0, v0, v8             /* F0 += s3 */
    vadd.vv  v2, v2, v10
    vadd.vv  v8, v4, v22            /* Fout[2m] = (Fm.r + s0.i, Fm.i - s0.r) */
    vsub.vv  v10, v6, v20
    vsub.vv  v4, v4, v22            /* Fout[m]  = (Fm.r - s0.i, Fm.i + s0.r) */
    vadd.vv  v6, v6, v20
    vsseg2e16.v  v0, (t2)
    add     a7, t2, t0
    vsseg2e16.v  v4, (a7)
    add     a7, a7, t0
    vsseg2e16.v  v8, (a7)
    slli    a7, a2, 2
    add     t2, t2, a7
    mul     a7, a2, t1
    add     t3, t3, a7
    sub     t4, t4, a2
    bnez    t4, .Lb3x_m2
    ret
    .size spx_kf_rvv_bfly3_i16, . - spx_kf_rvv_bfly3_i16

/* ------------------------------------------------------------------ */
/* S_MUL pair: v14 = sround(\va * \xa), v15 = sround(\vb * \xb). */
.macro kf_b5x_smul2 va, xa, vb, xb
    vwmul.vx v24, \va, \xa
    vwmul.vx v26, \vb, \xb
    vsetvli zero, zero, e32, m2, ta, ma
    vssra.vi v24, v24, 15
    vssra.vi v26, v26, 15
    vsetvli zero, zero, e16, m1, ta, ma
    vnsrl.wi v14, v24, 0
    vnsrl.wi v15, v26, 0
.endm

/* LMUL=2 S_MUL pair: v24 = sround(\va * \xa), v28 = sround(\vb * \xb)
   (in-place narrow into the low half of each wide group). */
.macro kf_b5x_smul2_m2 va, xa, vb, xb
    vwmul.vx v24, \va, \xa
    vwmul.vx v28, \vb, \xb
    vsetvli zero, zero, e32, m4, ta, ma
    vssra.vi v24, v24, 15
    vssra.vi v28, v28, 15
    vsetvli zero, zero, e16, m2, ta, ma
    vnsrl.wi v24, v24, 0
    vnsrl.wi v28, v28, 0
.endm

/* LMUL=2 C_MUL: (\dr,\di) = sround((v16,v18) * (v20,v22)); clobbers v24-v31. */
.macro kf_b5x_cmul_m2 dr, di
    kf_i16_cmulh_m2 vsub, \dr, v16, v20, v18, v22
    kf_i16_cmulh_m2 vadd, \di, v16, v22, v18, v20
.endm

/* LMUL=2 C_FIXDIV(.,5) of the leg in (v16,v18); clobbers a7, v24-v31. */
.macro kf_b5x_fixdiv_m2
    li       a7, 6553               /* SAMP_MAX/5 */
    vwmul.vx v24, v16, a7
    vwmul.vx v28, v18, a7
    vsetvli zero, zero, e32, m4, ta, ma
    vssra.vi v24, v24, 15
    vssra.vi v28, v28, 15
    vsetvli zero, zero, e16, m2, ta, ma
    vnsrl.wi v16, v24, 0
    vnsrl.wi v18, v28, 0
.endm

/* void spx_kf_rvv_bfly5_i16(kiss_fft_cpx *Fout, const kiss_fft_cpx *tw,
 *                           size_t fstride, int m, int inverse)
 * ya = tw[fstride*m] (a5/a6), yb = tw[2*fstride*m] (a3/a0); forward
 * applies C_FIXDIV(.,5). */
    .balign 4
    .globl  spx_kf_rvv_bfly5_i16
    .hidden spx_kf_rvv_bfly5_i16
    .type   spx_kf_rvv_bfly5_i16, @function
spx_kf_rvv_bfly5_i16:
    csrwi   vxrm, 0
    mul     t2, a2, a3
    slli    t2, t2, 2
    add     t3, a1, t2
    lh      a5, 0(t3)               /* ya.r */
    lh      a6, 2(t3)               /* ya.i */
    add     t3, t3, t2
    lh      t5, 0(t3)
    lh      t6, 2(t3)
    slli    t0, a3, 2               /* m*4 */
    slli    t1, a2, 2               /* fstride*4 */
    mv      t2, a0                  /* Fout run */
    mv      t4, a3                  /* rem = m */
    mv      a3, t5                  /* yb.r */
    mv      a0, t6                  /* yb.i */
    mv      t3, a1                  /* tw1 run */
    srli    t5, t0, 2               /* m */
    vsetvli a7, zero, e16, m1, ta, ma  /* a7 = VLMAX at m1 */
    bgeu    a7, t5, .Lb5x_strip     /* m fits one m1 register */
    slli    a7, a7, 2
    bgeu    a7, t5, .Lb5x_m2        /* <= 2 m2 strips: LMUL=2 wins (see
                                       kf_bfly3_i16 note) */
.Lb5x_strip:
    vsetvli a2, t4, e16, m1, ta, ma
    vlseg2e16.v  v0, (t2)           /* F0 (scratch[0] stays in v0v1) */
    add     a7, t2, t0
    vlseg2e16.v  v2, (a7)
    add     a7, a7, t0
    vlseg2e16.v  v4, (a7)
    add     a7, a7, t0
    vlseg2e16.v  v6, (a7)
    add     a7, a7, t0
    vlseg2e16.v  v8, (a7)
    vlsseg2e16.v v10, (t3), t1      /* T1 */
    sub     t5, t3, a1              /* d */
    add     t6, t3, t5              /* tw2 */
    slli    a7, t1, 1
    vlsseg2e16.v v12, (t6), a7      /* T2, stride 2*fstride */
    add     t6, t6, t5              /* tw3 */
    add     a7, a7, t1
    vlsseg2e16.v v14, (t6), a7      /* T3, stride 3*fstride */
    add     t6, t6, t5              /* tw4 */
    add     a7, a7, t1
    vlsseg2e16.v v16, (t6), a7      /* T4, stride 4*fstride */
    bnez    a4, 3f                  /* inverse: no C_FIXDIV */
    li      a7, 6553                /* SAMP_MAX/5 */
    vwmul.vx v20, v0, a7
    vwmul.vx v22, v1, a7
    vwmul.vx v24, v2, a7
    vwmul.vx v26, v3, a7
    vwmul.vx v28, v4, a7
    vwmul.vx v30, v5, a7
    vsetvli zero, zero, e32, m2, ta, ma
    vssra.vi v20, v20, 15
    vssra.vi v22, v22, 15
    vssra.vi v24, v24, 15
    vssra.vi v26, v26, 15
    vssra.vi v28, v28, 15
    vssra.vi v30, v30, 15
    vsetvli zero, zero, e16, m1, ta, ma
    vnsrl.wi v0, v20, 0
    vnsrl.wi v1, v22, 0
    vnsrl.wi v2, v24, 0
    vnsrl.wi v3, v26, 0
    vnsrl.wi v4, v28, 0
    vnsrl.wi v5, v30, 0
    vwmul.vx v20, v6, a7
    vwmul.vx v22, v7, a7
    vwmul.vx v24, v8, a7
    vwmul.vx v26, v9, a7
    vsetvli zero, zero, e32, m2, ta, ma
    vssra.vi v20, v20, 15
    vssra.vi v22, v22, 15
    vssra.vi v24, v24, 15
    vssra.vi v26, v26, 15
    vsetvli zero, zero, e16, m1, ta, ma
    vnsrl.wi v6, v20, 0
    vnsrl.wi v7, v22, 0
    vnsrl.wi v8, v24, 0
    vnsrl.wi v9, v26, 0
3:
    /* s_k = C_MUL(F_k, T_k) in place, k = 1..4 */
    vwmul.vv v20, v2, v10
    vwmul.vv v22, v3, v11
    vwmul.vv v24, v2, v11
    vwmul.vv v26, v3, v10
    vwmul.vv v28, v4, v12
    vwmul.vv v30, v5, v13
    vsetvli zero, zero, e32, m2, ta, ma
    vsub.vv  v20, v20, v22
    vadd.vv  v24, v24, v26
    vssra.vi v20, v20, 15
    vssra.vi v24, v24, 15
    vsetvli zero, zero, e16, m1, ta, ma
    vnsrl.wi v2, v20, 0             /* s1.r */
    vnsrl.wi v3, v24, 0             /* s1.i */
    vwmul.vv v20, v4, v13
    vwmul.vv v22, v5, v12
    vsetvli zero, zero, e32, m2, ta, ma
    vsub.vv  v28, v28, v30
    vadd.vv  v20, v20, v22
    vssra.vi v28, v28, 15
    vssra.vi v20, v20, 15
    vsetvli zero, zero, e16, m1, ta, ma
    vnsrl.wi v4, v28, 0             /* s2.r */
    vnsrl.wi v5, v20, 0             /* s2.i */
    vwmul.vv v20, v6, v14
    vwmul.vv v22, v7, v15
    vwmul.vv v24, v6, v15
    vwmul.vv v26, v7, v14
    vwmul.vv v28, v8, v16
    vwmul.vv v30, v9, v17
    vsetvli zero, zero, e32, m2, ta, ma
    vsub.vv  v20, v20, v22
    vadd.vv  v24, v24, v26
    vssra.vi v20, v20, 15
    vssra.vi v24, v24, 15
    vsetvli zero, zero, e16, m1, ta, ma
    vnsrl.wi v6, v20, 0             /* s3.r */
    vnsrl.wi v7, v24, 0             /* s3.i */
    vwmul.vv v20, v8, v17
    vwmul.vv v22, v9, v16
    vsetvli zero, zero, e32, m2, ta, ma
    vsub.vv  v28, v28, v30
    vadd.vv  v20, v20, v22
    vssra.vi v28, v28, 15
    vssra.vi v20, v20, 15
    vsetvli zero, zero, e16, m1, ta, ma
    vnsrl.wi v8, v28, 0             /* s4.r */
    vnsrl.wi v9, v20, 0             /* s4.i */
    /* combines (all wrapping e16, matching C's int arith + int16 store) */
    vadd.vv  v10, v2, v8            /* s7 = s1 + s4 */
    vadd.vv  v11, v3, v9
    vsub.vv  v2, v2, v8             /* s10 = s1 - s4 */
    vsub.vv  v3, v3, v9
    vadd.vv  v12, v4, v6            /* s8 = s2 + s3 */
    vadd.vv  v13, v5, v7
    vsub.vv  v4, v4, v6             /* s9 = s2 - s3 */
    vsub.vv  v5, v5, v7
    vadd.vv  v18, v0, v10           /* Fout += s7 + s8 */
    vadd.vv  v18, v18, v12
    vadd.vv  v19, v1, v11
    vadd.vv  v19, v19, v13
    vsseg2e16.v  v18, (t2)
    kf_b5x_smul2 v10, a5, v12, a3   /* ya.r*s7.r, yb.r*s8.r */
    vadd.vv  v20, v0, v14           /* s5.r = s0.r + . + . */
    vadd.vv  v20, v20, v15
    kf_b5x_smul2 v11, a5, v13, a3
    vadd.vv  v21, v1, v14           /* s5.i */
    vadd.vv  v21, v21, v15
    kf_b5x_smul2 v3, a6, v5, a0     /* ya.i*s10.i, yb.i*s9.i */
    vadd.vv  v22, v14, v15          /* s6.r */
    kf_b5x_smul2 v2, a6, v4, a0     /* ya.i*s10.r, yb.i*s9.r */
    vadd.vv  v23, v14, v15          /* s6.i = -(. + .) */
    vrsub.vi v23, v23, 0
    vsub.vv  v14, v20, v22          /* Fout[m] = s5 - s6 */
    vsub.vv  v15, v21, v23
    add     a7, t2, t0
    vsseg2e16.v  v14, (a7)
    vadd.vv  v14, v20, v22          /* Fout[4m] = s5 + s6 */
    vadd.vv  v15, v21, v23
    slli    a7, t0, 2
    add     a7, t2, a7
    vsseg2e16.v  v14, (a7)
    kf_b5x_smul2 v10, a3, v12, a5   /* yb.r*s7.r, ya.r*s8.r */
    vadd.vv  v20, v0, v14           /* s11.r */
    vadd.vv  v20, v20, v15
    kf_b5x_smul2 v11, a3, v13, a5
    vadd.vv  v21, v1, v14           /* s11.i */
    vadd.vv  v21, v21, v15
    kf_b5x_smul2 v3, a0, v5, a6     /* yb.i*s10.i, ya.i*s9.i */
    vsub.vv  v22, v15, v14          /* s12.r = -n1 + n2 */
    kf_b5x_smul2 v2, a0, v4, a6     /* yb.i*s10.r, ya.i*s9.r */
    vsub.vv  v23, v14, v15          /* s12.i = n1 - n2 */
    vadd.vv  v14, v20, v22          /* Fout[2m] = s11 + s12 */
    vadd.vv  v15, v21, v23
    slli    a7, t0, 1
    add     a7, t2, a7
    vsseg2e16.v  v14, (a7)
    vsub.vv  v14, v20, v22          /* Fout[3m] = s11 - s12 */
    vsub.vv  v15, v21, v23
    add     a7, a7, t0
    vsseg2e16.v  v14, (a7)
    slli    a7, a2, 2
    add     t2, t2, a7
    mul     a7, a2, t1
    add     t3, t3, a7
    sub     t4, t4, a2
    bnez    t4, .Lb5x_strip
    ret
/* LMUL=2 body, taken when m*2 > vlenb.  Wide temporaries are e32,m4, so
   only two wide groups (v24-27, v28-31) fit; each leg is loaded into
   (v16,v18)/(v20,v22), C_FIXDIVed (forward) and twiddle-multiplied into
   its s-slot before the next leg loads, and F0 is reloaded afterwards
   instead of staying live.  Element arithmetic is identical (bit-exact).
   s-slots: s1=(v0,v2) s2=(v4,v6) s3=(v8,v10) s4=(v12,v14). */
.Lb5x_m2:
    vsetvli a2, t4, e16, m2, ta, ma
    add     a7, t2, t0              /* leg 1 */
    vlseg2e16.v  v16, (a7)
    mv      t6, t3
    vlsseg2e16.v v20, (t6), t1
    bnez    a4, 11f
    kf_b5x_fixdiv_m2
11: kf_b5x_cmul_m2 v0, v2
    slli    a7, t0, 1               /* leg 2 */
    add     a7, t2, a7
    vlseg2e16.v  v16, (a7)
    sub     t5, t3, a1              /* d */
    add     t6, t6, t5
    slli    a7, t1, 1
    vlsseg2e16.v v20, (t6), a7      /* stride 2*fstride */
    bnez    a4, 12f
    kf_b5x_fixdiv_m2
12: kf_b5x_cmul_m2 v4, v6
    slli    a7, t0, 1               /* leg 3 */
    add     a7, a7, t0
    add     a7, t2, a7
    vlseg2e16.v  v16, (a7)
    add     t6, t6, t5
    slli    a7, t1, 1
    add     a7, a7, t1
    vlsseg2e16.v v20, (t6), a7      /* stride 3*fstride */
    bnez    a4, 13f
    kf_b5x_fixdiv_m2
13: kf_b5x_cmul_m2 v8, v10
    slli    a7, t0, 2               /* leg 4 */
    add     a7, t2, a7
    vlseg2e16.v  v16, (a7)
    add     t6, t6, t5
    slli    a7, t1, 2
    vlsseg2e16.v v20, (t6), a7      /* stride 4*fstride */
    bnez    a4, 14f
    kf_b5x_fixdiv_m2
14: kf_b5x_cmul_m2 v12, v14
    vlseg2e16.v  v16, (t2)          /* F0 (reload) */
    bnez    a4, 15f
    kf_b5x_fixdiv_m2
15:
    vadd.vv  v20, v0, v12           /* s7 = s1 + s4 */
    vadd.vv  v22, v2, v14
    vsub.vv  v0, v0, v12            /* s10 = s1 - s4 */
    vsub.vv  v2, v2, v14
    vadd.vv  v12, v4, v8            /* s8 = s2 + s3 */
    vadd.vv  v14, v6, v10
    vsub.vv  v4, v4, v8             /* s9 = s2 - s3 */
    vsub.vv  v6, v6, v10
    vadd.vv  v8, v16, v20           /* Fout += s7 + s8 */
    vadd.vv  v8, v8, v12
    vadd.vv  v10, v18, v22
    vadd.vv  v10, v10, v14
    vsseg2e16.v  v8, (t2)
    kf_b5x_smul2_m2 v20, a5, v12, a3   /* ya.r*s7.r, yb.r*s8.r */
    vadd.vv  v8, v24, v28
    vadd.vv  v8, v16, v8            /* s5.r */
    kf_b5x_smul2_m2 v22, a5, v14, a3
    vadd.vv  v10, v24, v28
    vadd.vv  v10, v18, v10          /* s5.i */
    kf_b5x_smul2_m2 v20, a3, v12, a5   /* yb.r*s7.r, ya.r*s8.r */
    vadd.vv  v20, v24, v28
    vadd.vv  v20, v16, v20          /* s11.r */
    kf_b5x_smul2_m2 v22, a3, v14, a5
    vadd.vv  v22, v24, v28
    vadd.vv  v22, v18, v22          /* s11.i */
    kf_b5x_smul2_m2 v2, a6, v6, a0     /* ya.i*s10.i, yb.i*s9.i */
    vadd.vv  v12, v24, v28          /* s6.r */
    kf_b5x_smul2_m2 v0, a6, v4, a0     /* ya.i*s10.r, yb.i*s9.r */
    vadd.vv  v14, v24, v28
    vrsub.vi v14, v14, 0            /* s6.i = -(. + .) */
    vsub.vv  v16, v8, v12           /* Fout[m] = s5 - s6 */
    vsub.vv  v18, v10, v14
    add     a7, t2, t0
    vsseg2e16.v  v16, (a7)
    vadd.vv  v16, v8, v12           /* Fout[4m] = s5 + s6 */
    vadd.vv  v18, v10, v14
    slli    a7, t0, 2
    add     a7, t2, a7
    vsseg2e16.v  v16, (a7)
    kf_b5x_smul2_m2 v2, a0, v6, a6     /* yb.i*s10.i, ya.i*s9.i */
    vsub.vv  v8, v28, v24           /* s12.r = -n1 + n2 */
    kf_b5x_smul2_m2 v0, a0, v4, a6     /* yb.i*s10.r, ya.i*s9.r */
    vsub.vv  v10, v24, v28          /* s12.i = n1 - n2 */
    vadd.vv  v16, v20, v8           /* Fout[2m] = s11 + s12 */
    vadd.vv  v18, v22, v10
    slli    a7, t0, 1
    add     a7, t2, a7
    vsseg2e16.v  v16, (a7)
    vsub.vv  v16, v20, v8           /* Fout[3m] = s11 - s12 */
    vsub.vv  v18, v22, v10
    add     a7, a7, t0
    vsseg2e16.v  v16, (a7)
    slli    a7, a2, 2
    add     t2, t2, a7
    mul     a7, a2, t1
    add     t3, t3, a7
    sub     t4, t4, a2
    bnez    t4, .Lb5x_m2
    ret
    .size spx_kf_rvv_bfly5_i16, . - spx_kf_rvv_bfly5_i16

#if defined(__ELF__)
/* Non-executable stack note */
    .section .note.GNU-stack, "", @progbits
#endif
