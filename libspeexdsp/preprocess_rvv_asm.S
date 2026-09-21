/* Copyright (C) 2003 Epic Games (written by Jean-Marc Valin)
 * Copyright (C) 2004-2006 Epic Games
 * Copyright (C) 2026 Tristan Matthews
 */
/*
   preprocess_rvv_asm.S -- hand-written RVV 1.0 kernels for the
   preprocessor (denoiser).

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

/* RVV 1.0 kernels (VLEN >= 128) for preprocess.c's per-bin loops, float
 * builds only (the fixed-point twins are branchy range-reduced divides
 * and polynomial approximations with no clean vector mapping). Most
 * kernels are element-wise over the N (or N+M) spectral bins and
 * strip-mine with LMUL=4/8; the two packed-spectrum kernels
 * (power_spectrum, apply_gain) walk the interior complex bins with
 * vlseg2/vsseg2 and leave the real-only edge bins to the C wrapper in
 * preprocess_rvv.h, and the two 3-tap smoothers (smooth_spectrum, zeta)
 * update the interior bins [1, N-2] with three offset streams.
 *
 * The Ephraim-Malah gain kernel vectorizes hypergeom_gain: trunc(2x)
 * (== floor, x >= 0) via vfcvt.rtz, the 21-entry table pair via two
 * vluxei32 gathers on a clamped index, and the ind < 0 / ind > 19 tails
 * via mask merges; sqrt and the divides map to vfsqrt/vfdiv. Float
 * results pair with a small checkasm tolerance (FMAs, 1/x refactoring,
 * reordered sums).
 *
 * `.option arch, +v` keeps the file self-contained (no -march needed).
 * Leaf functions, caller-saved registers only. The wrappers guarantee
 * N >= 4 and in-bounds loads. */

    .option arch, +v
    .text

/* ------------------------------------------------------------------ */
/* void spx_preproc_rvv_window_f32(float *frame, const float *window,
 *                                 int len)
 * frame[i] *= window[i] -- analysis/synthesis windowing. */
    .balign 4
    .globl  spx_preproc_rvv_window_f32
    .hidden spx_preproc_rvv_window_f32
    .type   spx_preproc_rvv_window_f32, @function
spx_preproc_rvv_window_f32:
1:  vsetvli t0, a2, e32, m8, ta, ma
    vle32.v v8, (a0)
    vle32.v v16, (a1)
    vfmul.vv v8, v8, v16
    vse32.v v8, (a0)
    sub     a2, a2, t0
    slli    t1, t0, 2
    add     a0, a0, t1
    add     a1, a1, t1
    bnez    a2, 1b
    ret
    .size spx_preproc_rvv_window_f32, . - spx_preproc_rvv_window_f32

/* ------------------------------------------------------------------ */
/* void spx_preproc_rvv_power_spectrum_f32(const float *ft, float *ps,
 *                                         int N)
 * ps[1..N-1] = re^2 + im^2 of the N-1 interior complex bins of the
 * packed 2N-point spectrum ft; the DC bin ps[0] stays in C. */
    .balign 4
    .globl  spx_preproc_rvv_power_spectrum_f32
    .hidden spx_preproc_rvv_power_spectrum_f32
    .type   spx_preproc_rvv_power_spectrum_f32, @function
spx_preproc_rvv_power_spectrum_f32:
    addi    a0, a0, 4              /* skip the DC bin */
    addi    a1, a1, 4
    addi    a2, a2, -1             /* N-1 pairs */
1:  vsetvli t0, a2, e32, m4, ta, ma
    vlseg2e32.v v8, (a0)           /* re v8, im v12 */
    vfmul.vv v16, v8, v8
    vfmacc.vv v16, v12, v12
    vse32.v v16, (a1)
    sub     a2, a2, t0
    slli    t1, t0, 3
    add     a0, a0, t1
    slli    t1, t0, 2
    add     a1, a1, t1
    bnez    a2, 1b
    ret
    .size spx_preproc_rvv_power_spectrum_f32, . - spx_preproc_rvv_power_spectrum_f32

/* ------------------------------------------------------------------ */
/* void spx_preproc_rvv_smooth_spectrum_f32(float *S, const float *ps,
 *                                          int N)
 * S[i] = .8*S[i] + .05*ps[i-1] + .1*ps[i] + .05*ps[i+1] over the
 * interior bins [1, N-2]; the edge bins stay in C. */
    .balign 4
    .globl  spx_preproc_rvv_smooth_spectrum_f32
    .hidden spx_preproc_rvv_smooth_spectrum_f32
    .type   spx_preproc_rvv_smooth_spectrum_f32, @function
spx_preproc_rvv_smooth_spectrum_f32:
    li      t1, 0x3f4ccccd
    fmv.w.x fa0, t1                /* .8f */
    li      t1, 0x3dcccccd
    fmv.w.x fa1, t1                /* .1f */
    li      t1, 0x3d4ccccd
    fmv.w.x fa2, t1                /* .05f */
    addi    a2, a2, -2             /* interior bin count */
    mv      a3, a1                 /* ps[i-1] */
    addi    a4, a1, 4              /* ps[i]   */
    addi    a5, a1, 8              /* ps[i+1] */
    addi    a0, a0, 4
1:  vsetvli t0, a2, e32, m4, ta, ma
    vle32.v v8, (a0)               /* S */
    vle32.v v12, (a3)
    vle32.v v16, (a4)
    vle32.v v20, (a5)
    vfmul.vf v8, v8, fa0
    vfmacc.vf v8, fa1, v16
    vfadd.vv v12, v12, v20
    vfmacc.vf v8, fa2, v12
    vse32.v v8, (a0)
    sub     a2, a2, t0
    slli    t1, t0, 2
    add     a0, a0, t1
    add     a3, a3, t1
    add     a4, a4, t1
    add     a5, a5, t1
    bnez    a2, 1b
    ret
    .size spx_preproc_rvv_smooth_spectrum_f32, . - spx_preproc_rvv_smooth_spectrum_f32

/* ------------------------------------------------------------------ */
/* void spx_preproc_rvv_min_track_swap_f32(float *Smin, float *Stmp,
 *                                         const float *S, int N)
 * Smin = min(Stmp, S); Stmp = S -- minimum statistics window restart. */
    .balign 4
    .globl  spx_preproc_rvv_min_track_swap_f32
    .hidden spx_preproc_rvv_min_track_swap_f32
    .type   spx_preproc_rvv_min_track_swap_f32, @function
spx_preproc_rvv_min_track_swap_f32:
1:  vsetvli t0, a3, e32, m8, ta, ma
    vle32.v v8, (a1)               /* Stmp */
    vle32.v v16, (a2)              /* S */
    vfmin.vv v24, v8, v16
    vse32.v v24, (a0)
    vse32.v v16, (a1)
    sub     a3, a3, t0
    slli    t1, t0, 2
    add     a0, a0, t1
    add     a1, a1, t1
    add     a2, a2, t1
    bnez    a3, 1b
    ret
    .size spx_preproc_rvv_min_track_swap_f32, . - spx_preproc_rvv_min_track_swap_f32

/* ------------------------------------------------------------------ */
/* void spx_preproc_rvv_min_track_f32(float *Smin, float *Stmp,
 *                                    const float *S, int N)
 * Smin = min(Smin, S); Stmp = min(Stmp, S) -- running minimum. */
    .balign 4
    .globl  spx_preproc_rvv_min_track_f32
    .hidden spx_preproc_rvv_min_track_f32
    .type   spx_preproc_rvv_min_track_f32, @function
spx_preproc_rvv_min_track_f32:
1:  vsetvli t0, a3, e32, m8, ta, ma
    vle32.v v16, (a2)              /* S */
    vle32.v v8, (a0)
    vfmin.vv v8, v8, v16
    vse32.v v8, (a0)
    vle32.v v24, (a1)
    vfmin.vv v24, v24, v16
    vse32.v v24, (a1)
    sub     a3, a3, t0
    slli    t1, t0, 2
    add     a0, a0, t1
    add     a1, a1, t1
    add     a2, a2, t1
    bnez    a3, 1b
    ret
    .size spx_preproc_rvv_min_track_f32, . - spx_preproc_rvv_min_track_f32

/* ------------------------------------------------------------------ */
/* void spx_preproc_rvv_update_prob_f32(const float *S, const float *Smin,
 *                                      int *update_prob, int N)
 * update_prob[i] = (.4*S[i] > Smin[i]) ? 1 : 0 */
    .balign 4
    .globl  spx_preproc_rvv_update_prob_f32
    .hidden spx_preproc_rvv_update_prob_f32
    .type   spx_preproc_rvv_update_prob_f32, @function
spx_preproc_rvv_update_prob_f32:
    li      t1, 0x3ecccccd
    fmv.w.x fa0, t1                /* .4f */
1:  vsetvli t0, a3, e32, m4, ta, ma
    vle32.v v8, (a0)
    vle32.v v12, (a1)
    vfmul.vf v8, v8, fa0
    vmflt.vv v0, v12, v8           /* Smin < .4*S */
    vmv.v.i v16, 0
    vmerge.vim v16, v16, 1, v0
    vse32.v v16, (a2)
    sub     a3, a3, t0
    slli    t1, t0, 2
    add     a0, a0, t1
    add     a1, a1, t1
    add     a2, a2, t1
    bnez    a3, 1b
    ret
    .size spx_preproc_rvv_update_prob_f32, . - spx_preproc_rvv_update_prob_f32

/* ------------------------------------------------------------------ */
/* void spx_preproc_rvv_noise_update_f32(const int *update_prob,
 *                                       const float *ps, float *noise,
 *                                       float beta, float beta_1, int N)
 * Where !update_prob[i] || ps[i] < noise[i]:
 *    noise[i] = max(0, beta_1*noise[i] + beta*ps[i]) */
    .balign 4
    .globl  spx_preproc_rvv_noise_update_f32
    .hidden spx_preproc_rvv_noise_update_f32
    .type   spx_preproc_rvv_noise_update_f32, @function
spx_preproc_rvv_noise_update_f32:
    fmv.w.x fa2, zero              /* 0.0f */
1:  vsetvli t0, a3, e32, m4, ta, ma
    vle32.v v8, (a0)               /* update_prob (int32) */
    vle32.v v12, (a1)              /* ps */
    vle32.v v16, (a2)              /* noise */
    vmseq.vx v0, v8, zero          /* !update_prob */
    vmflt.vv v1, v12, v16          /* ps < noise */
    vmor.mm v0, v0, v1
    vfmul.vf v20, v16, fa1
    vfmacc.vf v20, fa0, v12        /* beta_1*noise + beta*ps */
    vfmax.vf v20, v20, fa2
    vse32.v v20, (a2), v0.t
    sub     a3, a3, t0
    slli    t1, t0, 2
    add     a0, a0, t1
    add     a1, a1, t1
    add     a2, a2, t1
    bnez    a3, 1b
    ret
    .size spx_preproc_rvv_noise_update_f32, . - spx_preproc_rvv_noise_update_f32

/* ------------------------------------------------------------------ */
/* void spx_preproc_rvv_snr_f32(const float *ps, const float *noise,
 *                              const float *echo_noise,
 *                              const float *reverb_estimate,
 *                              const float *old_ps,
 *                              float *post, float *prior, int len)
 * tot   = 1 + noise + echo + reverb
 * post  = min(ps/tot - 1, 100)
 * gamma = .1 + .89*(old/(old+tot))^2
 * prior = min(gamma*max(0,post) + (1-gamma)*old/tot, 100)
 * ps/tot and old/tot go through r = 1/tot (one divide saved). */
    .balign 4
    .globl  spx_preproc_rvv_snr_f32
    .hidden spx_preproc_rvv_snr_f32
    .type   spx_preproc_rvv_snr_f32, @function
spx_preproc_rvv_snr_f32:
    li      t1, 0x3f800000
    fmv.w.x fa0, t1                /* 1.0f */
    li      t1, 0x42c80000
    fmv.w.x fa1, t1                /* 100.0f */
    li      t1, 0x3dcccccd
    fmv.w.x fa2, t1                /* .1f */
    li      t1, 0x3f63d70a
    fmv.w.x fa3, t1                /* .89f */
    fmv.w.x fa4, zero              /* 0.0f */
1:  vsetvli t0, a7, e32, m4, ta, ma
    vle32.v v4, (a0)               /* ps */
    vle32.v v8, (a1)               /* noise */
    vle32.v v12, (a2)              /* echo_noise */
    vle32.v v16, (a3)              /* reverb_estimate */
    vfadd.vv v8, v8, v12
    vfadd.vv v8, v8, v16
    vfadd.vf v8, v8, fa0           /* tot */
    vle32.v v12, (a4)              /* old_ps */
    vfrdiv.vf v16, v8, fa0         /* r = 1/tot */
    vfmul.vv v20, v4, v16
    vfsub.vf v20, v20, fa0         /* post = ps*r - 1 */
    vfmin.vf v20, v20, fa1
    vse32.v v20, (a5)
    vfadd.vv v24, v12, v8
    vfdiv.vv v24, v12, v24         /* old/(old+tot) */
    vfmul.vv v24, v24, v24
    vfmul.vf v24, v24, fa3
    vfadd.vf v24, v24, fa2         /* gamma */
    vfmax.vf v28, v20, fa4         /* max(0, post) */
    vfmul.vv v28, v24, v28
    vfrsub.vf v24, v24, fa0        /* 1-gamma */
    vfmul.vv v12, v12, v16         /* old*r */
    vfmacc.vv v28, v24, v12
    vfmin.vf v28, v28, fa1
    vse32.v v28, (a6)
    sub     a7, a7, t0
    slli    t1, t0, 2
    add     a0, a0, t1
    add     a1, a1, t1
    add     a2, a2, t1
    add     a3, a3, t1
    add     a4, a4, t1
    add     a5, a5, t1
    add     a6, a6, t1
    bnez    a7, 1b
    ret
    .size spx_preproc_rvv_snr_f32, . - spx_preproc_rvv_snr_f32

/* ------------------------------------------------------------------ */
/* void spx_preproc_rvv_zeta_f32(float *zeta, const float *prior, int N)
 * zeta[i] = .7*zeta[i] + .15*prior[i] + .075*prior[i-1] + .075*prior[i+1]
 * over the interior bins [1, N-2]; the edges stay in C. */
    .balign 4
    .globl  spx_preproc_rvv_zeta_f32
    .hidden spx_preproc_rvv_zeta_f32
    .type   spx_preproc_rvv_zeta_f32, @function
spx_preproc_rvv_zeta_f32:
    li      t1, 0x3f333333
    fmv.w.x fa0, t1                /* .7f */
    li      t1, 0x3e19999a
    fmv.w.x fa1, t1                /* .15f */
    li      t1, 0x3d99999a
    fmv.w.x fa2, t1                /* .075f */
    addi    a2, a2, -2             /* interior bin count */
    mv      a3, a1                 /* prior[i-1] */
    addi    a4, a1, 4              /* prior[i]   */
    addi    a5, a1, 8              /* prior[i+1] */
    addi    a0, a0, 4
1:  vsetvli t0, a2, e32, m4, ta, ma
    vle32.v v8, (a0)               /* zeta */
    vle32.v v12, (a3)
    vle32.v v16, (a4)
    vle32.v v20, (a5)
    vfmul.vf v8, v8, fa0
    vfmacc.vf v8, fa1, v16
    vfadd.vv v12, v12, v20
    vfmacc.vf v8, fa2, v12
    vse32.v v8, (a0)
    sub     a2, a2, t0
    slli    t1, t0, 2
    add     a0, a0, t1
    add     a3, a3, t1
    add     a4, a4, t1
    add     a5, a5, t1
    bnez    a2, 1b
    ret
    .size spx_preproc_rvv_zeta_f32, . - spx_preproc_rvv_zeta_f32

/* ------------------------------------------------------------------ */
/* void spx_preproc_rvv_em_gain_f32(const float *prior, const float *post,
 *                                  const float *ps, const float *gain_floor,
 *                                  float *gain, float *gain2,
 *                                  float *old_ps, int N)
 * Linear-frequency Ephraim-Malah gain loop:
 *   ratio = prior/(1+prior);  x = ratio*(1+post)
 *   MM    = hypergeom_gain(x)   (table gather + interpolation, with the
 *                                ind < 0 -> 1 and ind > 19 -> 1+.1296/x
 *                                tails mask-merged in)
 *   g     = min(1, ratio*MM);   g = 3*gain if .333*g > gain
 *   old_ps = .2*old_ps + .8*g^2*ps         (pre-floor g, as in C)
 *   gain  = max(g, gain_floor)
 *   gain2 = (p*sqrt(gain) + (1-p)*sqrt(gain_floor))^2,  p = gain2 in
 * LMUL=2: the strip needs 7 live register groups plus two gather index
 * temps around the table lookup. */
    .balign 4
    .globl  spx_preproc_rvv_em_gain_f32
    .hidden spx_preproc_rvv_em_gain_f32
    .type   spx_preproc_rvv_em_gain_f32, @function
spx_preproc_rvv_em_gain_f32:
    li      t1, 0x3f800000
    fmv.w.x fa0, t1                /* 1.0f */
    li      t1, 0x38d1b717
    fmv.w.x fa1, t1                /* 1e-4f */
    li      t1, 0x3e04b5dd
    fmv.w.x fa2, t1                /* .1296f */
    li      t1, 0x3eaa7efa
    fmv.w.x fa3, t1                /* .333f */
    li      t1, 0x40400000
    fmv.w.x fa4, t1                /* 3.0f */
    li      t1, 0x3e4ccccd
    fmv.w.x fa5, t1                /* .2f */
    li      t1, 0x3f4ccccd
    fmv.w.x fa6, t1                /* .8f */
    fmv.w.x fa7, zero              /* 0.0f */
    lla     t3, .Lhypergeom_table
    li      t4, 19
1:  vsetvli t0, a7, e32, m2, ta, ma
    vle32.v v4, (a0)               /* prior */
    vfadd.vf v8, v4, fa0
    vfdiv.vv v8, v4, v8            /* ratio = prior/(1+prior) */
    vle32.v v12, (a1)              /* post */
    vfadd.vf v12, v12, fa0
    vfmul.vv v12, v8, v12          /* x = ratio*(1+post) */
    vfadd.vv v16, v12, v12         /* 2x */
    vfcvt.rtz.x.f.v v20, v16       /* ind = trunc(2x) (== floor, x >= 0) */
    vfcvt.f.x.v v24, v20
    vfsub.vv v16, v16, v24         /* frac = 2x - ind */
    vmflt.vf v2, v12, fa7          /* tail: x < 0 -> MM = 1 */
    vmsgt.vx v3, v20, t4           /* tail: ind > 19 -> MM = 1+.1296/x */
    vmax.vx v20, v20, zero         /* clamp the gather index */
    vmin.vx v20, v20, t4
    vsll.vi v20, v20, 2
    vluxei32.v v24, (t3), v20      /* table[ind] */
    vadd.vi v20, v20, 4
    vluxei32.v v28, (t3), v20      /* table[ind+1] */
    vfsub.vv v28, v28, v24
    vfmacc.vv v24, v16, v28        /* lerp(table[ind], table[ind+1], frac) */
    vfadd.vf v16, v12, fa1
    vfsqrt.v v16, v16              /* sqrt(x + 1e-4) */
    vfdiv.vv v24, v24, v16
    vfrdiv.vf v16, v12, fa2
    vfadd.vf v16, v16, fa0         /* 1 + .1296/x */
    vmand.mm v0, v3, v3
    vmerge.vvm v24, v24, v16, v0
    vmand.mm v0, v2, v2
    vfmerge.vfm v24, v24, fa0, v0  /* MM */
    vfmul.vv v24, v8, v24
    vfmin.vf v24, v24, fa0         /* g = min(1, ratio*MM) */
    vle32.v v28, (a4)              /* gain (Bark-interpolated) */
    vfmul.vf v16, v24, fa3
    vmflt.vv v0, v28, v16          /* gain < .333*g */
    vfmul.vf v16, v28, fa4
    vmerge.vvm v24, v24, v16, v0   /* g = 3*gain on those lanes */
    vle32.v v8, (a6)               /* old_ps */
    vle32.v v12, (a2)              /* ps */
    vfmul.vv v16, v24, v24
    vfmul.vv v16, v16, v12         /* g^2*ps */
    vfmul.vf v8, v8, fa5
    vfmacc.vf v8, fa6, v16         /* .2*old_ps + .8*g^2*ps */
    vse32.v v8, (a6)
    vle32.v v20, (a3)              /* gain_floor */
    vfmax.vv v24, v24, v20         /* apply the gain floor */
    vse32.v v24, (a4)
    vfsqrt.v v24, v24
    vfsqrt.v v20, v20
    vle32.v v28, (a5)              /* p = gain2 in */
    vfmul.vv v24, v24, v28         /* p*sqrt(gain) */
    vfrsub.vf v28, v28, fa0
    vfmacc.vv v24, v28, v20        /* + (1-p)*sqrt(gain_floor) */
    vfmul.vv v24, v24, v24
    vse32.v v24, (a5)
    sub     a7, a7, t0
    slli    t1, t0, 2
    add     a0, a0, t1
    add     a1, a1, t1
    add     a2, a2, t1
    add     a3, a3, t1
    add     a4, a4, t1
    add     a5, a5, t1
    add     a6, a6, t1
    bnez    a7, 1b
    ret
    .size spx_preproc_rvv_em_gain_f32, . - spx_preproc_rvv_em_gain_f32

/* ------------------------------------------------------------------ */
/* void spx_preproc_rvv_apply_gain_f32(const float *gain2, float *ft,
 *                                     int N)
 * ft[2i-1] *= gain2[i]; ft[2i] *= gain2[i] for the N-1 interior complex
 * bins; the two real-only edge bins stay in C. */
    .balign 4
    .globl  spx_preproc_rvv_apply_gain_f32
    .hidden spx_preproc_rvv_apply_gain_f32
    .type   spx_preproc_rvv_apply_gain_f32, @function
spx_preproc_rvv_apply_gain_f32:
    addi    a0, a0, 4
    addi    a1, a1, 4
    addi    a2, a2, -1             /* N-1 pairs */
1:  vsetvli t0, a2, e32, m4, ta, ma
    vle32.v v4, (a0)               /* gain2 */
    vlseg2e32.v v8, (a1)           /* re v8, im v12 */
    vfmul.vv v8, v8, v4
    vfmul.vv v12, v12, v4
    vsseg2e32.v v8, (a1)
    sub     a2, a2, t0
    slli    t1, t0, 2
    add     a0, a0, t1
    slli    t1, t0, 3
    add     a1, a1, t1
    bnez    a2, 1b
    ret
    .size spx_preproc_rvv_apply_gain_f32, . - spx_preproc_rvv_apply_gain_f32

/* ------------------------------------------------------------------ */
/* void spx_preproc_rvv_word2int_sum_f32(spx_int16_t *x, const float *a,
 *                                       const float *b, int len)
 * x[i] = WORD2INT(a[i] + b[i]) (the overlap-add output loop), bit for
 * bit vs arch.h's float WORD2INT without the per-sample libm floor().
 * The a+b add must round like the C '+', so it runs under the caller's
 * rounding mode; only the +0.5 and the floor-convert run under RDN
 * (rounding x+0.5 down never crosses an integer boundary, so RDN there
 * is exact -- equivalence verified exhaustively over all floats), hence
 * the frm save/restore inside the strip loop. The clamp to
 * [-32768.0, 32766.5] collapses WORD2INT's two compare branches and
 * keeps the narrowing convert inside int16 range (never saturates). */
    .balign 4
    .globl  spx_preproc_rvv_word2int_sum_f32
    .hidden spx_preproc_rvv_word2int_sum_f32
    .type   spx_preproc_rvv_word2int_sum_f32, @function
spx_preproc_rvv_word2int_sum_f32:
    li      t1, 0x46FFFD00         /* 32766.5f */
    fmv.w.x fa0, t1
    li      t1, 0xC7000000         /* -32768.0f */
    fmv.w.x fa1, t1
    li      t1, 0x3F000000         /* 0.5f */
    fmv.w.x fa2, t1
1:  vsetvli t0, a3, e32, m8, ta, ma
    vle32.v v8, (a1)
    vle32.v v16, (a2)
    vfadd.vv v8, v8, v16           /* caller's frm, matches the C '+' */
    vfmin.vf v8, v8, fa0
    vfmax.vf v8, v8, fa1
    fsrmi   t6, 2                  /* frm := RDN, save the caller's */
    vfadd.vf v8, v8, fa2
    vsetvli zero, zero, e16, m4, ta, ma
    vfncvt.x.f.w v16, v8           /* floor, exact int16 range */
    fsrm    t6
    vse16.v v16, (a0)
    sub     a3, a3, t0
    slli    t1, t0, 2
    add     a1, a1, t1
    add     a2, a2, t1
    slli    t1, t0, 1
    add     a0, a0, t1
    bnez    a3, 1b
    ret
    .size spx_preproc_rvv_word2int_sum_f32, . - spx_preproc_rvv_word2int_sum_f32

/* ------------------------------------------------------------------ */
/* hypergeom_gain interpolation table (preprocess.c), 21 floats. */
    .section .rodata
    .balign 4
.Lhypergeom_table:
    .float 0.82157, 1.02017, 1.20461, 1.37534, 1.53363, 1.68092, 1.81865
    .float 1.94811, 2.07038, 2.18638, 2.29688, 2.40255, 2.50391, 2.60144
    .float 2.69551, 2.78647, 2.87458, 2.96015, 3.04333, 3.12431, 3.20326

#if defined(__ELF__)
/* Non-executable stack note */
    .section .note.GNU-stack, "", @progbits
#endif
