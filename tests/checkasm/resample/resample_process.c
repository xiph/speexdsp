#include <stdio.h>
#include "../internal.h"
#include "wrap.h"

/* Integration benchmark + correctness: time the *full* resampler pipeline across
 * a range of sample-rate conversions, and verify the SIMD pipeline against the
 * pure-C pipeline end to end. Two I/O paths are exercised per conversion:
 *
 *   - speex_resampler_process_float  (float in/out)
 *   - speex_resampler_process_int    (int16 in/out) -- additionally covers the
 *       NEON WORD2INT override (saturate_float_to_16bit), which the float path
 *       and the kernel tests never touch.
 *
 * "Certain sampling-rate ratios are more expensive than others": the cost is
 * driven by the reduced ratio num:den. A near-integer ratio reduces to a small
 * den_rate and uses a direct per-phase sinc table (cheap); a coprime ratio such
 * as 44100<->48000 (147:160) needs a large den_rate, forcing the interpolating
 * kernel (expensive). Each conversion is benchmarked separately (the name
 * encodes in_rate/out_rate/quality), C and SIMD side by side.
 *
 * The SIMD pipeline is whichever the build targets: NEON on aarch64, SSE/SSE2 on
 * x86. NEON only overrides inner_product_single (+ WORD2INT), so it speeds up the
 * small-den_rate (direct-single) conversions and matches C on the rest; SSE/SSE2
 * override all four kernels. The C-vs-SIMD comparison passes either way. */

#ifndef DISABLE_FLOAT_API

/* OUT_LEN / IN_LEN are PER CHANNEL. The interleaved (multi-channel) rows reuse
 * the same per-channel lengths, so the buffers are sized for MAX_CH channels:
 * IN_BUF = MAX_CH * IN_LEN, OUT_BUF = MAX_CH * OUT_LEN. */
enum { OUT_LEN = 256, IN_LEN = 2048, MAX_CH = 2, IN_BUF = MAX_CH * IN_LEN, OUT_BUF = MAX_CH * OUT_LEN };

static const struct { unsigned in_rate, out_rate; int quality; } conversions[] = {
    /* Cheap: ratio reduces to a small den_rate -> direct sinc table. */
    { 48000, 24000,  5 },   /* 2:1 downsample           */
    { 96000, 48000,  5 },   /* 2:1 downsample           */
    {  8000, 16000,  5 },   /* 1:2 upsample             */
    { 48000, 32000,  5 },   /* 3:2 downsample           */
    /* Expensive: coprime ratio -> large den_rate / interpolation. */
    { 44100, 48000,  5 },   /* 147:160 upsample         */
    { 48000, 44100,  5 },   /* 160:147 downsample       */
    { 44100, 32000,  5 },   /* 441:320 downsample       */
    { 16000, 44100,  5 },   /* 160:441 upsample         */
    /* Quality sweep on the canonical CD<->DAT conversion. */
    { 44100, 48000,  3 },
    { 44100, 48000,  8 },
    { 44100, 48000, 10 },
};

/* Stereo (interleaved) coverage: exercises the per-channel loop and, for int,
 * the deinterleave/reinterleave that the single-channel rows above skip. A small
 * representative spread -- one cheap direct-single ratio, one expensive
 * interpolating ratio, and a quality-10 (double-precision) variant. */
static const struct { unsigned in_rate, out_rate; int quality; } conversions_mc[] = {
    { 48000, 24000,  5 },   /* 2:1 down,    direct-single      */
    { 44100, 48000,  5 },   /* 147:160 up,  interpolate-single */
    { 44100, 48000, 10 },   /* 147:160 up,  interpolate-double */
};

/* Select this architecture's SIMD build (NEON and SSE are mutually exclusive).
 * The SSE translation unit carries both USE_SSE and USE_SSE2 on x86-64, so its
 * state covers single- and double-precision conversions alike. */
#if defined(USE_NEON)
#  define SIMD_CPU_FLAG                 SPEEXDSP_CPU_FLAG_NEON
#  define resample_make_state_simd      resample_make_state_neon
#  define resample_make_state_simd_ch   resample_make_state_neon_ch
#  define resample_process_simd         resample_process_neon
#  define resample_process_int_simd     resample_process_int_neon
#  define resample_process_il_simd      resample_process_il_neon
#  define resample_process_int_il_simd  resample_process_int_il_neon
#elif defined(USE_SSE) && !defined(FIXED_POINT)
#  define SIMD_CPU_FLAG                 SPEEXDSP_CPU_FLAG_SSE
#  define resample_make_state_simd      resample_make_state_sse
#  define resample_make_state_simd_ch   resample_make_state_sse_ch
#  define resample_process_simd         resample_process_sse
#  define resample_process_int_simd     resample_process_int_sse
#  define resample_process_il_simd      resample_process_il_sse
#  define resample_process_int_il_simd  resample_process_int_il_sse
#endif

void checkasm_check_resample_process(void)
{
    CHECKASM_ALIGN(float fin[IN_BUF]);
    CHECKASM_ALIGN(float fout_ref[OUT_BUF]);
    CHECKASM_ALIGN(float fout_new[OUT_BUF]);
    CHECKASM_ALIGN(spx_int16_t iin[IN_BUF]);
    CHECKASM_ALIGN(spx_int16_t iout_ref[OUT_BUF]);
    CHECKASM_ALIGN(spx_int16_t iout_new[OUT_BUF]);
    resample_fill_float(fin, IN_BUF);
    resample_fill_int16(iin, IN_BUF);

    for (size_t i = 0; i < sizeof(conversions) / sizeof(conversions[0]); i++) {
        const unsigned ir = conversions[i].in_rate;
        const unsigned orr = conversions[i].out_rate;
        const int q = conversions[i].quality;

        SpeexResamplerState *st_c = resample_make_state(ir, orr, q);
        if (!st_c) {
            fprintf(stderr, "resample_process: init failed %u->%u q%d\n", ir, orr, q);
            continue;
        }
#ifdef SIMD_CPU_FLAG
        SpeexResamplerState *st_s =
            (active_flags & SIMD_CPU_FLAG) ? resample_make_state_simd(ir, orr, q) : NULL;
#endif

        /* ---------------- float pipeline ---------------- */
        {
            checkasm_declare(int, SpeexResamplerState *, const float *, spx_uint32_t,
                             float *, spx_uint32_t);
            if (checkasm_check_func(resample_process_c, "resample_%u_%u_q%d", ir, orr, q))
                checkasm_bench_new(st_c, fin, IN_LEN, fout_new, OUT_LEN);
#ifdef SIMD_CPU_FLAG
            if (st_s && checkasm_check_func(resample_process_simd, "resample_%u_%u_q%d", ir, orr, q)) {
                /* Each wrapper zeroes the state, so st_c is a valid reference
                 * even after being benchmarked above. */
                int nc = resample_process_c(st_c, fin, IN_LEN, fout_ref, OUT_LEN);
                int nn = resample_process_simd(st_s, fin, IN_LEN, fout_new, OUT_LEN);
                /* Gate self-test: -DCHECKASM_FORCE_FAIL injects a deterministic
                 * divergence so the acceptance run can prove the gate fails on a
                 * mismatch. Inert otherwise. */
#ifdef CHECKASM_FORCE_FAIL
                if (nn > 0) fout_new[0] = fout_ref[0] + 1000.0f;
#endif
                if (nc != nn ||
                    !resample_float_within_tol(fout_ref, fout_new, (unsigned) nc, RESAMPLE_PROCESS_REL_TOL))
                    checkasm_fail();
                checkasm_bench_new(st_s, fin, IN_LEN, fout_new, OUT_LEN);
            }
#endif
        }

        /* ---------------- int16 pipeline (also covers WORD2INT / saturate) ---------------- */
        {
            checkasm_declare(int, SpeexResamplerState *, const spx_int16_t *, spx_uint32_t,
                             spx_int16_t *, spx_uint32_t);
            if (checkasm_check_func(resample_process_int_c, "resample_int_%u_%u_q%d", ir, orr, q))
                checkasm_bench_new(st_c, iin, IN_LEN, iout_new, OUT_LEN);
#ifdef SIMD_CPU_FLAG
            if (st_s && checkasm_check_func(resample_process_int_simd, "resample_int_%u_%u_q%d", ir, orr, q)) {
                int nc = resample_process_int_c(st_c, iin, IN_LEN, iout_ref, OUT_LEN);
                int nn = resample_process_int_simd(st_s, iin, IN_LEN, iout_new, OUT_LEN);
                if (nc != nn ||
                    !resample_int16_within_lsb(iout_ref, iout_new, (unsigned) nc, RESAMPLE_PROCESS_INT_MAX_LSB))
                    checkasm_fail();
                checkasm_bench_new(st_s, iin, IN_LEN, iout_new, OUT_LEN);
            }
#endif
        }

#ifdef SIMD_CPU_FLAG
        if (st_s)
            resample_destroy_state(st_s);
#endif
        resample_destroy_state(st_c);
    }

    /* ================= Stereo (interleaved) coverage ================= *
     * Same per-channel IN_LEN/OUT_LEN as the mono rows, run through the
     * interleaved API at MAX_CH channels. in_len/out_len are per channel; the
     * comparison spans nc * MAX_CH interleaved output values. */
    for (size_t i = 0; i < sizeof(conversions_mc) / sizeof(conversions_mc[0]); i++) {
        const unsigned ir = conversions_mc[i].in_rate;
        const unsigned orr = conversions_mc[i].out_rate;
        const int q = conversions_mc[i].quality;

        SpeexResamplerState *st_c = resample_make_state_ch(ir, orr, q, MAX_CH);
        if (!st_c) {
            fprintf(stderr, "resample_process: init failed %uch %u->%u q%d\n", MAX_CH, ir, orr, q);
            continue;
        }
#ifdef SIMD_CPU_FLAG
        SpeexResamplerState *st_s =
            (active_flags & SIMD_CPU_FLAG) ? resample_make_state_simd_ch(ir, orr, q, MAX_CH) : NULL;
#endif

        /* ---------------- interleaved float pipeline ---------------- */
        {
            checkasm_declare(int, SpeexResamplerState *, const float *, spx_uint32_t,
                             float *, spx_uint32_t);
            if (checkasm_check_func(resample_process_il_c, "resample_%dch_%u_%u_q%d", MAX_CH, ir, orr, q))
                checkasm_bench_new(st_c, fin, IN_LEN, fout_new, OUT_LEN);
#ifdef SIMD_CPU_FLAG
            if (st_s && checkasm_check_func(resample_process_il_simd, "resample_%dch_%u_%u_q%d", MAX_CH, ir, orr, q)) {
                int nc = resample_process_il_c(st_c, fin, IN_LEN, fout_ref, OUT_LEN);
                int nn = resample_process_il_simd(st_s, fin, IN_LEN, fout_new, OUT_LEN);
                if (nc != nn ||
                    !resample_float_within_tol(fout_ref, fout_new, (unsigned) nc * MAX_CH, RESAMPLE_PROCESS_REL_TOL))
                    checkasm_fail();
                checkasm_bench_new(st_s, fin, IN_LEN, fout_new, OUT_LEN);
            }
#endif
        }

        /* ---------------- interleaved int16 pipeline (deinterleave + WORD2INT) ---------------- */
        {
            checkasm_declare(int, SpeexResamplerState *, const spx_int16_t *, spx_uint32_t,
                             spx_int16_t *, spx_uint32_t);
            if (checkasm_check_func(resample_process_int_il_c, "resample_int_%dch_%u_%u_q%d", MAX_CH, ir, orr, q))
                checkasm_bench_new(st_c, iin, IN_LEN, iout_new, OUT_LEN);
#ifdef SIMD_CPU_FLAG
            if (st_s && checkasm_check_func(resample_process_int_il_simd, "resample_int_%dch_%u_%u_q%d", MAX_CH, ir, orr, q)) {
                int nc = resample_process_int_il_c(st_c, iin, IN_LEN, iout_ref, OUT_LEN);
                int nn = resample_process_int_il_simd(st_s, iin, IN_LEN, iout_new, OUT_LEN);
                if (nc != nn ||
                    !resample_int16_within_lsb(iout_ref, iout_new, (unsigned) nc * MAX_CH, RESAMPLE_PROCESS_INT_MAX_LSB))
                    checkasm_fail();
                checkasm_bench_new(st_s, iin, IN_LEN, iout_new, OUT_LEN);
            }
#endif
        }

#ifdef SIMD_CPU_FLAG
        if (st_s)
            resample_destroy_state(st_s);
#endif
        resample_destroy_state(st_c);
    }

    checkasm_report("resample_process");
}

#else /* DISABLE_FLOAT_API */

void checkasm_check_resample_process(void) { }

#endif /* !DISABLE_FLOAT_API */
