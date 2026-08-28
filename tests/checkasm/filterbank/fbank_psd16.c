#include <stdio.h>
#include <string.h>
#include "../internal.h"
#include "wrap.h"

/* filterbank_compute_psd16: unit test + micro-benchmark, C vs RVV
 * (float builds only -- the RVV kernel has no fixed-point twin). The
 * bank tables come from the real filterbank_new at the preprocessor's
 * M=24 bands, over the spectrum sizes preprocess.c produces for the
 * common rates plus a sub-vector-length size (5) and ragged tails
 * (33/257). The aliased layout preprocess.c uses (ps and mel disjoint
 * halves of one buffer, gain2 vs gain2+N) is checked too. */

#if !defined(FIXED_POINT) && defined(HAVE_RVV_FBANK)

static const struct { int rate, len; } configs[] = {
    {  8000,   5 },
    {  8000,   8 },
    {  8000,  16 },
    {  8000,  33 },
    {  8000, 128 },
    {  8000, 160 },
    { 16000, 257 },
    { 16000, 320 },
    { 32000, 640 },
    { 48000, 960 },
};
#define NUM_CONFIGS (sizeof(configs) / sizeof(configs[0]))

enum { MAX_LEN = 960, NBANDS = 24 };

void checkasm_check_fbank_psd16(void)
{
    checkasm_declare(void, const FilterBank *, const float *, float *);

    CHECKASM_ALIGN(float mel[NBANDS]);
    CHECKASM_ALIGN(float ps_ref[MAX_LEN]);
    CHECKASM_ALIGN(float ps_new[MAX_LEN]);
    /* preprocess.c layout: ps = gain2, mel = gain2 + N */
    CHECKASM_ALIGN(float aliased_ref[MAX_LEN + NBANDS]);
    CHECKASM_ALIGN(float aliased_new[MAX_LEN + NBANDS]);

    for (size_t c = 0; c < NUM_CONFIGS; c++) {
        const int rate = configs[c].rate;
        const int len = configs[c].len;
        FilterBank *bank = fbank_new_c(NBANDS, rate, len, 1);

        fbank_fill_gain(mel, NBANDS);

        if (checkasm_check_func(fbank_psd16_c, "fbank_psd16_%d_n%d", rate, len))
            checkasm_bench_new(bank, mel, ps_new);

        if (active_flags & SPEEXDSP_CPU_FLAG_RVV) {
            if (checkasm_check_func(fbank_psd16_rvv, "fbank_psd16_%d_n%d", rate, len)) {
                checkasm_call_ref(bank, mel, ps_ref);
                checkasm_call_new(bank, mel, ps_new);
                if (!fbank_buf_within_tol(ps_ref, ps_new, len, FBANK_PSD16_F32_REL_TOL))
                    checkasm_fail();

                fbank_fill_gain(aliased_ref + len, NBANDS);
                memcpy(aliased_new + len, aliased_ref + len, NBANDS * sizeof *aliased_new);
                checkasm_call_ref(bank, aliased_ref + len, aliased_ref);
                checkasm_call_new(bank, aliased_new + len, aliased_new);
                if (!fbank_buf_within_tol(aliased_ref, aliased_new, len, FBANK_PSD16_F32_REL_TOL))
                    checkasm_fail();

                checkasm_bench_new(bank, mel, ps_new);
            }
        }

        fbank_destroy_c(bank);
    }

    checkasm_report("fbank_psd16");
}

#else /* FIXED_POINT || !HAVE_RVV_FBANK */

void checkasm_check_fbank_psd16(void)
{
}

#endif
