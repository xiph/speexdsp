/* Pure-C reference build of the filterbank psd16 loop. #undef USE_RVV
 * keeps the loop scalar; this TU is the benchmark baseline, so it is
 * built with the no-autovec flags (checkasm_c_ref_args in
 * tests/meson.build). It also hosts the real (renamed) filterbank_new
 * the test uses to build authentic bank tables. */
#define CKA_PREFIX  ckafbc_
#define CKA_VARIANT c
#include "wrap_fbank_rename.h"
#include "wrap.h"

#undef USE_RVV

#include "../wrap_impl.h"

#ifndef FIXED_POINT

#include "wrap_fbank_shims.h"

FilterBank *fbank_new_c(int banks, float sampling, int len, int type)
{
    return filterbank_new(banks, sampling, len, type);
}

void fbank_destroy_c(FilterBank *bank)
{
    filterbank_destroy(bank);
}

#endif /* !FIXED_POINT */
