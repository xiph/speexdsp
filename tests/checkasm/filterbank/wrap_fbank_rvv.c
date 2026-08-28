/* RVV build of the filterbank psd16 loop for checkasm: base-ISA C, the
 * V instructions live in filterbank_rvv_asm.S (linked alongside).
 * FBANK_RVV_FORCE_ON pins the dispatch on so the asm is tested
 * deterministically, not per the build host's getauxval. */
#define CKA_PREFIX ckafbrvv_
#include "wrap_fbank_rename.h"
#include "wrap.h"

#ifdef USE_RVV
#  define FBANK_RVV_FORCE_ON
#endif

#include "wrap_fbank_impl.h"

#ifdef HAVE_RVV_FBANK

void fbank_psd16_rvv(const FilterBank *bank, const float *mel, float *ps)
{
    fbank_psd16(bank->bank_left, bank->bank_right, bank->filter_left,
                bank->filter_right, mel, ps, bank->len);
}

#endif /* HAVE_RVV_FBANK */
