/* External-linkage shim around filterbank.c's static inline fbank_psd16.
 * Included by wrap_fbank_c.c and wrap_fbank_rvv.c after the module, with
 * CKA_VARIANT naming the TU; wrap.h carries the matching prototype. */

void CKA_SHIM(fbank_psd16)(const FilterBank *bank, const float *mel, float *ps)
{
    fbank_psd16(bank->bank_left, bank->bank_right, bank->filter_left,
                bank->filter_right, mel, ps, bank->len);
}
