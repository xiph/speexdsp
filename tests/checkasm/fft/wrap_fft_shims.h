/* External-linkage shims around kiss_fft.c's static kf_bfly* stages and
 * the whole transform, one shared signature. Included by wrap_fft_c.c and
 * wrap_fft_rvv.c after the module, with CKA_VARIANT naming the TU; wrap.h
 * carries the matching prototypes. The radix-3 shim replicates kf_work's
 * per-sub-FFT loop (its C butterfly takes no N/mm). */

void CKA_SHIM(kf_bfly2)(kiss_fft_cfg cfg, kiss_fft_cpx *Fout, int fstride, int m, int N, int mm)
{
    kf_bfly2(Fout, (size_t) fstride, cfg, m, N, mm);
}

void CKA_SHIM(kf_bfly3)(kiss_fft_cfg cfg, kiss_fft_cpx *Fout, int fstride, int m, int N, int mm)
{
    int i;
    for (i = 0; i < N; i++)
        kf_bfly3(Fout + i * mm, (size_t) fstride, cfg, m);
}

void CKA_SHIM(kf_bfly4)(kiss_fft_cfg cfg, kiss_fft_cpx *Fout, int fstride, int m, int N, int mm)
{
    kf_bfly4(Fout, (size_t) fstride, cfg, m, N, mm);
}

void CKA_SHIM(kf_bfly5)(kiss_fft_cfg cfg, kiss_fft_cpx *Fout, int fstride, int m, int N, int mm)
{
    kf_bfly5(Fout, (size_t) fstride, cfg, m, N, mm);
}

void CKA_SHIM(fft)(kiss_fft_cfg cfg, const kiss_fft_cpx *fin, kiss_fft_cpx *fout)
{
    kiss_fft(cfg, fin, fout);
}
