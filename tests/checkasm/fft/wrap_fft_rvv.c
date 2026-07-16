/* RVV build of the kf_bfly* butterflies for checkasm: base-ISA C, the V
 * instructions live in kiss_fft_rvv_asm.S (linked alongside).
 * KISS_FFT_RVV_FORCE_ON pins the dispatch on so the asm is tested
 * deterministically, not per the build host's getauxval. */
#define CKA_PREFIX ckafrvv_
#include "wrap_fft_rename.h"
#include "wrap.h"

#ifdef USE_RVV
#  define KISS_FFT_RVV_FORCE_ON
#endif

#include "wrap_fft_impl.h"

#ifdef USE_RVV

void kf_bfly2_rvv(kiss_fft_cfg cfg, kiss_fft_cpx *Fout, int fstride, int m, int N, int mm)
{
    kf_bfly2(Fout, (size_t) fstride, cfg, m, N, mm);
}

void kf_bfly3_rvv(kiss_fft_cfg cfg, kiss_fft_cpx *Fout, int fstride, int m, int N, int mm)
{
    int i;
    for (i = 0; i < N; i++)
        kf_bfly3(Fout + i * mm, (size_t) fstride, cfg, m);
}

void kf_bfly4_rvv(kiss_fft_cfg cfg, kiss_fft_cpx *Fout, int fstride, int m, int N, int mm)
{
    kf_bfly4(Fout, (size_t) fstride, cfg, m, N, mm);
}

void kf_bfly5_rvv(kiss_fft_cfg cfg, kiss_fft_cpx *Fout, int fstride, int m, int N, int mm)
{
    kf_bfly5(Fout, (size_t) fstride, cfg, m, N, mm);
}

void fft_rvv(kiss_fft_cfg cfg, const kiss_fft_cpx *fin, kiss_fft_cpx *fout)
{
    kiss_fft(cfg, fin, fout);
}

#endif /* USE_RVV */
