/* Copyright (C) 2026 Tristan Matthews */
/*
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
