/* Pure-C reference build of the kf_bfly* butterflies plus the shared cfg and
 * stage helpers. #undef USE_RVV keeps the butterflies scalar; this TU is the
 * benchmark baseline, so it is built with the no-autovec flags
 * (checkasm_c_ref_args in tests/meson.build). */
#define CKA_PREFIX  ckafc_
#define CKA_VARIANT c
#include "wrap_fft_rename.h"
#include "wrap.h"

#undef USE_RVV

#include "../wrap_impl.h"
#include "wrap_fft_shims.h"

/* ------------- cfg + stage helpers ------------- */

kiss_fft_cfg fft_make_cfg(int nfft, int inverse)
{
    return kiss_fft_alloc(nfft, inverse, NULL, NULL);
}

void fft_free_cfg(kiss_fft_cfg cfg)
{
    speex_free(cfg);
}

/* Mirrors kf_work's traversal: level i has fstride = N = product of the
 * earlier radices and mm = the parent stage length (1 for the outermost
 * level, where N==1 leaves it unused). */
int fft_stages(int nfft, struct fft_stage *out, int max)
{
    int facbuf[2 * MAXFACTORS];
    int n = 0, fstride = 1, prev_m = 1, i;

    kf_factor(nfft, facbuf);
    for (i = 0; n < max; i++) {
        int p = facbuf[2 * i], m = facbuf[2 * i + 1];
        out[n].p       = p;
        out[n].m       = m;
        out[n].fstride = fstride;
        out[n].N       = fstride;
        out[n].mm      = (i == 0) ? 1 : prev_m;
        n++;
        if (m == 1)
            break;
        prev_m = m;
        fstride *= p;
    }
    return n;
}
