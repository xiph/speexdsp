/* RVV build of the kf_bfly* butterflies for checkasm: base-ISA C, the V
 * instructions live in kiss_fft_rvv_asm.S (linked alongside).
 * KISS_FFT_RVV_FORCE_ON pins the dispatch on so the asm is tested
 * deterministically, not per the build host's getauxval. */
#define CKA_PREFIX  ckafrvv_
#define CKA_VARIANT rvv
#include "wrap_fft_rename.h"
#include "wrap.h"

#ifdef USE_RVV
#  define KISS_FFT_RVV_FORCE_ON
#endif

#include "../wrap_impl.h"

#ifdef USE_RVV
#include "wrap_fft_shims.h"
#endif
