/* RVV build of the filterbank psd16 loop for checkasm: base-ISA C, the
 * V instructions live in filterbank_rvv_asm.S (linked alongside).
 * FBANK_RVV_FORCE_ON pins the dispatch on so the asm is tested
 * without reading the state's arch. checkasm still runs this variant only
 * when spx_rvv_detect() passes. */
#define CKA_PREFIX  ckafbrvv_
#define CKA_VARIANT rvv
#include "wrap_fbank_rename.h"
#include "wrap.h"

#ifdef USE_RVV
#  define FBANK_RVV_FORCE_ON
#endif

#include "../wrap_impl.h"

#ifdef HAVE_RVV_FBANK
#include "wrap_fbank_shims.h"
#endif
