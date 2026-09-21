/* RVV build of the preprocess.c per-bin kernels for checkasm: base-ISA C,
 * the V instructions live in preprocess_rvv_asm.S (linked alongside).
 * PREPROC_RVV_FORCE_ON pins the dispatch on so the asm is tested
 * deterministically, not per the build host's getauxval. */
#define CKA_PREFIX  ckaprvv_
#define CKA_VARIANT rvv
#include "wrap_preproc_rename.h"
#include "wrap.h"

#ifdef USE_RVV
#  define PREPROC_RVV_FORCE_ON
#endif

#include "../wrap_impl.h"
#include "wrap_preproc_stubs.h"

#ifdef HAVE_RVV_PREPROC
#include "wrap_preproc_shims.h"
#endif
