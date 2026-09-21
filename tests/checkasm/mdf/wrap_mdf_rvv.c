/* RVV build of the mdf.c spectral kernels for checkasm: base-ISA C, the V
 * instructions live in mdf_rvv_asm.S (linked alongside). MDF_RVV_FORCE_ON
 * pins the dispatch on so the asm is tested deterministically, not per the
 * build host's getauxval. */
#define CKA_PREFIX  ckamrvv_
#define CKA_VARIANT rvv
#include "wrap_mdf_rename.h"
#include "wrap.h"

#ifdef USE_RVV
#  define MDF_RVV_FORCE_ON
#endif

#include "../wrap_impl.h"
#include "wrap_mdf_stubs.h"

#ifdef HAVE_RVV_MDF
#include "wrap_mdf_shims.h"
#endif
