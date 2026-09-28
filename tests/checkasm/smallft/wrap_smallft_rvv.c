/* RVV build of the smallft radix stages for checkasm: base-ISA C, the V
 * instructions live in smallft_rvv_asm.S (linked alongside).
 * SMALLFT_RVV_FORCE_ON pins the dispatch on so the asm is tested
 * without reading the state's arch. checkasm still runs this variant only
 * when spx_rvv_detect() passes. */
#define CKA_PREFIX  ckasrvv_
#define CKA_VARIANT rvv
#include "wrap_smallft_rename.h"
#include "wrap.h"

#ifdef USE_RVV
#  define SMALLFT_RVV_FORCE_ON
#endif

#include "../wrap_impl.h"

#ifdef HAVE_RVV_SMALLFT
#include "wrap_smallft_shims.h"
#endif
