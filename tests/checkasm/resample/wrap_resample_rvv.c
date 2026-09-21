/* RVV build of the resampler kernels for checkasm. Base-ISA C; the V
 * instructions live only in resample_rvv_asm.S (linked in alongside). Keeps the
 * native USE_RVV so resample.c pulls in resample_rvv.h.
 *
 * The library dispatches at runtime via spx_rvv_enabled; the harness instead
 * sets RESAMPLE_RVV_FORCE_ON (=> SPX_RVV_ON == 1) so checkasm deterministically
 * tests the asm kernel, not whatever the build host's getauxval reports.
 * See ../wrap_impl.h for the rename mechanics. */
#define CKA_PREFIX  ckarvv_
#define CKA_VARIANT rvv
#include "wrap_resample_rename.h"
#include "wrap.h"

#ifdef USE_RVV
#  define RESAMPLE_RVV_FORCE_ON
#endif

#include "../wrap_impl.h"

#ifdef USE_RVV
#include "wrap_resample_shims.h"

CKA_RESAMPLER_BASIC_SHIM(direct_single)
CKA_RESAMPLER_BASIC_SHIM(interpolate_single)
#ifndef FIXED_POINT
CKA_RESAMPLER_BASIC_SHIM(direct_double)
CKA_RESAMPLER_BASIC_SHIM(interpolate_double)
#endif
#endif /* USE_RVV */
