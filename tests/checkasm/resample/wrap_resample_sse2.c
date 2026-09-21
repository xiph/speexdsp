/* SSE2 build of the double-precision resampler functions. Compiled only when
 * has_sse2 && !fixed-point. The native USE_SSE2 makes resample_sse.h define
 * OVERRIDE_INNER_PRODUCT_DOUBLE / OVERRIDE_INTERPOLATE_PRODUCT_DOUBLE, so the
 * double kernels use the SSE2 intrinsics. No full-pipeline wrappers (the SSE
 * TU's state already runs the SSE2 double kernels), hence
 * CKA_KERNEL_SHIMS_ONLY. See ../wrap_impl.h for the HAVE_CONFIG_H/rename mechanics. */
#define CKA_PREFIX  ckasse2_
#define CKA_VARIANT sse2
#include "wrap_resample_rename.h"
#include "wrap.h"
#include "../wrap_impl.h"

#define CKA_KERNEL_SHIMS_ONLY
#include "wrap_resample_shims.h"

CKA_RESAMPLER_BASIC_SHIM(direct_double)
CKA_RESAMPLER_BASIC_SHIM(interpolate_double)
