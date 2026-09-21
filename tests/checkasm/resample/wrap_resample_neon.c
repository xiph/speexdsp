/* NEON build of the resampler functions. Compiled only when has_neon. Keeps the
 * native USE_NEON so resample.c pulls in resample_neon.h. Only the kernels that
 * header overrides (the HAVE_NEON_* gates in wrap.h) differ from the C
 * reference and are exposed; the double-precision ones would be byte-identical
 * to C. See ../wrap_impl.h for the HAVE_CONFIG_H/rename mechanics. */
#define CKA_PREFIX  ckaneon_
#define CKA_VARIANT neon
#include "wrap_resample_rename.h"
#include "wrap.h"

#ifdef USE_NEON
#  if !defined(__aarch64__) && !defined(__ARM_NEON)
#    error "wrap_resample_neon.c built without aarch64 / ARM NEON ISA; check toolchain flags"
#  endif
#endif

#include "../wrap_impl.h"
#include "wrap_resample_shims.h"

#ifdef HAVE_NEON_DIRECT_SINGLE
CKA_RESAMPLER_BASIC_SHIM(direct_single)
#endif
#ifdef HAVE_NEON_INTERPOLATE_SINGLE
CKA_RESAMPLER_BASIC_SHIM(interpolate_single)
#endif
