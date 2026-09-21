/* SSE build of the single-precision resampler functions. Compiled only when
 * has_sse && !fixed-point. Keeps the native USE_SSE so resample.c pulls in
 * resample_sse.h, overriding inner_product_single / interpolate_product_single.
 * The TU also keeps USE_SSE2, so its full-pipeline state runs SSE single- and
 * SSE2 double-precision kernels, i.e. the full library pipeline. See
 * ../wrap_impl.h for the HAVE_CONFIG_H/rename mechanics. */
#define CKA_PREFIX  ckasse_
#define CKA_VARIANT sse
#include "wrap_resample_rename.h"
#include "wrap.h"
#include "../wrap_impl.h"
#include "wrap_resample_shims.h"

CKA_RESAMPLER_BASIC_SHIM(direct_single)
CKA_RESAMPLER_BASIC_SHIM(interpolate_single)
