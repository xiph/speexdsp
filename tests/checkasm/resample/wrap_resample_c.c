/* Pure-C reference build of the four resampler_basic_* functions, plus the
 * shared state helpers. This is the only TU that owns SpeexResamplerState
 * construction/inspection (it has both the four static function addresses and
 * the public init API in scope after #including resample.c).
 *
 * We #undef the native USE_SSE/USE_SSE2/USE_NEON/USE_RVV before pulling in
 * resample.c so its #ifdef USE_SSE/USE_NEON/USE_RVV skip the SIMD headers and
 * the inner-product kernels stay the generic C fallbacks. ../wrap_impl.h then
 * handles the HAVE_CONFIG_H/rename mechanics and includes resample.c.
 *
 * This TU is the scalar baseline every benchmark compares against, so
 * auto-vectorization is disabled for it (checkasm_c_ref_args in
 * tests/meson.build) -- otherwise a toolchain whose default -march carries V
 * could vectorize the baseline and the benchmark would compare RVV against RVV. */
#define CKA_PREFIX  ckac_
#define CKA_VARIANT c
#include "wrap_resample_rename.h"
#include "wrap.h"

#undef USE_SSE
#undef USE_SSE2
#undef USE_NEON
#undef USE_RVV

#include "../wrap_impl.h"
#include "wrap_resample_shims.h"

/* ------------- State helpers ------------- */

void resample_destroy_state(SpeexResamplerState *st)
{
    speex_resampler_destroy(st);
}

enum resample_kind resample_kind(const SpeexResamplerState *st)
{
    if (st->resampler_ptr == resampler_basic_direct_single)
        return RESAMPLE_KIND_DIRECT_SINGLE;
    if (st->resampler_ptr == resampler_basic_interpolate_single)
        return RESAMPLE_KIND_INTERPOLATE_SINGLE;
#ifndef FIXED_POINT
    if (st->resampler_ptr == resampler_basic_direct_double)
        return RESAMPLE_KIND_DIRECT_DOUBLE;
    if (st->resampler_ptr == resampler_basic_interpolate_double)
        return RESAMPLE_KIND_INTERPOLATE_DOUBLE;
#endif
    return RESAMPLE_KIND_OTHER;
}

unsigned resample_filt_len(const SpeexResamplerState *st)   { return st->filt_len; }
unsigned resample_den_rate(const SpeexResamplerState *st)   { return st->den_rate; }
unsigned resample_oversample(const SpeexResamplerState *st) { return st->oversample; }

/* ------------- Functions under test ------------- */

CKA_RESAMPLER_BASIC_SHIM(direct_single)
CKA_RESAMPLER_BASIC_SHIM(interpolate_single)
#ifndef FIXED_POINT
CKA_RESAMPLER_BASIC_SHIM(direct_double)
CKA_RESAMPLER_BASIC_SHIM(interpolate_double)
#endif
