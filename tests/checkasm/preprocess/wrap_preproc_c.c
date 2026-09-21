/* Pure-C reference build of the preprocess.c per-bin kernels. #undef
 * USE_RVV keeps the loops scalar; this TU is the benchmark baseline, so
 * it is built with the no-autovec flags (checkasm_c_ref_args in
 * tests/meson.build). */
#define CKA_PREFIX  ckapc_
#define CKA_VARIANT c
#include "wrap_preproc_rename.h"
#include "wrap.h"

#undef USE_RVV

#include "../wrap_impl.h"
#include "wrap_preproc_stubs.h"

#ifndef FIXED_POINT
#include "wrap_preproc_shims.h"
#endif
