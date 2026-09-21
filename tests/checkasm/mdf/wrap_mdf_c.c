/* Pure-C reference build of the mdf.c spectral kernels. #undef USE_RVV
 * keeps the loops scalar; this TU is the benchmark baseline, so it is
 * built with the no-autovec flags (checkasm_c_ref_args in
 * tests/meson.build). */
#define CKA_PREFIX  ckamc_
#define CKA_VARIANT c
#include "wrap_mdf_rename.h"
#include "wrap.h"

#undef USE_RVV

#include "../wrap_impl.h"
#include "wrap_mdf_stubs.h"
#include "wrap_mdf_shims.h"
