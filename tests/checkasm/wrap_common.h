#ifndef SPEEXDSP_TESTS_CHECKASM_WRAP_COMMON_H
#define SPEEXDSP_TESTS_CHECKASM_WRAP_COMMON_H

/* Shared plumbing for the wrap_<mod>_<variant>.c translation units. Each
 * one #includes a whole libspeexdsp module (see wrap_impl.h) to get a C or
 * SIMD build of the module's static kernels, so the module's public API is
 * renamed to a per-TU prefix (the module's wrap_<mod>_rename.h) and the
 * shims exposing the kernels carry a per-TU suffix (wrap_<mod>_shims.h).
 *
 *   #define CKA_PREFIX  ckamc_        unique token, before the rename header
 *   #define CKA_VARIANT c             c / sse / neon / rvv, before the shims
 */

#ifndef CKA_PREFIX
#  error "define CKA_PREFIX (a unique token) before including the wrap headers"
#endif

#define CKA_CAT2(a, b) a ## b
#define CKA_CAT(a, b)  CKA_CAT2(a, b)

/* CKA_SHIM(foo) -> foo_c, foo_rvv, ... per CKA_VARIANT. The name is pasted
 * in the outer macro so a kernel name that the module's RVV header has
 * #define'd to its dispatch wrapper is not expanded first. */
#define CKA_SHIM(name)    CKA_SHIM_1(name ## _, CKA_VARIANT)
#define CKA_SHIM_1(a, b)  CKA_CAT2(a, b)

#endif /* SPEEXDSP_TESTS_CHECKASM_WRAP_COMMON_H */
