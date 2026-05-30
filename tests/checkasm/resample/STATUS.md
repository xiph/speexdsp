# speexdsp SIMD Optimization Status

## NEON

| Function | Data Type | Status |
|----------|-----------|--------|
| `inner_product_single` | Fixed | ✓ Optimized |
| `inner_product_single` | Float | ✓ Optimized |
| `inner_product_double` | Float | ✗ No optimization |
| `interpolate_product_single` | Fixed | ✗ No optimization |
| `interpolate_product_single` | Float | ✗ No optimization |
| `interpolate_product_double` | Float | ✗ No optimization |

## SSE/SSE2

| Function | Data Type | Architecture | Status |
|----------|-----------|--------------|--------|
| `inner_product_single` | Float | SSE¹ | ✓ Optimized |
| `interpolate_product_single` | Float | SSE¹ | ✓ Optimized |
| `inner_product_double` | Float | SSE2 | ✓ Optimized |
| `interpolate_product_double` | Float | SSE2 | ✓ Optimized |

¹ The `_single` routines use plain SSE intrinsics (`<xmmintrin.h>`,
`resample_sse.h`) with no `#ifdef USE_SSE2` guard, so there is no distinct
SSE2 implementation to test — the same code simply also runs on SSE2-capable
CPUs. The checkasm test for `inner_product_single` therefore exercises the SSE
variant only (gated on `SPEEXDSP_CPU_FLAG_SSE`). Only the `_double` routines
genuinely require SSE2 (`<emmintrin.h>`, double-precision accumulation).

## Do Not Optimize

- FIXED_POINT SSE: Not needed.
- FIXED_POINT NEON: Only `_single` functions used; `_double` variants not applicable

| Function | Data Type | Architecture | Status |
|----------|-----------|--------------|--------|
| `inner_product_double` | Float | NEON | Do not optimize |
| `interpolate_product_double` | Float | NEON | Do not optimize |
| `inner_product_single` | Fixed | SSE | Do not optimize |
| `interpolate_product_single` | Fixed | SSE | Do not optimize |
| `inner_product_double` | Fixed | SSE | Do not optimize |
| `interpolate_product_double` | Fixed | SSE | Do not optimize |