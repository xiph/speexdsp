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
| `inner_product_double` | Float | SSE2² | ✓ Optimized |
| `interpolate_product_single` | Float | SSE¹ | ✓ Optimized |
| `interpolate_product_double` | Float | SSE2² | ✓ Optimized |

¹ The `SSE` handles both `SSE` and `SSE2`<BR>
² No `SSE` for `_double`

## Do Not Optimize

- `FIXED_POINT` `SSE` not needed
- `FIXED_POINT` `_double` not needed

| Function | Data Type | Architecture | Status |
|----------|-----------|--------------|--------|
| `inner_product_single` | Fixed | SSE | Do not optimize |
| `inner_product_double` | Fixed | SSE | Do not optimize |
| `inner_product_double` | Fixed | NEON | Do not optimize |
| `interpolate_product_single` | Fixed | SSE | Do not optimize |
| `interpolate_product_double` | Fixed | SSE | Do not optimize |
| `interpolate_product_double` | Fixed | NEON | Do not optimize |
