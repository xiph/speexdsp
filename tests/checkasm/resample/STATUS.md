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

---

## Test harness coverage

### `resample_kernels` — per-kernel unit test + micro-benchmark

The four kernels (`resampler_basic_*`) are tested over multiple (**rate, 
quality**) configurations. The configs are chosen so each kernel lands in the 
expected `(method, precision)` (asserted at runtime; mismatches are skipped 
with a note) while producing **different filter lengths**. 
Randomized input is compared against the C reference with a
relative tolerance (float) or 1-LSB tolerance (fixed point).  
Because `filt_len` is always a multiple of 8 but the aarch64 NEON
`inner_product_single` strides 16 with a `len % 16` remainder loop, the
direct-single configs deliberately include `filt_len % 16 == 8` cases (the 3:2
downsamples → 72, 120, and the Q0 upsample → 8, which is all-remainder) as well
as `% 16 == 0` cases (8000/16000→…, 80, 160), exercising both the main and
remainder paths.  
The SSE/SSE2 kernels instead stride 8 (or 2 for the
interpolating variants), and `filt_len` is always a multiple of 8, so they have
no remainder path — on x86 the same configs just exercise a range of
accumulation lengths.  

### `resample_process` — full-pipeline integration benchmark

Benchmarks the public resampler API end to end across a table of conversions 
spanning cost classes — cheap near-integer ratios (small den_rate, direct sinc 
table) vs. expensive coprime ratios (large den_rate, interpolation, e.g. 
44100↔48000 = 147:160) — plus a quality sweep. Each conversion is run through 
both I/O APIs: speex_resampler_process_float (rows resample_*) and speex_resampler_process_int (rows resample_int_*). The int16 path additionally 
covers the NEON WORD2INT override (saturate_float_to_16bit, and 
saturate_32bit_to_16bit in fixed point), which the float path and the kernel 
tests never reach. For each, the C and SIMD pipelines are timed side by side 
and verified against each other end to end. The SIMD build is whichever the 
platform targets: on aarch64 it is NEON (a speedup on the direct-single 
conversions, matching C on the interpolating/double ones, which have no 
override yet); on x86 it is SSE/SSE2 (the SSE translation unit carries both, 
so all four kernels are accelerated).

## CI

`checkasm` is registered as a Meson `test()`, so every `.meson` CI job
(`fixed-point`, `no-float`, `fftw3`, `no-simd`, …) runs the *correctness* pass
through `meson test` — a pass/fail gate that fails the pipeline on any
C-vs-SIMD mismatch. The dedicated `meson tests-checkasm` job additionally runs
`--bench` for performance numbers.
