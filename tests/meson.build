if not get_option('tests').enabled()
    subdir_done()
endif

checkasm_sources = files(
    'checkasm/checkasm.c',
    'checkasm/cpu.c',
    'checkasm/resample/resample_kernels.c',
    'checkasm/resample/resampler_basic_direct_single.c',
    'checkasm/resample/resampler_basic_direct_double.c',
    'checkasm/resample/resampler_basic_interpolate_single.c',
    'checkasm/resample/resampler_basic_interpolate_double.c',
    'checkasm/resample/resample_process.c',
    'checkasm/fft/fft_bfly.c',
    'checkasm/fft/fft_transform.c',
)

# The C-reference TU is the scalar baseline every benchmark compares against, so
# disable auto-vectorization for it. If the toolchain's default -march carries V
# (vendor toolchains, -Dc_args=-march=rv64gcv, zig cross files with +v), the
# autovectorizer would turn the baseline into vector code and the benchmark would
# silently compare RVV against RVV. -mno-implicit-float additionally blocks the
# vector code clang leaks via aggregate-zeroing lowering (it emits vse64.v for the
# accum[4] zeroing even with the -fno-*vectorize flags set);
# --param=riscv-autovec-preference=none is the gcc-side autovec switch, and
# has_argument() drops whichever flag the compiler rejects.
checkasm_c_ref_args = []
if host_cpu_family.startswith('riscv')
    foreach arg : ['-fno-tree-vectorize', '-fno-tree-slp-vectorize',
                   '-fno-vectorize', '-fno-slp-vectorize',
                   '-mno-implicit-float',
                   '--param=riscv-autovec-preference=none']
        if cc.has_argument(arg)
            checkasm_c_ref_args += arg
        endif
    endforeach
endif

if has_neon
    checkasm_sources += files('checkasm/resample/wrap_resample_neon.c')
endif
if has_sse and not get_option('fixed-point').enabled()
    checkasm_sources += files('checkasm/resample/wrap_resample_sse.c')
endif
if has_sse2 and not get_option('fixed-point').enabled()
    checkasm_sources += files('checkasm/resample/wrap_resample_sse2.c')
endif

checkasm_incs = [
    include_directories('../libspeexdsp'),
    config_inc,
    speexdsp_incs,
]

# Each wrapper #includes a renamed full copy of resample.c, so each needs its own
# static_library: the C-ref lib carries the no-autovec c_args; the RVV lib is
# base-ISA C and links resample_rvv_asm.S to supply the kernels.
checkasm_extra_libs = [
    static_library('checkasm_resample_c',
        files('checkasm/resample/wrap_resample_c.c'),
        speexdsp_config_types_h,
        include_directories: checkasm_incs,
        c_args: checkasm_c_ref_args,
        dependencies: [checkasm_dep],
    ),
    static_library('checkasm_fft_c',
        files('checkasm/fft/wrap_fft_c.c'),
        speexdsp_config_types_h,
        include_directories: checkasm_incs,
        c_args: checkasm_c_ref_args,
        dependencies: [checkasm_dep],
    ),
]
if has_rvv_asm
    checkasm_extra_libs += static_library('checkasm_resample_rvv',
        files('checkasm/resample/wrap_resample_rvv.c',
              '../libspeexdsp/resample_rvv_asm.S'),
        speexdsp_config_types_h,
        include_directories: checkasm_incs,
        dependencies: [checkasm_dep],
    )
    checkasm_extra_libs += static_library('checkasm_fft_rvv',
        files('checkasm/fft/wrap_fft_rvv.c',
              '../libspeexdsp/kiss_fft_rvv_asm.S'),
        speexdsp_config_types_h,
        include_directories: checkasm_incs,
        dependencies: [checkasm_dep],
    )
endif

checkasm = executable('checkasm',
    checkasm_sources,
    speexdsp_config_types_h,
    include_directories: checkasm_incs,
    link_with: checkasm_extra_libs,
    dependencies: [
        checkasm_dep,
    ],
)

# SIMD Correctness gate: meson test -C build
test('checkasm', checkasm, timeout: 300)

# SIMD Benchmark: meson test --benchmark -C build
benchmark('checkasm', checkasm, timeout: 0, verbose: true, 
    args: ['--bench', '--duration=30000', '--repeat=5'])
