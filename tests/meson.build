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
    'checkasm/mdf/mdf_kernels.c',
    'checkasm/preprocess/preproc_kernels.c',
    'checkasm/smallft/smallft_drft.c',
    'checkasm/fftwrap/fftwrap_scale.c',
    'checkasm/filterbank/fbank_psd16.c',
)

checkasm_c_ref_args = []
if get_option('checkasm-no-autovec').enabled()
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
# static_library: the C-ref lib carries checkasm_c_ref_args (scalar-forcing flags
# only under -Dcheckasm-no-autovec=enabled); the RVV lib is
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
    static_library('checkasm_mdf_c',
        files('checkasm/mdf/wrap_mdf_c.c'),
        speexdsp_config_types_h,
        include_directories: checkasm_incs,
        c_args: checkasm_c_ref_args,
        dependencies: [checkasm_dep],
    ),
    static_library('checkasm_preproc_c',
        files('checkasm/preprocess/wrap_preproc_c.c'),
        speexdsp_config_types_h,
        include_directories: checkasm_incs,
        c_args: checkasm_c_ref_args,
        dependencies: [checkasm_dep],
    ),
    static_library('checkasm_smallft_c',
        files('checkasm/smallft/wrap_smallft_c.c'),
        speexdsp_config_types_h,
        include_directories: checkasm_incs,
        c_args: checkasm_c_ref_args,
        dependencies: [checkasm_dep],
    ),
    static_library('checkasm_fbank_c',
        files('checkasm/filterbank/wrap_fbank_c.c'),
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
    checkasm_extra_libs += static_library('checkasm_mdf_rvv',
        files('checkasm/mdf/wrap_mdf_rvv.c',
              '../libspeexdsp/mdf_rvv_asm.S'),
        speexdsp_config_types_h,
        include_directories: checkasm_incs,
        dependencies: [checkasm_dep],
    )
    checkasm_extra_libs += static_library('checkasm_preproc_rvv',
        files('checkasm/preprocess/wrap_preproc_rvv.c',
              '../libspeexdsp/preprocess_rvv_asm.S'),
        speexdsp_config_types_h,
        include_directories: checkasm_incs,
        dependencies: [checkasm_dep],
    )
    checkasm_extra_libs += static_library('checkasm_smallft_rvv',
        files('checkasm/smallft/wrap_smallft_rvv.c',
              '../libspeexdsp/smallft_rvv_asm.S'),
        speexdsp_config_types_h,
        include_directories: checkasm_incs,
        dependencies: [checkasm_dep],
    )
    checkasm_extra_libs += static_library('checkasm_fbank_rvv',
        files('checkasm/filterbank/wrap_fbank_rvv.c',
              '../libspeexdsp/filterbank_rvv_asm.S'),
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
